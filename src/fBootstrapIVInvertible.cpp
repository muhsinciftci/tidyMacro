#include "fBootstrapIVInvertible.h"
#include "fGetBands.h"
#include "fOLS.h"
#include "fVAR.h"
#include "fgenerateVARdata.h"
#include "flagmakerMatrix.h"
#include "fwoldIRF.h"
#include <RcppArmadillo.h>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]

static arma::vec clean_instr_inv(const arma::vec& instr_raw, int p,
                                  const arma::mat& y) {
    int T = static_cast<int>(y.n_rows);
    arma::mat instr_mat(instr_raw.n_elem, 1);
    instr_mat.col(0) = instr_raw;
    arma::vec z0 = instr_raw.rows(p, T - 1);
    arma::mat z0_mat(z0.n_elem, 1);
    z0_mat.col(0) = z0;
    arma::mat lag_i   = flagmakerMatrix(instr_mat, p);
    arma::mat lag_X   = flagmakerMatrix(y, p);
    arma::mat X_clean = arma::join_horiz(lag_i, lag_X);
    return fOLS_cpp(z0_mat, X_clean, 1, 0, 0, 1).err.col(0);
}

arma::cube fBootstrapIVInvertible_cpp(
    const arma::mat& y,
    const arma::vec& instr,
    const VARResult& var_result,
    int nboot, int p, int c, int hor,
    const arma::ivec& cumu)
{
    const arma::mat& beta = var_result.beta;
    const arma::mat& eps  = var_result.residuals;

    const int N = static_cast<int>(y.n_cols);
    const int T = static_cast<int>(y.n_rows);
    const int m = static_cast<int>(eps.n_rows);
    const int H = hor + 1;

    arma::cube bootirf(N, H, nboot, arma::fill::zeros);

    // Pre-generate all Rademacher signs sequentially (Armadillo RNG is not thread-safe)
    arma::mat all_rad(m, nboot);
    for (int b = 0; b < nboot; ++b) {
        arma::vec u = arma::randu<arma::vec>(m);
        for (int i = 0; i < m; ++i) all_rad(i, b) = (u(i) > 0.5) ? 1.0 : -1.0;
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int b = 0; b < nboot; ++b) {

        // 1. Wild bootstrap: apply pre-generated Rademacher signs to VAR residuals
        const arma::vec rad = all_rad.col(b);

        arma::mat eps_wild(m, N, arma::fill::none);
        for (int i = 0; i < m; ++i) eps_wild.row(i) = eps.row(i) * rad(i);
        arma::mat y_new = fgenerateVARdata(y, p, c, beta, eps_wild);

        // 2. Apply same Rademacher signs to instrument (indices p..T-1)
        arma::vec z_new = instr;
        for (int i = p; i < T; ++i) z_new(i) = instr(i) * rad(i - p);

        // 3. Clean instrument
        arma::vec z_b = clean_instr_inv(z_new, p, y_new);

        // 4. Estimate VAR and Wold IRF
        VARResult     vr_b = fVAR_cpp(y_new, p, c, R_NilValue);
        WoldIRFResult wr_b = fwoldIRF_cpp(vr_b, hor);
        const arma::mat& eps_b    = vr_b.residuals;
        const arma::mat  Sigma_inv = arma::inv_sympd(vr_b.sigma);
        const arma::cube& wold_b  = wr_b.irfwold;

        // 5. Re-estimate structural impact vector via IV (lightweight OLS: beta only)
        arma::mat z_b_mat(z_b.n_elem, 1);
        z_b_mat.col(0) = z_b;
        arma::mat s0 = fOLS_cpp(arma::mat(eps_b), z_b_mat, 0, 0, 0, 1).beta; // 1 x N
        double s0_1 = s0(0, 0);
        arma::vec s_b = s0.t() / s0_1;
        double scaler_b = std::sqrt(arma::dot(s_b, Sigma_inv * s_b));

        // 6. Compute IRF, cumulate, sign-normalise
        arma::mat irf_b(N, H, arma::fill::none);
        for (int h = 0; h < H; ++h)
            irf_b.col(h) = wold_b.slice(h) * s_b / scaler_b;

        for (int i = 0; i < static_cast<int>(cumu.n_elem); ++i) {
            int v = cumu(i) - 1;
            for (int h = 1; h < H; ++h) irf_b(v, h) += irf_b(v, h - 1);
        }

        double sgn = (irf_b(0, 0) >= 0.0) ? 1.0 : -1.0;
        bootirf.slice(b) = irf_b * sgn;
    }

    return bootirf;
}

//' Bootstrap IRFs for the Invertible Case (Wild Bootstrap)
//'
//' Implements the wild bootstrap with Rademacher signs for standard SVAR-IV.
//' The same Rademacher vector is applied to both the VAR residuals and the
//' raw instrument; at each draw the structural impact vector is re-estimated
//' by IV regression.
//'
//' @param y T x N matrix of endogenous variables.
//' @param instr Length-T raw instrument vector (before cleaning).
//' @param var_result List from \code{fVAR}.
//' @param nboot Number of bootstrap replications.
//' @param p VAR lag order.
//' @param c Integer (0/1): include intercept.
//' @param hor Impulse-response horizon.
//' @param cumu 1-indexed integer vector of variable positions to cumulate.
//' @param prc Primary confidence level (e.g. 90). Default 90.
//' @param prc2 Secondary confidence level (e.g. 68). Default 68.
//'
//' @return List with N x (hor+1) matrices: upper, lower (at \code{prc}),
//'   upper2, lower2 (at \code{prc2}), median.
//'
//' @seealso \code{\link{fBootstrapIVRecover}}, \code{\link{fGetBands}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fBootstrapIVInvertible(
    const arma::mat& y,
    const arma::vec& instr,
    const Rcpp::List& var_result,
    int nboot, int p, int c, int hor,
    const arma::ivec& cumu,
    double prc = 90.0, double prc2 = 68.0)
{
    VARResult vr;
    vr.beta      = Rcpp::as<arma::mat>(var_result["beta"]);
    vr.residuals = Rcpp::as<arma::mat>(var_result["residuals"]);
    vr.sigma     = Rcpp::as<arma::mat>(var_result["sigma"]);
    vr.p         = Rcpp::as<int>(var_result["p"]);
    vr.c         = Rcpp::as<int>(var_result["c"]);
    vr.n_exog    = var_result.containsElementNamed("n_exog") ?
                   Rcpp::as<int>(var_result["n_exog"]) : 0;

    arma::cube bootirf = fBootstrapIVInvertible_cpp(
        y, instr, vr, nboot, p, c, hor, cumu);

    // Sort once, extract both confidence bands in a single pass
    FGetBands2Result bands = fGetBands2_cpp(bootirf, prc, prc2);

    return Rcpp::List::create(
        Rcpp::Named("upper")   = bands.upper,
        Rcpp::Named("lower")   = bands.lower,
        Rcpp::Named("upper2")  = bands.upper2,
        Rcpp::Named("lower2")  = bands.lower2,
        Rcpp::Named("median")  = bands.median
    );
}
