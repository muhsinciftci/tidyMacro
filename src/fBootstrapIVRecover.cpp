#include "fBootstrapIVRecover.h"
#include "fGetBands.h"
#include "fOLS.h"
#include "fVAR.h"
#include "fgenerateVARdata.h"
#include "flagmakerMatrix.h"
#include "fwoldIRF.h"
#include "fPolyConvolve.h"
#include <RcppArmadillo.h>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]

static arma::vec clean_instr_rec(const arma::vec& instr_raw, int p,
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

arma::cube fBootstrapIVRecover_cpp(
    const arma::mat& y,
    const arma::vec& instr,
    const VARResult& var_result,
    const arma::vec& noise,
    const arma::vec& delta,
    int nboot, int p, int c, int r, int hor,
    const arma::ivec& cumu)
{
    const arma::mat& beta = var_result.beta;
    const arma::mat& eps  = var_result.residuals;

    const int N       = static_cast<int>(y.n_cols);
    const int T       = static_cast<int>(y.n_rows);
    const int m       = static_cast<int>(eps.n_rows);
    const int t_shock = m - r;
    const int H       = hor + 1;

    arma::cube bootirf(N, H, nboot, arma::fill::zeros);

    // Pre-generate all bootstrap indices sequentially (Armadillo RNG is not thread-safe)
    arma::imat all_idx(t_shock, nboot);
    for (int b = 0; b < nboot; ++b)
        all_idx.col(b) = arma::randi<arma::ivec>(
            t_shock, arma::distr_param(0, t_shock - 1));

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int b = 0; b < nboot; ++b) {

        // 1. iid joint resample of VAR residuals and noise
        const arma::ivec idx = all_idx.col(b);

        arma::mat eps_boot(t_shock, N, arma::fill::none);
        arma::vec v_boot(t_shock, arma::fill::none);
        for (int i = 0; i < t_shock; ++i) {
            eps_boot.row(i) = eps.row(idx(i));
            v_boot(i)       = noise(idx(i));
        }
        // Append final r rows of original eps (fixed final conditions)
        arma::mat eps_new = arma::join_vert(eps_boot, eps.rows(m - r, m - 1));

        // 2. Generate new y
        arma::mat y_new = fgenerateVARdata(y, p, c, beta, eps_new);

        // 3. Regenerate instrument from DGP: z(t) = delta(F)*e(t) + v(t)
        arma::mat eps_f(t_shock, N * (r + 1), arma::fill::none);
        for (int j = 0; j <= r; ++j)
            eps_f.cols(j * N, (j + 1) * N - 1) =
                eps_new.rows(j, j + t_shock - 1);
        arma::mat X_boot = arma::join_horiz(arma::ones(t_shock, 1), eps_f);

        arma::vec z_new = instr;
        z_new.rows(p, T - r - 1) = X_boot * delta + v_boot;

        // 4. Clean instrument
        arma::vec z_b = clean_instr_rec(z_new, p, y_new);
        const int m_b = static_cast<int>(z_b.n_elem);

        // 5. Estimate VAR and Wold IRF
        VARResult     vr_b = fVAR_cpp(y_new, p, c, R_NilValue);
        WoldIRFResult wr_b = fwoldIRF_cpp(vr_b, hor);
        const arma::mat& eps_b    = vr_b.residuals;
        const arma::mat  Sigma_inv = arma::inv_sympd(vr_b.sigma);

        // 6. Re-estimate psi(L): regress eps on [1, z_b, lags(z_b)]
        arma::mat z_b_mat(z_b.n_elem, 1);
        z_b_mat.col(0) = z_b;
        arma::mat lag_zb = flagmakerMatrix(z_b_mat, r);
        arma::mat Z_b = arma::join_horiz(
            arma::join_horiz(arma::ones(m_b - r, 1), z_b.rows(r, m_b - 1)),
            lag_zb);

        // Lightweight OLS: beta only (mode=1 skips r2/varbhat/F-stat)
        arma::mat psi_b = fOLS_cpp(
            arma::mat(eps_b.rows(r, m_b - 1)), Z_b, 0, 0, 0, 1).beta;

        arma::cube psiL_b(N, 1, r + 1, arma::fill::zeros);
        double a2_b = 0.0;
        for (int i = 0; i <= r; ++i) {
            arma::vec psi_temp = psi_b.row(1 + i).t();
            psiL_b.slice(i).col(0) = psi_temp;
            a2_b += arma::dot(psi_temp, Sigma_inv * psi_temp);
        }

        // 7. Compute IRF, cumulate, sign-normalise
        arma::cube irf_raw =
            fPolyConvolve_cpp(wr_b.irfwold, psiL_b, H) / std::sqrt(a2_b);

        arma::mat irf_mb(N, H, arma::fill::none);
        for (int h = 0; h < H; ++h) irf_mb.col(h) = irf_raw.slice(h).col(0);

        for (int i = 0; i < static_cast<int>(cumu.n_elem); ++i) {
            int v = cumu(i) - 1;
            for (int h = 1; h < H; ++h) irf_mb(v, h) += irf_mb(v, h - 1);
        }

        double sgn = (irf_mb(0, 0) >= 0.0) ? 1.0 : -1.0;
        bootirf.slice(b) = irf_mb * sgn;
    }

    return bootirf;
}

//' Bootstrap IRFs for the Recoverable Non-Invertible Case
//'
//' Implements the iid joint-resampling bootstrap for generalized SVAR-IV when
//' the shock is non-invertible but recoverable (Forni, Gambetti & Ricco 2022).
//' At each draw the VAR residuals and the instrument noise are resampled jointly,
//' the instrument is regenerated from the DGP
//' \eqn{z(t) = \delta(F)e(t) + v(t)}, and \eqn{\Psi(L)} is re-estimated.
//'
//' @param y T x N matrix of endogenous variables.
//' @param instr Length-T raw instrument vector (before cleaning).
//' @param var_result List from \code{fVAR}.
//' @param noise Length-(T-p-r) noise residuals \eqn{v_t} from the
//'   invertibility regression.
//' @param delta Length-\eqn{N(r+1)+1} coefficient vector (intercept + current
//'   and r future VAR residual coefficients) from the invertibility regression.
//' @param nboot Number of bootstrap replications.
//' @param p VAR lag order.
//' @param c Integer (0/1): include intercept.
//' @param r Number of future residuals used in the invertibility test.
//' @param hor Impulse-response horizon.
//' @param cumu 1-indexed integer vector of variable positions to cumulate.
//' @param prc Primary confidence level (e.g. 90). Default 90.
//' @param prc2 Secondary confidence level (e.g. 68). Default 68.
//'
//' @return List with N x (hor+1) matrices: upper, lower (at \code{prc}),
//'   upper2, lower2 (at \code{prc2}), median.
//'
//' @seealso \code{\link{fBootstrapIVInvertible}}, \code{\link{fGetBands}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fBootstrapIVRecover(
    const arma::mat& y,
    const arma::vec& instr,
    const Rcpp::List& var_result,
    const arma::vec& noise,
    const arma::vec& delta,
    int nboot, int p, int c, int r, int hor,
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

    arma::cube bootirf = fBootstrapIVRecover_cpp(
        y, instr, vr, noise, delta, nboot, p, c, r, hor, cumu);

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
