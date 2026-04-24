// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// Prevent Armadillo from spawning its own OpenMP threads inside our parallel region
#define ARMA_DONT_USE_OPENMP

#include "fhdIV.h"
#include "fGetShock.h"
#include "fcompanionMatrix.h"
#include "fOLS.h"
#include "fVAR.h"
#include "fgenerateVARdata.h"
#include "fmbb_var.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// =========================================================================
// Internal helper: companion-form HD recursion
//
// Given unit-normalised structural impact vector s and VAR parameters,
// recovers the structural shock via Stock-Watson (2018) and propagates it
// through the companion form to obtain the T_eff x N HD matrix.
//
// HD_big[:, t] = F * HD_big[:, t-1]  +  (s * V[t-1]) padded with zeros
// HDshock      = HD_big[0:N-1, 1:T_eff].t()   →  T_eff x N
// =========================================================================
static arma::mat hd_recursion_iv(const arma::mat& residuals,
                                  const arma::mat& sigma,
                                  const arma::vec& s,
                                  const arma::mat& beta,
                                  int c, int p) {

    const int N       = static_cast<int>(residuals.n_cols);
    const int T_eff   = static_cast<int>(residuals.n_rows);
    const int nvarXeq = N * p;

    // Structural shock series (unit normalisation, shockSize = 1)
    arma::vec V = fGetShock_cpp(residuals, sigma, s, 1.0, "unit");

    // Companion matrix F  (nvarXeq x nvarXeq)
    CompanionMatrixResult comp = fcompanionMatrix_cpp(beta, c, p);
    const arma::mat& F = comp.comp;

    // Companion-form recursion
    arma::mat HD_big(nvarXeq, T_eff + 1, arma::fill::zeros);

    for (int i = 1; i <= T_eff; ++i) {
        HD_big.col(i)         = F * HD_big.col(i - 1);
        HD_big.col(i).head(N) += s * V(i - 1);
    }

    // Return T_eff x N  (drop the t=0 initialisation column, keep top N rows)
    return HD_big.head_rows(N).cols(1, T_eff).t();
}


// =========================================================================
// fhdIV_cpp  —  point-estimate companion-form HD for IV-identified shock
// =========================================================================
HDIVResult fhdIV_cpp(const arma::mat& residuals,
                     const arma::mat& sigma,
                     const arma::vec& s,
                     const arma::mat& beta,
                     int c, int p) {
    HDIVResult result;
    result.HDshock = hd_recursion_iv(residuals, sigma, s, beta, c, p);
    return result;
}


// =========================================================================
// fbootstrapHDIV_cpp  —  MBB bootstrap uncertainty bands for the IV HD
//
// For each bootstrap draw:
//   1. Resample residuals and instrument via MBB
//   2. Generate and re-estimate bootstrap VAR
//   3. Re-identify the IV shock (first-stage OLS, second-stage projection)
//   4. Run companion-form HD recursion
//
// Quantiles are recentered around the point estimate following:
//   Kaenzig (2021) / Gertler-Karadi MATLAB convention:
//   upper = q(1-a/2) - median(boot) + HDshock_point
// =========================================================================
BootstrapHDIVResult fbootstrapHDIV_cpp(const arma::mat& y,
                                        const VARResult& var_result,
                                        const arma::mat& Z,
                                        const arma::vec& s,
                                        int nboot, int blocksize,
                                        const arma::ivec& adjustZ,
                                        const arma::ivec& adjustu,
                                        int policyvar,
                                        double prc,
                                        int n_threads) {

    const int N      = static_cast<int>(y.n_cols);
    const int T_eff  = static_cast<int>(var_result.residuals.n_rows);
    const int p      = var_result.p;
    const int c      = var_result.c;
    const int pv_idx = policyvar - 1;   // 0-based

    // Convert 1-based R indices to 0-based C++ indices
    const int Z_start = adjustZ(0) - 1;
    const int Z_end   = adjustZ(1) - 1;
    const int u_start = adjustu(0) - 1;
    const int u_end   = adjustu(1) - 1;

    const arma::mat Zstar   = Z.rows(Z_start, Z_end);
    const arma::mat Epsstar = var_result.residuals.rows(u_start, u_end);

    // Point-estimate HD (used for recentering)
    arma::mat HDshock_pt = hd_recursion_iv(var_result.residuals,
                                            var_result.sigma,
                                            s, var_result.beta, c, p);

    // Storage for bootstrap HDs: T_eff x N x nboot
    arma::cube bootHD(T_eff, N, nboot, arma::fill::none);

    // OpenMP setup
    int actual_threads = 1;
#ifdef _OPENMP
    actual_threads = (n_threads == 0)
        ? std::max(1, omp_get_max_threads() - 1)
        : n_threads;
    omp_set_num_threads(actual_threads);
    Rprintf("Using %d thread(s) for bootstrap HD computation...\n", actual_threads);
#else
    Rprintf("OpenMP not available. Running in single-threaded mode.\n");
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int b = 0; b < nboot; ++b) {
        // 1. MBB resample
        MBBVARResult mbb = fmbb_var_cpp(Epsstar, p, blocksize, Zstar);

        // 2. Pad residuals outside the proxy window with originals
        arma::mat res_boot = var_result.residuals;
        res_boot.rows(u_start, u_end) = mbb.eps_boot;

        // 3. Generate bootstrap VAR data and re-estimate
        arma::mat varboot  = fgenerateVARdata(y, p, c, var_result.beta, res_boot);
        VARResult var_loop = fVAR_cpp(varboot, p, c, R_NilValue);

        // 4. IV identification on bootstrap residuals
        arma::vec u_p_b = var_loop.residuals.col(pv_idx).rows(u_start, u_end);
        arma::mat u_q_b = var_loop.residuals.rows(u_start, u_end);
        u_q_b.shed_col(pv_idx);

        OLSResult ols_b     = fOLS_cpp(u_p_b, mbb.M_boot, 0, 0, 0, 0);
        const arma::mat& uhat_b = ols_b.fitted_partial;
        arma::mat sq_sp_b   = arma::solve(uhat_b.t() * uhat_b,
                                          uhat_b.t() * u_q_b);

        arma::vec s_b(N, arma::fill::zeros);
        int counter = 0;
        for (int i = 0; i < N; ++i) {
            s_b(i) = (i == pv_idx) ? 1.0 : sq_sp_b(counter++);
        }

        // 5. Bootstrap HD
        bootHD.slice(b) = hd_recursion_iv(var_loop.residuals, var_loop.sigma,
                                           s_b, var_loop.beta, c, p);
    }

    // Compute percentiles and recenter (Kaenzig/Gertler-Karadi convention)
    const double up_pct  = 50.0 + prc * 0.5;
    const double low_pct = 50.0 - prc * 0.5;

    arma::mat upper(T_eff, N, arma::fill::none);
    arma::mat lower(T_eff, N, arma::fill::none);
    arma::mat boot_med(T_eff, N, arma::fill::none);

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int t = 0; t < T_eff; ++t) {
        for (int n = 0; n < N; ++n) {
            // Extract bootstrap values for this (t, n) without arma::tube view
            arma::vec vals(nboot);
            for (int b = 0; b < nboot; ++b) {
                vals(b) = bootHD(t, n, b);
            }
            arma::vec sorted = arma::sort(vals);

            // Linear interpolation at requested percentiles
            auto interp = [&](double pct) -> double {
                double idx = (pct / 100.0) * (nboot - 1);
                int    fl  = static_cast<int>(idx);
                double fr  = idx - fl;
                if (fl < nboot - 1)
                    return sorted(fl) * (1.0 - fr) + sorted(fl + 1) * fr;
                return sorted(nboot - 1);
            };

            boot_med(t, n) = interp(50.0);
            upper(t, n)    = interp(up_pct)  - boot_med(t, n) + HDshock_pt(t, n);
            lower(t, n)    = interp(low_pct) - boot_med(t, n) + HDshock_pt(t, n);
        }
    }

    BootstrapHDIVResult result;
    result.HDshock = HDshock_pt;
    result.upper   = upper;
    result.lower   = lower;
    return result;
}


// =========================================================================
// R wrappers
// =========================================================================

//' Companion-Form Historical Decomposition for IV-Identified Shock
//'
//' Computes the companion-form historical decomposition attributing each
//' observation of all VAR variables to the IV-identified structural shock,
//' following Stock and Watson (2018).
//'
//' @param residuals (T-p) x N matrix of VAR residuals from \code{fVAR}.
//' @param sigma N x N residual covariance matrix (e.g. \code{var_result$sigma}).
//' @param s N x 1 unit-normalised structural impact vector (first element = 1).
//' @param beta (Np+c) x N VAR coefficient matrix from \code{fVAR}.
//' @param c Integer intercept indicator (1 = include, 0 = exclude).
//' @param p Integer VAR lag order.
//'
//' @return A list with element:
//'   \item{HDshock}{(T-p) x N matrix. Column \code{j} is the contribution
//'     of the IV shock to variable \code{j} at each point in time.}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fhdIV(const arma::mat& residuals,
                 const arma::mat& sigma,
                 const arma::vec& s,
                 const arma::mat& beta,
                 int c, int p) {

    HDIVResult result = fhdIV_cpp(residuals, sigma, s, beta, c, p);
    return Rcpp::List::create(Rcpp::Named("HDshock") = result.HDshock);
}


//' Bootstrap Confidence Bands for the IV Historical Decomposition
//'
//' Computes MBB bootstrap uncertainty bands for the companion-form historical
//' decomposition of an IV-identified structural shock. In each replication
//' the VAR is re-estimated, the IV identification is re-run, and the
//' companion-form HD recursion is applied to the bootstrap draw. Bands are
//' recentered around the point estimate following Kaenzig (2021).
//'
//' @param y T x N matrix of original endogenous variables.
//' @param var_result List from \code{fVAR}.
//' @param Z Instrument matrix aligned to the proxy sample.
//' @param s N x 1 unit-normalised point-estimate structural impact vector.
//' @param nboot Number of bootstrap replications.
//' @param blocksize MBB block size.
//' @param adjustZ Integer vector \code{c(start, end)} selecting the proxy-sample
//'   rows of \code{Z} (1-based).
//' @param adjustu Integer vector \code{c(start, end)} selecting the proxy-sample
//'   rows of the residuals (1-based).
//' @param policyvar Integer (1-based) index of the IV policy variable. Default 1.
//' @param prc Confidence level in percent (e.g. 90 for 90\% CI). Default 90.
//' @param n_threads OpenMP threads. 0 = all cores minus one. Default 0.
//'
//' @return A list with elements:
//'   \item{HDshock}{(T-p) x N point-estimate HD matrix.}
//'   \item{upper}{(T-p) x N upper confidence bands.}
//'   \item{lower}{(T-p) x N lower confidence bands.}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fbootstrapHDIV(const arma::mat& y,
                           const Rcpp::List& var_result,
                           const arma::mat& Z,
                           const arma::vec& s,
                           int nboot, int blocksize,
                           const arma::ivec& adjustZ,
                           const arma::ivec& adjustu,
                           int policyvar = 1,
                           double prc = 90.0,
                           int n_threads = 0) {

    VARResult var_struct;
    var_struct.beta      = Rcpp::as<arma::mat>(var_result["beta"]);
    var_struct.residuals = Rcpp::as<arma::mat>(var_result["residuals"]);
    var_struct.sigma     = Rcpp::as<arma::mat>(var_result["sigma"]);
    var_struct.p         = Rcpp::as<int>(var_result["p"]);
    var_struct.c         = Rcpp::as<int>(var_result["c"]);
    var_struct.n_exog    = var_result.containsElementNamed("n_exog")
                               ? Rcpp::as<int>(var_result["n_exog"]) : 0;

    BootstrapHDIVResult result = fbootstrapHDIV_cpp(
        y, var_struct, Z, s, nboot, blocksize, adjustZ, adjustu,
        policyvar, prc, n_threads);

    return Rcpp::List::create(
        Rcpp::Named("HDshock") = result.HDshock,
        Rcpp::Named("upper")   = result.upper,
        Rcpp::Named("lower")   = result.lower);
}
