// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// CRITICAL: Define this BEFORE including RcppArmadillo to prevent conflicts
#define ARMA_DONT_USE_OPENMP

#include "fbootstrapBQ.h"
#include "fVAR.h"
#include "fbootstrapVAR.h"
#include "fbqIRF.h"
#include "fwoldIRF.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstdio>
#include <algorithm>
#include <vector>

// Internal C++ function (called from other C++ code)
// Calls _cpp versions of internal functions - no list overhead!
BootstrapBQResult
fbootstrapBQ_cpp(const arma::mat &y, const VARResult &var_result, int nboot,
                 int horizon, double prc, double prc2, const std::string &bootscheme,
                 const arma::uvec &cumulate,
                 Rcpp::Nullable<arma::vec> scaling, int n_threads = 0) {

  const int p      = var_result.p;
  const int c      = var_result.c;
  const int T      = y.n_rows;
  const int N      = y.n_cols;
  const int H      = horizon + 1;
  const int n_coef = var_result.beta.n_rows;
  const int slice_sz = N * N * H;

  const double df = static_cast<double>(T - 1 - p - N * p);
  if (df <= 0.0)
    Rcpp::stop("Degrees of freedom T-1-p-N*p = %.0f <= 0.", df);

  arma::mat bootbq_flat(slice_sz, nboot, arma::fill::zeros);
  arma::cube boot_beta(N, n_coef, nboot, arma::fill::zeros);

  // Hoist the Rcpp::Nullable<arma::vec> scaling parse out of the parallel loop.
  bool   has_scaling = scaling.isNotNull();
  int    scale_idx   = 0;
  double scale_size  = 0.0;
  if (has_scaling) {
    arma::vec s = Rcpp::as<arma::vec>(scaling);
    scale_idx   = static_cast<int>(s(0)) - 1;
    scale_size  = s(1);
  }

  int actual_threads = 1;

#ifdef _OPENMP
  if (n_threads == 0) {
    int max_threads = omp_get_max_threads();
    actual_threads = std::max(1, max_threads - 1);
  } else {
    actual_threads = n_threads;
  }
  omp_set_num_threads(actual_threads);
  std::printf("Using %d thread(s) for parallel bootstrap computation...\n",
              actual_threads);
#else
  std::printf("OpenMP not available. Running in single-threaded mode.\n");
  actual_threads = 1;
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int b = 0; b < nboot; ++b) {

    // Bootstrap a new dataset
    BootstrapVARResult boot_data = fbootstrapVAR_cpp(y, var_result, bootscheme);
    const arma::mat &varboot = boot_data.ynext;

    // Estimate VAR on bootstrapped data
    VARResult var_loop = fVAR_cpp(varboot, p, c, R_NilValue);
    boot_beta.slice(b) = var_loop.beta.t();

    // Compute Wold IRFs
    WoldIRFResult wold_result = fwoldIRF_cpp(var_loop, horizon);
    const arma::cube &wold_loop = wold_result.irfwold;

    // BQ identification: C1 = sum of Wold slices, D1 = lower chol, K = C1\D1
    // sigma from fVAR_cpp uses DF-corrected denominator: (T-p) - c - N*p - n_exog
    arma::mat C1(N, N, arma::fill::zeros);
    for (int h = 0; h < H; ++h) {
      C1 += wold_loop.slice(h);
    }

    arma::mat D1;
    arma::mat LR_cov = C1 * var_loop.sigma * C1.t();
    if (!arma::chol(D1, LR_cov, "lower")) {
      // Regularise if not positive definite
      LR_cov = 0.5 * (LR_cov + LR_cov.t());
      LR_cov.diag() += 1e-8 * arma::trace(LR_cov) / N;
      arma::chol(D1, LR_cov, "lower");
    }

    arma::mat K_loop = arma::solve(C1, D1);

    // Compute BQ IRFs with optional scaling (pre-parsed — no Rcpp::as here)
    arma::cube bqirf_loop =
        fbqIRF_cpp(wold_loop, K_loop, has_scaling, scale_idx, scale_size);

    // Cumulate selected variables along horizon (0-based indices in cumulate)
    // Extract to temporary vec first — arma::cumsum does not accept subviews
    for (arma::uword ci = 0; ci < cumulate.n_elem; ++ci) {
      arma::uword row_idx = cumulate(ci);
      for (int i = 0; i < N; ++i) {
        arma::vec tmp = bqirf_loop.tube(row_idx, i);
        bqirf_loop.tube(row_idx, i) = arma::cumsum(tmp);
      }
    }

    bootbq_flat.col(b) = arma::vectorise(bqirf_loop);
  }

  // Compute percentile bands
  const double up_pct   = 50.0 + prc  * 0.5;
  const double low_pct  = 50.0 - prc  * 0.5;
  const double up_pct2  = 50.0 + prc2 * 0.5;
  const double low_pct2 = 50.0 - prc2 * 0.5;

  arma::vec upper_vec(slice_sz);
  arma::vec lower_vec(slice_sz);
  arma::vec upper2_vec(slice_sz);
  arma::vec lower2_vec(slice_sz);

  auto nth_pct = [](std::vector<double>& v, double pct) -> double {
    const int n = static_cast<int>(v.size());
    const double raw = (pct / 100.0) * (n - 1);
    const int lo = static_cast<int>(raw);
    const double frac = raw - lo;
    std::nth_element(v.begin(), v.begin() + lo, v.end());
    const double lo_val = v[lo];
    if (frac < 1e-12 || lo + 1 >= n) return lo_val;
    return lo_val * (1.0 - frac) +
           *std::min_element(v.begin() + lo + 1, v.end()) * frac;
  };

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < slice_sz; ++i) {
    arma::rowvec row_rv = bootbq_flat.row(i);
    std::vector<double> sv(row_rv.begin(), row_rv.end());
    upper_vec(i)  = nth_pct(sv, up_pct);
    lower_vec(i)  = nth_pct(sv, low_pct);
    upper2_vec(i) = nth_pct(sv, up_pct2);
    lower2_vec(i) = nth_pct(sv, low_pct2);
  }

  arma::cube upper (upper_vec.memptr(),  N, N, H, false, true);
  arma::cube lower (lower_vec.memptr(),  N, N, H, false, true);
  arma::cube upper2(upper2_vec.memptr(), N, N, H, false, true);
  arma::cube lower2(lower2_vec.memptr(), N, N, H, false, true);

  BootstrapBQResult result;
  result.bootbq_flat = bootbq_flat;
  result.upper       = upper;
  result.lower       = lower;
  result.upper2      = upper2;
  result.lower2      = lower2;
  result.boot_beta   = boot_beta;
  result.N           = N;
  result.H           = H;
  return result;
}

//' Bootstrap Blanchard-Quah Long-Run Identified Impulse Response Functions
//'
//' @param y T x N matrix of original endogenous variables.
//' @param var_result List from \code{fVAR()}.
//' @param nboot Number of bootstrap replications.
//' @param horizon Maximum IRF horizon.
//' @param prc Confidence level in percent (e.g. 68).
//' @param bootscheme \code{"residual"} or \code{"wild"}.
//' @param cumulate Integer vector (1-based) of variable indices whose IRFs
//'   should be cumulated along the horizon. Typically used when the VAR is
//'   estimated in first differences and level responses are required.
//' @param scaling Optional numeric vector of length 2. First element is the
//'   variable index (1-based) for normalisation; second is the shock size.
//'   Default NULL (no normalisation).
//' @param n_threads Number of OpenMP threads. 0 = all cores minus one.
//'   Default 0.
//'
//' @return List with elements:
//'   \describe{
//'     \item{bootbq}{N x N x (horizon+1) x nboot array of bootstrapped BQ IRFs}
//'     \item{upper}{N x N x (horizon+1) upper confidence bands}
//'     \item{lower}{N x N x (horizon+1) lower confidence bands}
//'     \item{boot_beta}{N x (Np+c) x nboot array of bootstrapped coefficients}
//'   }
//'
//' @details
//' Implements the residual or wild bootstrap for BQ long-run identified VARs.
//' In each replication the full BQ identification is re-computed:
//' \eqn{C_1 = \sum_h \Psi_h}, \eqn{D_1 = \text{chol}(C_1 \Sigma C_1')^\top},
//' \eqn{K = C_1^{-1} D_1}. Selected variables are then cumulated if the VAR
//' is estimated in differences.
//'
//' @references
//' Blanchard, O. J., & Quah, D. (1989). The dynamic effects of aggregate
//' demand and supply disturbances. \emph{American Economic Review}, 79(4),
//' 655--673.
//'
//' Gali, J. (1999). Technology, employment, and the business cycle: Do
//' technology shocks explain aggregate fluctuations?
//' \emph{American Economic Review}, 89(1), 249--271.
//'
//' @examples
//' \dontrun{
//' var_result <- fVAR(y, p = 2, c = 1)
//' wold       <- fwoldIRF(var_result, horizon = 40)
//' Sigma      <- var_result$sigma
//' C1         <- apply(wold, c(1, 2), sum)
//' D1         <- t(chol(C1 %*% Sigma %*% t(C1)))
//' K          <- solve(C1, D1)
//' point_irf  <- fbqIRF(wold, K)
//'
//' # Bootstrap with cumulation of both variables (VAR in first differences)
//' boot <- fbootstrapBQ(y, var_result, nboot = 1000, horizon = 40,
//'                      prc = 68, bootscheme = "residual",
//'                      cumulate = c(1L, 2L))
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fbootstrapBQ(const arma::mat &y, const Rcpp::List &var_result,
                        int nboot, int horizon, double prc = 90.0, double prc2 = 68.0,
                        const std::string &bootscheme = "residual",
                        Rcpp::IntegerVector cumulate = Rcpp::IntegerVector(),
                        Rcpp::Nullable<arma::vec> scaling = R_NilValue,
                        int n_threads = 0) {

  // Unpack R list into VARResult struct
  VARResult var_result_struct;
  var_result_struct.beta       = Rcpp::as<arma::mat>(var_result["beta"]);
  var_result_struct.residuals  = Rcpp::as<arma::mat>(var_result["residuals"]);
  var_result_struct.sigma = Rcpp::as<arma::mat>(var_result["sigma"]);
  var_result_struct.p          = Rcpp::as<int>(var_result["p"]);
  var_result_struct.c          = Rcpp::as<int>(var_result["c"]);
  var_result_struct.n_exog     = var_result.containsElementNamed("n_exog")
                                   ? Rcpp::as<int>(var_result["n_exog"])
                                   : 0;

  // Convert cumulate from 1-based (R) to 0-based (C++)
  arma::uvec cumulate_cpp = (cumulate.size() == 0)
      ? arma::uvec()
      : Rcpp::as<arma::uvec>(cumulate) - 1;

  BootstrapBQResult result = fbootstrapBQ_cpp(y, var_result_struct, nboot,
                                              horizon, prc, prc2, bootscheme,
                                              cumulate_cpp, scaling, n_threads);

  Rcpp::NumericVector bootbq_out(result.bootbq_flat.begin(),
                                 result.bootbq_flat.end());
  bootbq_out.attr("dim") = Rcpp::IntegerVector::create(result.N, result.N,
                                                        result.H, nboot);

  return Rcpp::List::create(Rcpp::Named("bootbq")    = bootbq_out,
                            Rcpp::Named("upper")     = result.upper,
                            Rcpp::Named("lower")     = result.lower,
                            Rcpp::Named("upper2")    = result.upper2,
                            Rcpp::Named("lower2")    = result.lower2,
                            Rcpp::Named("boot_beta") = result.boot_beta);
}
