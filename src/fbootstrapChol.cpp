// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// CRITICAL: Define this BEFORE including RcppArmadillo to prevent conflicts
#define ARMA_DONT_USE_OPENMP

#include "fbootstrapChol.h"
#include "fVAR.h"
#include "fbootstrapVAR.h"
#include "fcholeskyIRF.h"
#include "fwoldIRF.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstdio>
#include <algorithm>
#include <vector>

// Internal C++ function (called from other C++ code)
BootstrapCholResult
fbootstrapChol_cpp(const arma::mat &y, const VARResult &var_result, int nboot,
                   int horizon, double prc, double prc2,
                   const std::string &bootscheme,
                   Rcpp::Nullable<arma::mat> exog, int n_threads = 0) {

  const int p      = var_result.p;
  const int c      = var_result.c;
  const int n_exog = var_result.n_exog;
  const int N      = y.n_cols;
  const int H      = horizon + 1;
  const int n_coef = var_result.beta.n_rows;
  const int slice_sz = N * N * H;

  if (n_exog > 0 && exog.isNull())
    Rcpp::stop("Original VAR used exogenous variables. You must provide the 'exog' parameter.");
  if (n_exog == 0 && exog.isNotNull())
    Rcpp::stop("Original VAR did not use exogenous variables. Do not provide the 'exog' parameter.");

  arma::mat bootchol_flat(slice_sz, nboot, arma::fill::zeros);
  arma::cube boot_beta(N, n_coef, nboot, arma::fill::zeros);

  // Convert exog once outside the parallel region (avoids repeated Rcpp::as).
  const bool has_exog = exog.isNotNull();
  arma::mat exog_mat;
  if (has_exog) exog_mat = Rcpp::as<arma::mat>(exog);

  int actual_threads = 1;
#ifdef _OPENMP
  actual_threads = (n_threads == 0)
                       ? std::max(1, omp_get_max_threads() - 1)
                       : n_threads;
  omp_set_num_threads(actual_threads);
  std::printf("Using %d thread(s) for parallel bootstrap computation...\n",
              actual_threads);
#else
  std::printf("OpenMP not available. Running in single-threaded mode.\n");
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int b = 0; b < nboot; ++b) {

    // 1. Bootstrap a new dataset
    BootstrapVARResult boot_data = fbootstrapVAR_cpp(y, var_result, bootscheme);

    // 2. Estimate VAR on bootstrapped data
    VARResult var_loop = has_exog
                           ? fVAR_cpp(boot_data.ynext, p, c, exog_mat)
                           : fVAR_cpp(boot_data.ynext, p, c, R_NilValue);
    boot_beta.slice(b) = var_loop.beta.t();

    // 3. Compute Wold IRFs
    WoldIRFResult wold_result = fwoldIRF_cpp(var_loop, horizon);
    const arma::cube &wold_loop = wold_result.irfwold;

    // 4. Lower Cholesky factor of sigma
    arma::mat S_loop = arma::chol(var_loop.sigma, "lower");

    // 5. Cholesky IRFs
    arma::cube cholirf_loop = fcholeskyIRF(wold_loop, S_loop);

    bootchol_flat.col(b) = arma::vectorise(cholirf_loop);
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
    arma::rowvec row_rv = bootchol_flat.row(i);
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

  BootstrapCholResult result;
  result.bootchol_flat = bootchol_flat;
  result.upper         = upper;
  result.lower         = lower;
  result.upper2        = upper2;
  result.lower2        = lower2;
  result.boot_beta     = boot_beta;
  result.N             = N;
  result.H             = H;
  return result;
}

//' Bootstrap Cholesky Identified Impulse Response Functions
//'
//' @param y T x N matrix of original endogenous variables.
//' @param var_result List from \code{fVAR()}.
//' @param nboot Number of bootstrap replications.
//' @param horizon Maximum IRF horizon.
//' @param prc Outer confidence level in percent (e.g. 90).
//' @param prc2 Inner confidence level in percent (e.g. 68).
//' @param bootscheme \code{"residual"} or \code{"wild"}.
//' @param exog Optional T x M matrix of exogenous variables. Must be provided
//'   if the original VAR was estimated with exogenous variables; must be
//'   omitted otherwise.
//' @param n_threads Number of OpenMP threads. 0 = all cores minus one.
//'   Default 0.
//'
//' @return List with elements:
//'   \describe{
//'     \item{bootchol}{N x N x (horizon+1) x nboot array of bootstrapped Cholesky IRFs}
//'     \item{upper}{N x N x (horizon+1) outer upper confidence bands}
//'     \item{lower}{N x N x (horizon+1) outer lower confidence bands}
//'     \item{upper2}{N x N x (horizon+1) inner upper confidence bands}
//'     \item{lower2}{N x N x (horizon+1) inner lower confidence bands}
//'     \item{boot_beta}{N x (Np+c) x nboot array of bootstrapped coefficients}
//'   }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fbootstrapChol(const arma::mat &y, const Rcpp::List &var_result,
                          int nboot, int horizon, double prc = 90.0,
                          double prc2 = 68.0,
                          const std::string &bootscheme = "residual",
                          Rcpp::Nullable<arma::mat> exog = R_NilValue,
                          int n_threads = 0) {

  VARResult var_result_struct;
  var_result_struct.beta      = Rcpp::as<arma::mat>(var_result["beta"]);
  var_result_struct.residuals = Rcpp::as<arma::mat>(var_result["residuals"]);
  var_result_struct.sigma     = Rcpp::as<arma::mat>(var_result["sigma"]);
  var_result_struct.p         = Rcpp::as<int>(var_result["p"]);
  var_result_struct.c         = Rcpp::as<int>(var_result["c"]);
  var_result_struct.n_exog    = var_result.containsElementNamed("n_exog")
                                  ? Rcpp::as<int>(var_result["n_exog"])
                                  : 0;

  BootstrapCholResult result = fbootstrapChol_cpp(y, var_result_struct, nboot,
                                                  horizon, prc, prc2,
                                                  bootscheme, exog, n_threads);

  Rcpp::NumericVector bootchol_out(result.bootchol_flat.begin(),
                                   result.bootchol_flat.end());
  bootchol_out.attr("dim") = Rcpp::IntegerVector::create(result.N, result.N,
                                                          result.H, nboot);

  return Rcpp::List::create(Rcpp::Named("bootchol")   = bootchol_out,
                            Rcpp::Named("upper")      = result.upper,
                            Rcpp::Named("lower")      = result.lower,
                            Rcpp::Named("upper2")     = result.upper2,
                            Rcpp::Named("lower2")     = result.lower2,
                            Rcpp::Named("boot_beta")  = result.boot_beta);
}
