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

// Internal C++ function (called from other C++ code)
// Calls _cpp versions of internal functions - no list overhead!
BootstrapCholResult
fbootstrapChol_cpp(const arma::mat &y, const VARResult &var_result, int nboot,
                   int horizon, double prc, const std::string &bootscheme,
                   Rcpp::Nullable<arma::mat> exog, int n_threads = 0) {

  const arma::mat &beta = var_result.beta;
  const int p = var_result.p;
  const int c = var_result.c;
  const int n_exog = var_result.n_exog;

  if (n_exog > 0 && exog.isNull()) {
    Rcpp::stop("Original VAR estimation used exogenous variables. You must "
               "provide the 'exog' parameter.");
  }

  if (n_exog == 0 && exog.isNotNull()) {
    Rcpp::stop("Original VAR estimation did not use exogenous variables. Do "
               "not provide the 'exog' parameter.");
  }

  const int T = y.n_rows;
  const int N = y.n_cols;
  const int H = horizon + 1;
  const int n_coef = beta.n_rows;
  const int slice_sz = N * N * H;

  const double df = static_cast<double>(T - 1 - p - N * p);
  if (df <= 0.0)
    Rcpp::stop("Degrees of freedom T-1-p-N*p = %.0f <= 0.", df);

  arma::mat bootchol_flat(slice_sz, nboot, arma::fill::zeros);
  arma::cube boot_beta(N, n_coef, nboot, arma::fill::zeros);

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

    BootstrapVARResult boot_data = fbootstrapVAR_cpp(y, var_result, bootscheme);
    const arma::mat &varboot = boot_data.ynext;

    VARResult var_loop;
    if (exog.isNotNull()) {
      var_loop = fVAR_cpp(varboot, p, c, exog);
    } else {
      var_loop = fVAR_cpp(varboot, p, c, R_NilValue);
    }

    const arma::mat &betaloop = var_loop.beta;
    const arma::mat &err_loop = var_loop.residuals;

    boot_beta.slice(b) = betaloop.t();

    arma::mat sigma_loop = (err_loop.t() * err_loop) / df;

    arma::mat S_loop;
    if (!arma::chol(S_loop, sigma_loop, "lower")) {
      sigma_loop = 0.5 * (sigma_loop + sigma_loop.t());
      sigma_loop.diag() += 1e-8 * arma::trace(sigma_loop) / N;
      arma::chol(S_loop, sigma_loop, "lower");
    }

    WoldIRFResult wold_result = fwoldIRF_cpp(var_loop, horizon);
    const arma::cube &wold_loop = wold_result.irfwold;

    arma::cube cholirf_loop = fcholeskyIRF(wold_loop, S_loop);

    bootchol_flat.col(b) = arma::vectorise(cholirf_loop);
  }

  const double up_pct = 50.0 + prc * 0.5;
  const double low_pct = 50.0 - prc * 0.5;

  arma::vec upper_vec(slice_sz);
  arma::vec lower_vec(slice_sz);

  auto pctile = [](const arma::vec &sv, double pct) -> double {
    const int n = sv.n_elem;
    double idx = (pct / 100.0) * (n - 1);
    int lo = static_cast<int>(idx);
    double frac = idx - lo;
    return (lo < n - 1) ? sv(lo) * (1.0 - frac) + sv(lo + 1) * frac : sv(n - 1);
  };

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < slice_sz; ++i) {
    arma::vec sv = arma::sort(bootchol_flat.row(i).t());
    upper_vec(i) = pctile(sv, up_pct);
    lower_vec(i) = pctile(sv, low_pct);
  }

  arma::cube upper(upper_vec.memptr(), N, N, H, false, true);
  arma::cube lower(lower_vec.memptr(), N, N, H, false, true);

  BootstrapCholResult result;
  result.bootchol_flat = bootchol_flat;
  result.upper = upper;
  result.lower = lower;
  result.boot_beta = boot_beta;
  result.N = N;
  result.H = H;
  return result;
}

//' Bootstrap Cholesky Impulse Response Functions
//'
//' @param y T x N matrix of original endogenous variables.
//' @param var_result List from \code{fVAR}.
//' @param nboot Number of bootstrap replications.
//' @param horizon Maximum IRF horizon.
//' @param prc Confidence level in percent (e.g. 68).
//' @param bootscheme \code{"residual"} or \code{"wild"}.
//' @param exog Optional T x M exogenous matrix. Default NULL.
//' @param n_threads OpenMP threads. 0 = all cores minus one.
//'
//' @return List: \code{bootchol}, \code{upper}, \code{lower}, \code{boot_beta}.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fbootstrapChol(const arma::mat &y, const Rcpp::List &var_result,
                          int nboot, int horizon, double prc,
                          const std::string &bootscheme,
                          Rcpp::Nullable<arma::mat> exog = R_NilValue,
                          int n_threads = 0) {

  VARResult var_result_struct;
  var_result_struct.beta = Rcpp::as<arma::mat>(var_result["beta"]);
  var_result_struct.residuals = Rcpp::as<arma::mat>(var_result["residuals"]);
  var_result_struct.sigma_full = Rcpp::as<arma::mat>(var_result["sigma_full"]);
  var_result_struct.p = Rcpp::as<int>(var_result["p"]);
  var_result_struct.c = Rcpp::as<int>(var_result["c"]);

  if (var_result.containsElementNamed("n_exog")) {
    var_result_struct.n_exog = Rcpp::as<int>(var_result["n_exog"]);
  } else {
    var_result_struct.n_exog = 0;
  }

  BootstrapCholResult result = fbootstrapChol_cpp(
      y, var_result_struct, nboot, horizon, prc, bootscheme, exog, n_threads);

  Rcpp::NumericVector bootchol_out(result.bootchol_flat.begin(),
                                   result.bootchol_flat.end());
  bootchol_out.attr("dim") = Rcpp::IntegerVector::create(result.N, result.N, result.H, nboot);

  return Rcpp::List::create(Rcpp::Named("bootchol") = bootchol_out,
                            Rcpp::Named("upper") = result.upper,
                            Rcpp::Named("lower") = result.lower,
                            Rcpp::Named("boot_beta") = result.boot_beta);
}
