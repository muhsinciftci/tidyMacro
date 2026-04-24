// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include "fbootstrapCholCorrected.h"
#include "fbootstrapChol.h"
#include "fVAR.h"
#include "fbootstrapVAR.h"
#include "fcompanionMatrix.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstdio>
#include <algorithm>

// Bias-correct beta using the Kilian shortcut.
// beta      : (n_coef x N) coefficient matrix (as stored in VARResult)
// boot_mean : N x n_coef matrix — mean of the first-pass bootstrap estimates
//             (streamed online in Pass 1; no full cube needed)
// Beta_t    : output — bias-corrected beta transposed (N x n_coef)
// corrections: output — number of shrinkage iterations applied
static void bias_correct_chol(const arma::mat &beta, int c, int p,
                               const arma::mat &boot_mean,
                               arma::mat &Beta_t, int &corrections) {
  const arma::mat beta_t = beta.t();
  const arma::mat bias   = boot_mean - beta_t;

  // Check stability of original estimates
  CompanionMatrixResult comp0 = fcompanionMatrix_cpp(beta, c, p);
  arma::cx_vec ev0 = arma::eig_gen(comp0.comp);
  double max_ev    = arma::max(arma::abs(ev0));

  corrections = 1;

  if (max_ev >= 1.0) {
    // Original estimates already unstable — no correction
    Beta_t = beta_t;
    return;
  }

  Beta_t = beta_t - bias;

  {
    CompanionMatrixResult comp1 = fcompanionMatrix_cpp(Beta_t.t(), c, p);
    arma::cx_vec ev1 = arma::eig_gen(comp1.comp);
    max_ev = arma::max(arma::abs(ev1));
  }

  // Shrink bias towards zero until corrected estimates are stable
  double delta = 1.0;
  while (max_ev >= 1.0) {
    delta -= 0.01;
    corrections += 1;
    if (delta < 0.0 || corrections > 200) {
      Beta_t = beta_t;
      break;
    }
    Beta_t = beta_t - bias * delta;
    CompanionMatrixResult comp_loop = fcompanionMatrix_cpp(Beta_t.t(), c, p);
    arma::cx_vec ev_loop = arma::eig_gen(comp_loop.comp);
    max_ev = arma::max(arma::abs(ev_loop));
  }
}

// Internal C++ function (called from other C++ code)
BootstrapCholResult
fbootstrapCholCorrected_cpp(const arma::mat &y, const VARResult &var_result,
                             int nboot1, int nboot2, int horizon, double prc,
                             double prc2, const std::string &bootscheme,
                             Rcpp::Nullable<arma::mat> exog, int n_threads = 0) {

  const int p      = var_result.p;
  const int c      = var_result.c;
  const int n_exog = var_result.n_exog;
  const int N      = static_cast<int>(y.n_cols);
  const int n_coef = static_cast<int>(var_result.beta.n_rows);

  if (n_exog > 0 && exog.isNull())
    Rcpp::stop("Original VAR used exogenous variables. You must provide the 'exog' parameter.");
  if (n_exog == 0 && exog.isNotNull())
    Rcpp::stop("Original VAR did not use exogenous variables. Do not provide the 'exog' parameter.");

  int actual_threads = 1;
#ifdef _OPENMP
  actual_threads = (n_threads == 0)
                       ? std::max(1, omp_get_max_threads() - 1)
                       : n_threads;
  omp_set_num_threads(actual_threads);
  std::printf("[Pass 1] Bias estimation: %d reps, %d thread(s)\n",
              nboot1, actual_threads);
#else
  std::printf("[Pass 1] Bias estimation: %d reps, single-threaded\n", nboot1);
#endif

  // ------------------------------------------------------------------ //
  // Pass 1: stream the coefficient mean online — no full cube
  // ------------------------------------------------------------------ //
  const bool has_exog = exog.isNotNull();
  arma::mat exog_mat;
  if (has_exog) exog_mat = Rcpp::as<arma::mat>(exog);

  arma::mat boot_mean(N, n_coef, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp parallel
  {
    arma::mat thread_sum(N, n_coef, arma::fill::zeros);
#pragma omp for schedule(dynamic) nowait
    for (int b = 0; b < nboot1; ++b) {
      BootstrapVARResult boot_data = fbootstrapVAR_cpp(y, var_result, bootscheme);
      VARResult var_loop = has_exog
                             ? fVAR_cpp(boot_data.ynext, p, c, exog_mat)
                             : fVAR_cpp(boot_data.ynext, p, c, R_NilValue);
      thread_sum += var_loop.beta.t();
    }
#pragma omp critical
    boot_mean += thread_sum;
  }
#else
  for (int b = 0; b < nboot1; ++b) {
    BootstrapVARResult boot_data = fbootstrapVAR_cpp(y, var_result, bootscheme);
    VARResult var_loop = has_exog
                           ? fVAR_cpp(boot_data.ynext, p, c, exog_mat)
                           : fVAR_cpp(boot_data.ynext, p, c, R_NilValue);
    boot_mean += var_loop.beta.t();
  }
#endif
  boot_mean /= static_cast<double>(nboot1);

  // ------------------------------------------------------------------ //
  // Bias correction (Kilian shortcut)
  // ------------------------------------------------------------------ //
  arma::mat Beta_t;
  int corrections = 1;
  bias_correct_chol(var_result.beta, c, p, boot_mean, Beta_t, corrections);

  if (corrections > 1)
    std::printf("[Bias correction] %d shrinkage iteration(s)\n", corrections);
  else
    std::printf("[Bias correction] Full correction applied\n");

  // ------------------------------------------------------------------ //
  // Pass 2: bootstrap with bias-corrected beta
  // ------------------------------------------------------------------ //
#ifdef _OPENMP
  std::printf("[Pass 2] Bias-corrected bootstrap: %d reps, %d thread(s)\n",
              nboot2, actual_threads);
#else
  std::printf("[Pass 2] Bias-corrected bootstrap: %d reps, single-threaded\n", nboot2);
#endif

  VARResult corrected_var = var_result;
  corrected_var.beta      = Beta_t.t();

  return fbootstrapChol_cpp(y, corrected_var, nboot2, horizon, prc, prc2,
                             bootscheme, exog, actual_threads);
}

//' Bootstrap Bias-Corrected Cholesky Identified Impulse Response Functions
//'
//' @param y T x N matrix of original endogenous variables.
//' @param var_result List from \code{fVAR()}.
//' @param nboot1 Number of first-pass replications for bias estimation.
//' @param nboot2 Number of second-pass replications for confidence bands.
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
//'     \item{bootchol}{N x N x (horizon+1) x nboot2 array of bootstrapped Cholesky IRFs}
//'     \item{upper}{N x N x (horizon+1) outer upper confidence bands}
//'     \item{lower}{N x N x (horizon+1) outer lower confidence bands}
//'     \item{upper2}{N x N x (horizon+1) inner upper confidence bands}
//'     \item{lower2}{N x N x (horizon+1) inner lower confidence bands}
//'     \item{boot_beta}{N x (Np+c) x nboot2 array of bootstrapped coefficients}
//'   }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fbootstrapCholCorrected(const arma::mat &y, const Rcpp::List &var_result,
                                    int nboot1, int nboot2, int horizon,
                                    double prc = 90.0, double prc2 = 68.0,
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

  BootstrapCholResult result = fbootstrapCholCorrected_cpp(y, var_result_struct,
                                                            nboot1, nboot2, horizon,
                                                            prc, prc2, bootscheme,
                                                            exog, n_threads);

  Rcpp::NumericVector bootchol_out(result.bootchol_flat.begin(),
                                   result.bootchol_flat.end());
  bootchol_out.attr("dim") = Rcpp::IntegerVector::create(result.N, result.N,
                                                          result.H, nboot2);

  return Rcpp::List::create(Rcpp::Named("bootchol")   = bootchol_out,
                            Rcpp::Named("upper")      = result.upper,
                            Rcpp::Named("lower")      = result.lower,
                            Rcpp::Named("upper2")     = result.upper2,
                            Rcpp::Named("lower2")     = result.lower2,
                            Rcpp::Named("boot_beta")  = result.boot_beta);
}
