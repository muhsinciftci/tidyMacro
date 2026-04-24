#ifndef FBOOTSTRAPCHOL_H
#define FBOOTSTRAPCHOL_H

#include "fVAR.h"
#include <RcppArmadillo.h>

// Struct to hold bootstrap Cholesky results
struct BootstrapCholResult {
  arma::mat bootchol_flat; // slice_sz x nboot matrix (flattened)
  arma::cube upper;        // N x N x (horizon+1) upper confidence bands (prc)
  arma::cube lower;        // N x N x (horizon+1) lower confidence bands (prc)
  arma::cube upper2;       // N x N x (horizon+1) upper confidence bands (prc2)
  arma::cube lower2;       // N x N x (horizon+1) lower confidence bands (prc2)
  arma::cube boot_beta;    // N x n_coef x nboot bootstrapped coefficients
  int N;                   // Number of variables
  int H;                   // Horizon + 1
};

// Internal C++ function (for use in other C++ code)
BootstrapCholResult
fbootstrapChol_cpp(const arma::mat &y, const VARResult &var_result, int nboot,
                   int horizon, double prc, double prc2,
                   const std::string &bootscheme,
                   Rcpp::Nullable<arma::mat> exog, int n_threads);

// R wrapper function (for calling from R)
Rcpp::List fbootstrapChol(const arma::mat &y, const Rcpp::List &var_result,
                          int nboot, int horizon, double prc, double prc2,
                          const std::string &bootscheme,
                          Rcpp::Nullable<arma::mat> exog, int n_threads);

#endif // FBOOTSTRAPCHOL_H
