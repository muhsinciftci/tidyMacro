#ifndef FBOOTSTRAPCHOL_H
#define FBOOTSTRAPCHOL_H

#include "fVAR.h" // Include to access VARResult struct
#include <RcppArmadillo.h>

// Struct to hold bootstrap Cholesky results
struct BootstrapCholResult {
  arma::mat bootchol_flat; // slice_sz x nboot matrix (flattened)
  arma::cube upper;        // N x N x (horizon+1) upper confidence bands
  arma::cube lower;        // N x N x (horizon+1) lower confidence bands
  arma::cube boot_beta;    // N x n_coef x nboot bootstrapped coefficients
  int N;                   // Number of variables
  int H;                   // Horizon + 1
};

// Internal C++ function (for use in other C++ code)
// Accepts VARResult struct directly for efficient C++ to C++ calls
BootstrapCholResult
fbootstrapChol_cpp(const arma::mat &y, const VARResult &var_result, int nboot,
                   int horizon, double prc, const std::string &bootscheme,
                   Rcpp::Nullable<arma::mat> exog, int n_threads);

// R wrapper function (for calling from R)
Rcpp::List fbootstrapChol(const arma::mat &y, const Rcpp::List &var_result,
                          int nboot, int horizon, double prc,
                          const std::string &bootscheme,
                          Rcpp::Nullable<arma::mat> exog, int n_threads);

#endif // FBOOTSTRAPCHOL_H
