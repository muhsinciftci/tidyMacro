#ifndef FBOOTSTRAPBQ_H
#define FBOOTSTRAPBQ_H

#include "fVAR.h"
#include <RcppArmadillo.h>

// Struct to hold bootstrap BQ results
struct BootstrapBQResult {
  arma::mat bootbq_flat; // slice_sz x nboot matrix (flattened)
  arma::cube upper;      // N x N x (horizon+1) upper confidence bands
  arma::cube lower;      // N x N x (horizon+1) lower confidence bands
  arma::cube boot_beta;  // N x n_coef x nboot bootstrapped coefficients
  int N;                 // Number of variables
  int H;                 // Horizon + 1
};

// Internal C++ function (for use in other C++ code)
BootstrapBQResult
fbootstrapBQ_cpp(const arma::mat &y, const VARResult &var_result, int nboot,
                 int horizon, double prc, const std::string &bootscheme,
                 const arma::uvec &cumulate,
                 Rcpp::Nullable<arma::vec> scaling, int n_threads);

// R wrapper function (for calling from R)
Rcpp::List fbootstrapBQ(const arma::mat &y, const Rcpp::List &var_result,
                        int nboot, int horizon, double prc,
                        const std::string &bootscheme,
                        const arma::uvec &cumulate,
                        Rcpp::Nullable<arma::vec> scaling, int n_threads);

#endif // FBOOTSTRAPBQ_H
