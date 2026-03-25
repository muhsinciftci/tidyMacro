#ifndef FBOOTSTRAPMAX_H
#define FBOOTSTRAPMAX_H

#include "fVAR.h"
#include <RcppArmadillo.h>

struct BootstrapMaxResult {
    arma::mat bootmax_flat; // (N*(horizon+1)) x nboot flattened
    arma::mat upper;        // N x (horizon+1) upper bands
    arma::mat lower;        // N x (horizon+1) lower bands
    arma::cube boot_beta;   // N x n_coef x nboot bootstrapped coefficients
    int N;
    int H;
};

BootstrapMaxResult
fbootstrapMax_cpp(const arma::mat& y, const VARResult& var_result,
                  int nboot, int horizon, int var_idx, double prc,
                  const arma::uvec& cumulate,
                  Rcpp::Nullable<arma::vec> scaling,
                  int n_threads);

Rcpp::List fbootstrapMax(const arma::mat& y, const Rcpp::List& var_result,
                         int nboot, int horizon, int var_idx, double prc,
                         const arma::uvec& cumulate,
                         Rcpp::Nullable<arma::vec> scaling,
                         int n_threads);

#endif // FBOOTSTRAPMAX_H
