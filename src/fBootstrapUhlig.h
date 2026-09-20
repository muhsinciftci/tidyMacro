#ifndef FBOOTSTRAPUHLIG_H
#define FBOOTSTRAPUHLIG_H

#include "fVAR.h"
#include <RcppArmadillo.h>

struct BootstrapUhligResult {
    arma::mat bootuhlig_flat; // (N*(horizon+1)) x nboot flattened
    arma::mat upper;          // N x (horizon+1) upper bands (conf)
    arma::mat lower;          // N x (horizon+1) lower bands (conf)
    arma::mat upper2;         // N x (horizon+1) upper bands (conf2)
    arma::mat lower2;         // N x (horizon+1) lower bands (conf2)
    arma::cube boot_beta;     // N x n_coef x nboot bootstrapped coefficients
    int N;
    int H;
};

BootstrapUhligResult
fBootstrapUhlig_cpp(const arma::mat& y, const VARResult& var_result,
                    int nboot, int horizon, int idx, double conf, double conf2,
                    const arma::uvec& cumulate,
                    Rcpp::Nullable<arma::mat> exog,
                    int n_threads);

Rcpp::List fBootstrapUhlig(const arma::mat& y, const Rcpp::List& var_result,
                            int nboot, int horizon, int idx, double conf, double conf2,
                            const arma::uvec& cumulate,
                            Rcpp::Nullable<arma::mat> exog,
                            int n_threads);

#endif // FBOOTSTRAPUHLIG_H
