#ifndef FBOOTSTRAPUHLIG_H
#define FBOOTSTRAPUHLIG_H

#include "fVAR.h"
#include <RcppArmadillo.h>

struct BootstrapUhligResult {
    arma::mat bootuhlig_flat; // (N*(horizon+1)) x nboot flattened
    arma::mat upper;          // N x (horizon+1) upper bands (prc)
    arma::mat lower;          // N x (horizon+1) lower bands (prc)
    arma::mat upper2;         // N x (horizon+1) upper bands (prc2)
    arma::mat lower2;         // N x (horizon+1) lower bands (prc2)
    arma::cube boot_beta;     // N x n_coef x nboot bootstrapped coefficients
    int N;
    int H;
};

BootstrapUhligResult
fbootstrapUhlig_cpp(const arma::mat& y, const VARResult& var_result,
                    int nboot, int horizon, int idx, double prc, double prc2,
                    const arma::uvec& cumulate,
                    Rcpp::Nullable<arma::mat> exog,
                    int n_threads);

Rcpp::List fbootstrapUhlig(const arma::mat& y, const Rcpp::List& var_result,
                            int nboot, int horizon, int idx, double prc, double prc2,
                            const arma::uvec& cumulate,
                            Rcpp::Nullable<arma::mat> exog,
                            int n_threads);

#endif // FBOOTSTRAPUHLIG_H
