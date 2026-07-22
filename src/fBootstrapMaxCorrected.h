#ifndef FBOOTSTRAPMAXCORRECTED_H
#define FBOOTSTRAPMAXCORRECTED_H

#include "fBootstrapMax.h"
#include "fVAR.h"
#include <RcppArmadillo.h>

BootstrapMaxResult
fBootstrapMaxCorrected_cpp(const arma::mat& y, const VARResult& var_result,
                           int nboot1, int nboot2, int horizon, int var_idx,
                           double prc, double prc2, const arma::uvec& cumulate,
                           Rcpp::Nullable<arma::vec> scaling,
                           Rcpp::Nullable<arma::mat> exog,
                           int n_threads);

Rcpp::List fBootstrapMaxCorrected(const arma::mat& y, const Rcpp::List& var_result,
                                   int nboot1, int nboot2, int horizon, int var_idx,
                                   double prc, double prc2, const arma::uvec& cumulate,
                                   Rcpp::Nullable<arma::vec> scaling,
                                   Rcpp::Nullable<arma::mat> exog,
                                   int n_threads);

#endif // FBOOTSTRAPMAXCORRECTED_H
