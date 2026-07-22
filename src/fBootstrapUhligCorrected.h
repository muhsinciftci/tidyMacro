#ifndef FBOOTSTRAPUHLIGCORRECTED_H
#define FBOOTSTRAPUHLIGCORRECTED_H

#include "fBootstrapUhlig.h"
#include "fVAR.h"
#include <RcppArmadillo.h>

BootstrapUhligResult
fBootstrapUhligCorrected_cpp(const arma::mat& y, const VARResult& var_result,
                              int nboot1, int nboot2, int horizon, int idx,
                              double prc, double prc2, const arma::uvec& cumulate,
                              Rcpp::Nullable<arma::mat> exog,
                              int n_threads);

Rcpp::List fBootstrapUhligCorrected(const arma::mat& y, const Rcpp::List& var_result,
                                     int nboot1, int nboot2, int horizon, int idx,
                                     double prc, double prc2, const arma::uvec& cumulate,
                                     Rcpp::Nullable<arma::mat> exog,
                                     int n_threads);

#endif // FBOOTSTRAPUHLIGCORRECTED_H
