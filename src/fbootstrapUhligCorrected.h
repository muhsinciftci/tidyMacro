#ifndef FBOOTSTRAPUHLIGCORRECTED_H
#define FBOOTSTRAPUHLIGCORRECTED_H

#include "fbootstrapUhlig.h"
#include "fVAR.h"
#include <RcppArmadillo.h>

BootstrapUhligResult
fbootstrapUhligCorrected_cpp(const arma::mat& y, const VARResult& var_result,
                              int nboot1, int nboot2, int horizon, int idx,
                              double prc, const arma::uvec& cumulate,
                              int n_threads);

Rcpp::List fbootstrapUhligCorrected(const arma::mat& y, const Rcpp::List& var_result,
                                     int nboot1, int nboot2, int horizon, int idx,
                                     double prc, const arma::uvec& cumulate,
                                     int n_threads);

#endif // FBOOTSTRAPUHLIGCORRECTED_H
