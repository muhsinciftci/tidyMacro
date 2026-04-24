#ifndef FBOOTSTRAPHOLCORRECTED_H
#define FBOOTSTRAPHOLCORRECTED_H

#include "fbootstrapChol.h"
#include "fVAR.h"
#include <RcppArmadillo.h>

// Internal C++ function (for use in other C++ code)
BootstrapCholResult
fbootstrapCholCorrected_cpp(const arma::mat &y, const VARResult &var_result,
                             int nboot1, int nboot2, int horizon, double prc,
                             double prc2, const std::string &bootscheme,
                             Rcpp::Nullable<arma::mat> exog, int n_threads);

// R wrapper function (for calling from R)
Rcpp::List fbootstrapCholCorrected(const arma::mat &y, const Rcpp::List &var_result,
                                    int nboot1, int nboot2, int horizon, double prc,
                                    double prc2, const std::string &bootscheme,
                                    Rcpp::Nullable<arma::mat> exog, int n_threads);

#endif // FBOOTSTRAPHOLCORRECTED_H
