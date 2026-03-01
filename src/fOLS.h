#ifndef FOLS_H
#define FOLS_H

#include <RcppArmadillo.h>

// Struct to hold OLS estimation results
struct OLSResult {
  arma::mat beta;           // Coefficient estimates
  arma::mat fitted;         // Fitted values
  arma::mat err;            // Residuals
  double r2;                // R-squared statistic
  arma::mat fitted_partial; // Fitted values excluding intercept
};

// Internal C++ function (for use in other C++ code)
OLSResult fOLS_cpp(arma::mat y, arma::mat X, int c);

// R wrapper function (for calling from R)
Rcpp::List fOLS(arma::mat y, arma::mat X, int c);

#endif // FOLS_H
