#ifndef FVAR_H
#define FVAR_H

#include <RcppArmadillo.h>

// Struct to hold VAR estimation results
struct VARResult {
  arma::mat beta;        // Coefficient matrix (Np + c + M) x N
  arma::mat residuals;   // (T-p) x N matrix of OLS residuals
  arma::mat sigma_full;  // N x N residual covariance matrix
  int p;                 // Lag order
  int c;                 // Intercept indicator
  int n_exog;            // Number of exogenous variables
};

// Internal C++ function (for use in other C++ code)
VARResult fVAR_cpp(const arma::mat& y, int p, int c,
                   Rcpp::Nullable<arma::mat> exog);

// R wrapper function (for calling from R)
Rcpp::List fVAR(const arma::mat& y, int p, int c,
                Rcpp::Nullable<arma::mat> exog);

#endif // FVAR_H
