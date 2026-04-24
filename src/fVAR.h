#ifndef FVAR_H
#define FVAR_H

#include <RcppArmadillo.h>

// Struct to hold VAR estimation results
struct VARResult {
  arma::mat beta;        // Coefficient matrix (Np + c + M) x N
  arma::mat residuals;   // (T-p) x N matrix of OLS residuals
  arma::mat sigma;       // N x N residual covariance matrix (DF-corrected)
  int p;                 // Lag order
  int c;                 // Intercept indicator
  int n_exog;            // Number of exogenous variables
};

// Internal C++ function (for use in other C++ code) — Nullable-exog overload
VARResult fVAR_cpp(const arma::mat& y, int p, int c,
                   Rcpp::Nullable<arma::mat> exog);

// Internal C++ function — pre-converted exog overload (thread-safe, no Rcpp::as)
// Use this from bootstrap loops where exog has already been converted once
// outside the parallel region.
VARResult fVAR_cpp(const arma::mat& y, int p, int c,
                   const arma::mat& exog);


// R wrapper function (for calling from R)
Rcpp::List fVAR(const arma::mat& y, int p, int c,
                Rcpp::Nullable<arma::mat> exog);

#endif // FVAR_H
