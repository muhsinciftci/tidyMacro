#ifndef FVARX_H
#define FVARX_H

#include <RcppArmadillo.h>

// Struct to hold VARX estimation results
struct VARXResult {
  arma::mat beta;      // Coefficient matrix (Np + M + Mp + c) x N
  arma::mat residuals; // (T-p) x N residual matrix
};

// Internal C++ function (for use in other C++ code)
VARXResult fVARX_cpp(const arma::mat& y, const arma::mat& ex, int p, int c);

// R wrapper function (for calling from R)
Rcpp::List fVARX(const arma::mat& y, const arma::mat& ex, int p, int c);

#endif // FVARX_H
