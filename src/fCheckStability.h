#ifndef FCHECKSTABILITY_H
#define FCHECKSTABILITY_H

#include <RcppArmadillo.h>
#include "fVAR.h"  // Include to access VARResult struct
// [[Rcpp::depends(RcppArmadillo)]]

// Struct to hold VAR stability results
struct StabilityResult {
  double max_eig;  // Largest eigenvalue modulus max_i |lambda_i|
  bool stable;     // true if max_eig < 1
};

// Internal C++ function (for use in other C++ code)
// Accepts VARResult struct directly for efficient C++ to C++ calls
StabilityResult fCheckStability_cpp(const VARResult& var_result);

// Overload for callers that hold raw coefficients rather than a VARResult
StabilityResult fCheckStability_cpp(const arma::mat& beta, int c, int p);

// R wrapper function (for calling from R) — prints, returns nothing
void fCheckStability(const Rcpp::List& var_model);

#endif // FCHECKSTABILITY_H
