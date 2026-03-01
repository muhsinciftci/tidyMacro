#ifndef FBOOTSTRAPVAR_H
#define FBOOTSTRAPVAR_H

#include <RcppArmadillo.h>
#include "fVAR.h"  // Include to access VARResult struct

// Struct to hold bootstrap VAR results
struct BootstrapVARResult {
  arma::mat ynext;      // T x N bootstrapped data matrix
  arma::vec rademacher; // T_eff-length vector of Rademacher signs
};

// Internal C++ function (for use in other C++ code)
// Accepts VARResult struct directly for efficient C++ to C++ calls
BootstrapVARResult fbootstrapVAR_cpp(const arma::mat& y, 
                                     const VARResult& var_result,
                                     const std::string& bootscheme);

// R wrapper function (for calling from R)
Rcpp::List fbootstrapVAR(const arma::mat& y, 
                         const Rcpp::List& fVAR_result, 
                         const std::string& bootscheme);

#endif
