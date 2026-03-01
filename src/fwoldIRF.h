#ifndef FWOLDIRF_H
#define FWOLDIRF_H

#include <RcppArmadillo.h>
#include "fVAR.h"  // Include to access VARResult struct

// Struct to hold Wold IRF results
struct WoldIRFResult {
  arma::cube irfwold; // N x N x (horizon+1) cube of impulse response functions
};

// Internal C++ function (for use in other C++ code)
// Accepts VARResult struct directly for efficient C++ to C++ calls
WoldIRFResult fwoldIRF_cpp(const VARResult& var_result, int horizon);

// R wrapper function (for calling from R)
arma::cube fwoldIRF(const Rcpp::List& fVAR, int horizon);

#endif // FWOLDIRF_H
