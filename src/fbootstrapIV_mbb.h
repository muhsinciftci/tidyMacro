#ifndef FBOOTSTRAPIV_MBB_H
#define FBOOTSTRAPIV_MBB_H

#include "fVAR.h" // Include to access VARResult struct
#include <RcppArmadillo.h>

// Struct to hold bootstrap IV MBB results
struct BootstrapIVMBBResult {
  arma::mat upper;     // N x (horizon+1) matrix of upper confidence bands
  arma::mat lower;     // N x (horizon+1) matrix of lower confidence bands
  arma::mat meanirf;   // N x (horizon+1) matrix of mean impulse responses
  arma::mat medianirf; // N x (horizon+1) matrix of median impulse responses
};

// Internal C++ function (for use in other C++ code)
// Accepts VARResult struct directly for efficient C++ to C++ calls
BootstrapIVMBBResult
fbootstrapIV_mbb_cpp(const arma::mat &y, const VARResult &var_result,
                     const arma::mat &Z, int nboot, int blocksize,
                     const arma::ivec &adjustZ, const arma::ivec &adjustu,
                     int policyvar, int horizon, double prc,
                     Rcpp::Nullable<arma::mat> exog, int n_threads);

// R wrapper function (for calling from R)
Rcpp::List fbootstrapIV_mbb(const arma::mat &y, const Rcpp::List &var_result,
                            const arma::mat &Z, int nboot, int blocksize,
                            const arma::ivec &adjustZ,
                            const arma::ivec &adjustu, int policyvar,
                            int horizon, double prc,
                            Rcpp::Nullable<arma::mat> exog, int n_threads);

#endif // FBOOTSTRAPIV_MBB_H
