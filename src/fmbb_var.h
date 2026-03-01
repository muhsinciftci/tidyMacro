#ifndef FMBB_VAR_H
#define FMBB_VAR_H

#include <RcppArmadillo.h>

// Struct to hold MBB VAR results
struct MBBVARResult {
  arma::mat eps_boot;  // (T-p) x N matrix of bootstrapped residuals (centered)
  arma::mat M_boot;    // (T-p) x K matrix of bootstrapped instruments (centered)
};

// Internal C++ function - Nullable version (for use when M might be null)
MBBVARResult fmbb_var_cpp(const arma::mat& eps, 
                          int lags, 
                          int BlockSize, 
                          Rcpp::Nullable<arma::mat> M);

// Internal C++ function - Direct matrix version (for use in loops, thread-safe)
// This overload is more efficient and thread-safe (no Rcpp::wrap/as needed)
MBBVARResult fmbb_var_cpp(const arma::mat& eps, 
                          int lags, 
                          int BlockSize, 
                          const arma::mat& M);

// R wrapper function (for calling from R)
Rcpp::List fmbb_var(const arma::mat& eps, 
                    int lags, 
                    int BlockSize, 
                    Rcpp::Nullable<arma::mat> M);

#endif // FMBB_VAR_H
