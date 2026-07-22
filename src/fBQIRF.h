#ifndef FBQIRF_H
#define FBQIRF_H

#include <RcppArmadillo.h>

// Function declaration — default value defined in fBQIRF.cpp, not here
arma::cube fBQIRF(const arma::cube& wold, const arma::mat& K,
                  Rcpp::Nullable<arma::vec> scaling);

// Pre-parsed scaling overload (thread-safe, no Rcpp::as).
// has_scaling=false => no normalisation; otherwise divide entire cube by
// bqirf(idx_0based, idx_0based, 0) / shock_size.
arma::cube fBQIRF_cpp(const arma::cube& wold, const arma::mat& K,
                      bool has_scaling, int idx_0based, double shock_size);

#endif // FBQIRF_H
