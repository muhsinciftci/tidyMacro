#ifndef FEVD_IV_H
#define FEVD_IV_H

#include <RcppArmadillo.h>

// Struct to hold FEVD IV results
struct FEVDIVResult {
  arma::vec scaler;              // N x 1 vector of scaling factors
  arma::mat ivirf_scaled;        // N x (hor+1) matrix of unit variance impulse responses
  arma::rowvec oil_news_unit_var; // 1 x T1 vector of unit variance shock series
  arma::mat check_uv;            // Scalar variance check
  arma::mat denom;               // N x (hor+1) matrix of total forecast error variances
  arma::mat fevd_iv;             // N x (hor+1) matrix of FEVD shares
};

// Internal C++ function (for use in other C++ code)
FEVDIVResult fevd_iv_cpp(arma::vec s, arma::mat S, arma::cube wold, int N, int hor,
                         arma::mat sigma, arma::mat u, int T1, int p);

// R wrapper function (for calling from R)
Rcpp::List fevd_iv(arma::vec s, arma::mat S, arma::cube wold, int N, int hor,
                   arma::mat sigma, arma::mat u, int T1, int p);

#endif // FEVD_IV_H
