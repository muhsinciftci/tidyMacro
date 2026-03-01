#ifndef FHISTDEC_H
#define FHISTDEC_H

#include <RcppArmadillo.h>
#include "fVAR.h"

// Struct to hold historical decomposition results
struct HistDecResult {
    arma::mat histdec;  // (T-p) x N matrix: contribution of each shock to the series
    arma::vec ystar;    // (T-p) x 1 demeaned series
};

// Internal C++ function (for use in other C++ code)
HistDecResult fhistdec_cpp(const arma::mat& y,
                           const VARResult& var_result,
                           const arma::mat& K,
                           int series);

// R wrapper
Rcpp::List fhistdec(const arma::mat& y,
                    const Rcpp::List& fVAR,
                    const arma::mat& K,
                    int series);

#endif // FHISTDEC_H
