#ifndef FBQIRF_H
#define FBQIRF_H

#include <RcppArmadillo.h>

// Function declaration — default value defined in fbqIRF.cpp, not here
arma::cube fbqIRF(const arma::cube& wold, const arma::mat& K,
                  Rcpp::Nullable<arma::vec> scaling);

#endif // FBQIRF_H
