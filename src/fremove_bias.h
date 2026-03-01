#ifndef FREMOVE_BIAS_H
#define FREMOVE_BIAS_H

#include <RcppArmadillo.h>

// Struct to hold remove bias results
struct RemoveBiasResult {
  arma::mat Beta;       // N x (Np+c+M) bias-corrected coefficient matrix
  int corrections;      // Integer count of iterative shrinkage steps
};

// Internal C++ function (for use in other C++ code)
RemoveBiasResult fremove_bias_cpp(const arma::mat& beta,
                                  int c,
                                  int p,
                                  const arma::cube& boot_beta);

// R wrapper function (for calling from R)
Rcpp::List fremove_bias(const arma::mat& beta,
                        int c,
                        int p,
                        const arma::cube& boot_beta);

#endif // FREMOVE_BIAS_H
