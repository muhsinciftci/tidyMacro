#ifndef FEVD_CHOL_H
#define FEVD_CHOL_H

#include <RcppArmadillo.h>

// Struct to hold Cholesky FEVD results
struct FEVDCholResult {
  arma::cube fevd;  // N x N x (horizon+1) forecast error variance decomposition
                    // fevd(i, j, h) = share of variance of variable i at horizon h
                    // explained by shock j
};

// Internal C++ function (for use in other C++ code)
FEVDCholResult fevd_chol_cpp(const arma::cube& chol_irf);

// R wrapper function (for calling from R)
// shock: 1-indexed. If 0 (default), returns full 3D cube; otherwise returns N x horizon matrix.
Rcpp::List fevd_chol(const arma::cube& chol_irf, int shock);

#endif // FEVD_CHOL_H
