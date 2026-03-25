#ifndef FUHLIG_MAXSHARE_H
#define FUHLIG_MAXSHARE_H

#include <RcppArmadillo.h>

// Internal C++ function
// wold: N x N x H cube (H = horizon+1), S: N x N lower Cholesky, idx: 0-based
arma::vec fuhlig_maxshare_cpp(const arma::cube& wold, const arma::mat& S, int idx);

// R wrapper (idx is 1-based)
arma::vec fuhlig_maxshare(const arma::cube& wold, const arma::mat& S, int idx);

#endif // FUHLIG_MAXSHARE_H
