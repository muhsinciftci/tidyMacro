
#include "fcholeskyIRF.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Compute Cholesky Impulse Response Functions
//'
//' @param wold Wold representation cube (N x N x horizon)
//' @param S Cholesky factor matrix (N x N), lower triangular
//'
//' @return A cube containing:
//'   Cholesky structural impulse response functions (N x N x horizon)
//'
//' @details
//' This function computes structural impulse response functions using the
//' Cholesky decomposition identification scheme. For each horizon h, the
//' structural IRF is computed as:
//' \deqn{IRF_h = \Psi_h \cdot S}
//' where \eqn{\Psi_h} is the Wold representation at horizon h and S is the
//' Cholesky factor (lower triangular) of the covariance matrix of reduced-form
//' residuals.
//'
//' The Cholesky identification imposes a recursive structure on the
//' contemporaneous relationships between variables, with the ordering
//' of variables determining the causal structure.
//'
//' @examples
//' \dontrun{
//' # Generate sample Wold representation
//' N <- 3
//' horizon <- 20
//' wold <- array(rnorm(N * N * horizon), dim = c(N, N, horizon))
//'
//' # Compute Cholesky factor
//' Sigma <- matrix(c(1, 0.5, 0.3,
//'                   0.5, 1, 0.4,
//'                   0.3, 0.4, 1), 3, 3)
//' S <- t(chol(Sigma))  # Lower triangular
//'
//' # Compute Cholesky IRF
//' chol_irf <- fcholeskyIRF(wold, S)
//' }
//'
//' @export
// [[Rcpp::export]]
arma::cube fcholeskyIRF(const arma::cube &wold, const arma::mat &S) {

  int N = wold.n_rows;
  int horizon = wold.n_slices;

  // Pre-allocate output cube
  arma::cube chol(N, N, horizon);

  // Compute Cholesky IRF for each horizon
  for (int h = 0; h < horizon; h++) {
    chol.slice(h) = wold.slice(h) * S;
  }

  return chol;
}
