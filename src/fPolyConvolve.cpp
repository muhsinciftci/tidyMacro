#include "fPolyConvolve.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::cube fPolyConvolve_cpp(const arma::cube& A, const arma::cube& B, int nlags) {
    int A1    = A.n_rows;
    int A3    = A.n_slices;
    int B2    = B.n_cols;
    int B3    = B.n_slices;
    int total = A3 + B3 - 1;

    arma::cube C(A1, B2, total, arma::fill::zeros);
    for (int a = 0; a < A3; a++)
        for (int b = 0; b < B3; b++)
            C.slice(a + b) += A.slice(a) * B.slice(b);

    return C.slices(0, nlags - 1);
}

//' Polynomial Convolution of Two Matrix Lag Polynomials
//'
//' Computes \eqn{C(L) = A(L) \cdot B(L)} truncated to \code{nlags} terms.
//' Both polynomials are stored as 3-D arrays where the third dimension indexes
//' lag order (slice 1 = lag 0, slice 2 = lag 1, etc.).
//'
//' @param A 3-D array of dimension \eqn{d_1 \times d_2 \times A_3}
//'   (e.g. Wold IRF coefficients from \code{fwoldIRF}).
//' @param B 3-D array of dimension \eqn{d_2 \times K \times B_3}
//'   (e.g. \eqn{\Psi(L)} coefficients).
//' @param nlags Integer. Number of output slices to return
//'   (must be \eqn{\leq A_3 + B_3 - 1}).
//'
//' @return A 3-D array of dimension \eqn{d_1 \times K \times} \code{nlags}.
//'
//' @seealso \code{\link{fwoldIRF}}, \code{\link{fcholeskyIRF}}
//'
//' @export
// [[Rcpp::export]]
arma::cube fPolyConvolve(const arma::cube& A, const arma::cube& B, int nlags) {
    return fPolyConvolve_cpp(A, B, nlags);
}
