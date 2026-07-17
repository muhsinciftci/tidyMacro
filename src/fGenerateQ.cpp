#include <RcppArmadillo.h>
#include "fGenerateQ.h"
// [[Rcpp::depends(RcppArmadillo)]]

void fGenerateQ_inplace(arma::mat& Q,
                        arma::mat& R,
                        arma::mat& G,
                        const arma::uword N) {
    if (N == 0u) {
        Rcpp::stop("N must be positive.");
    }

    G.randn(N, N);
    arma::qr_econ(Q, R, G);

    for (arma::uword j = 0; j < N; ++j) {
        if (R(j, j) < 0.0) {
            Q.col(j) *= -1.0;
        }
    }
}

arma::mat fGenerateQ_cpp(int N) {
    if (N <= 0) {
        Rcpp::stop("N must be positive.");
    }

    arma::mat Q;
    arma::mat R;
    arma::mat G;
    fGenerateQ_inplace(Q, R, G, static_cast<arma::uword>(N));
    return Q;
}

//' Draw a Random Orthonormal Matrix (RZW 2010)
//'
//' Generates an N x N orthonormal matrix Q (QQ' = Q'Q = I) by applying QR
//' decomposition to a matrix of i.i.d. standard normals and normalising the
//' diagonal of R to be positive, yielding a unique decomposition.
//'
//' @param N Integer. Dimension of the square matrix.
//'
//' @return An N x N orthonormal matrix.
//'
//' @details
//' Draw M ~ N(0, I_{N x N}), compute the QR decomposition M = QR, then set
//' Q[,i] = -Q[,i] whenever R[i,i] < 0.  The resulting Q is uniformly
//' distributed on the Stiefel manifold (Haar measure), as required for
//' sign-restriction identification.
//'
//' @references
//' Rubio-Ramirez, J. F., Waggoner, D. F., & Zha, T. (2010). Structural vector
//' autoregressions: Theory of identification and algorithms for inference.
//' \emph{Review of Economic Studies}, 77(2), 665--696.
//'
//' @examples
//' \dontrun{
//' Q <- fGenerateQ(3)
//' round(t(Q) %*% Q, 10)  # should be identity
//' }
//'
//' @export
// [[Rcpp::export]]
arma::mat fGenerateQ(int N) {
    return fGenerateQ_cpp(N);
}
