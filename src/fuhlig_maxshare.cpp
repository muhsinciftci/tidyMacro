// [[Rcpp::depends(RcppArmadillo)]]
#include "fuhlig_maxshare.h"
#include <RcppArmadillo.h>

// Replicates uhlig_maxshare.m (Barsky & Sims 2012 / Uhlig 2004 max-share criterion).
//
// Builds Omega = sum_{h=0}^{H-1} cumsum_h  where cumsum_h = sum_{k=0}^{h} D_k'D_k
// and D_k = wold.slice(k).row(idx) * S  (1 x N).
// This equals sum_{k=0}^{H-1} (H-k) * D_k'*D_k — a horizon-weighted FEV matrix.
// The leading eigenvector of the lower-right (N-1)x(N-1) block gives the free
// parameters; h2 = [0; eigenvec] imposes the zero-impact restriction.

arma::vec fuhlig_maxshare_cpp(const arma::cube& wold, const arma::mat& S, int idx) {
    const int H = static_cast<int>(wold.n_slices);
    const int N = static_cast<int>(wold.n_rows);

    arma::mat omega(N, N, arma::fill::zeros);
    arma::mat temp(N, N, arma::fill::zeros);

    for (int h = 0; h < H; ++h) {
        arma::rowvec D = wold.slice(h).row(idx) * S;  // 1 x N
        temp += D.t() * D;                             // N x N, cumulative
        omega += temp;                                  // sum of cumulative sums
    }

    // Leading eigenvector of lower-right (N-1) x (N-1) submatrix
    arma::mat sub = omega.submat(1, 1, N - 1, N - 1);
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    arma::eig_sym(eigenvalues, eigenvectors, sub);

    // eig_sym sorts ascending — last column is the leading eigenvector
    arma::vec h2(N, arma::fill::zeros);
    h2.rows(1, N - 1) = eigenvectors.col(N - 2);

    return h2;
}

//' @export
// [[Rcpp::export]]
arma::vec fuhlig_maxshare(const arma::cube& wold, const arma::mat& S, int idx) {
    return fuhlig_maxshare_cpp(wold, S, idx - 1);  // 1-based -> 0-based
}
