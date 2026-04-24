// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Closed-form LR-max structural IRF (point estimate companion to fbootstrapMax).
// Maximises wold.slice(H-1).row(var_idx-1) * S * h2 subject to h2(0)=0, ||h2||=1.
// Sign: ensures first variable's last-horizon non-structural response is non-negative.

//' @export
// [[Rcpp::export]]
arma::mat fmaxIRF(const arma::cube& wold, const arma::mat& S, int var_idx) {
    const int H = static_cast<int>(wold.n_slices);
    const int N = static_cast<int>(wold.n_rows);

    // Optimal h2: normalised last column of (wold(H-1, var_idx-1, :) * S).cols(1..N-1)
    arma::rowvec M_full = wold.slice(H - 1).row(var_idx - 1) * S;
    arma::rowvec M_sub  = M_full.cols(1, N - 1);
    double M_norm = arma::norm(M_sub);

    arma::vec h2(N, arma::fill::zeros);
    if (M_norm > 1e-14)
        h2.rows(1, N - 1) = M_sub.t() / M_norm;

    // Precompute impact = S * h2 once (used for both sign check and IRF loop)
    arma::vec impact = S * h2;

    // Sign: first variable's last-horizon non-structural response non-negative
    // Use a row-dot product with the precomputed impact instead of forming
    // the full (wold.slice(H-1) * S) matrix.
    if (arma::dot(wold.slice(H - 1).row(0), impact) < 0.0) impact = -impact;

    // Structural IRFs: N x H matrix
    arma::mat irf(N, H, arma::fill::none);
    for (int h = 0; h < H; ++h)
        irf.col(h) = wold.slice(h) * impact;

    return irf;
}
