// [[Rcpp::depends(RcppArmadillo)]]
#include "fuhlig_maxshare.h"
#include <RcppArmadillo.h>

// Uhlig max-share structural IRF (point estimate companion to fbootstrapUhlig).
// Uses fuhlig_maxshare_cpp to find the FEV-maximising h2, then applies the
// same sign convention as the bootstrap: first variable's last-horizon
// non-structural response is non-negative.

//' @export
// [[Rcpp::export]]
arma::mat fuhligIRF(const arma::cube& wold, const arma::mat& S, int idx) {
    const int H = static_cast<int>(wold.n_slices);
    const int N = static_cast<int>(wold.n_rows);

    arma::vec h2 = fuhlig_maxshare_cpp(wold, S, idx - 1);  // 1-based -> 0-based

    // Precompute impact = S * h2 once (used for both sign check and IRF loop)
    arma::vec impact = S * h2;

    // Sign: first variable's last-horizon non-structural response non-negative
    if (arma::dot(wold.slice(H - 1).row(0), impact) < 0.0) impact = -impact;

    // Structural IRFs: N x H matrix
    arma::mat irf(N, H, arma::fill::none);
    for (int h = 0; h < H; ++h)
        irf.col(h) = wold.slice(h) * impact;

    return irf;
}
