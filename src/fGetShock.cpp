#include <RcppArmadillo.h>
#include "fGetShock.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Internal C++ function
// Stock-Watson (2018) formula, two normalisation modes:
//
//   normalize = "unit":
//     shock = shockSize * (U * Sigma^{-1} * s) / (s' * Sigma^{-1} * s)
//     s must be unit-normalised (first element = 1); shockSize rescales
//     (e.g. shockSize = 10 for a 10% oil price shock)
//
//   normalize = "sd":
//     b1 = s / sqrt(s' * Sigma^{-1} * s)  →  shock has unit variance
//     shock = U * Sigma^{-1} * s / sqrt(s' * Sigma^{-1} * s)
//     shockSize ignored
//
// residuals : T x N matrix of VAR residuals
// sigma_full: N x N covariance matrix (full sample)
// s         : N x 1 structural impact vector
arma::vec fGetShock_cpp(const arma::mat& residuals, const arma::mat& sigma_full,
                        const arma::vec& s, double shockSize,
                        std::string normalize) {
    arma::vec sig_inv_s = arma::solve(sigma_full, s);
    arma::vec num       = residuals * sig_inv_s;

    double quad = arma::as_scalar(s.t() * sig_inv_s);   // s' * Sigma^{-1} * s

    if (normalize == "unit") {
        return shockSize * num / quad;
    }
    // "sd": divide by sqrt(s' * Sigma^{-1} * s) so shock has unit variance
    return num / std::sqrt(quad);
}

//' Recover Structural Shock Series (Stock-Watson 2018)
//'
//' @param residuals A T x N matrix of VAR residuals.
//' @param sigma_full An N x N residual covariance matrix (full sample).
//' @param s An N x 1 structural impact vector.
//' @param shockSize Scalar normalisation applied in \code{"unit"} mode (default 1).
//'   Set to 10 to express shocks in "10\% impact" units, consistent with IRF
//'   scaling. Ignored when \code{normalize = "sd"}.
//' @param normalize Character string selecting the normalisation mode.
//'   \itemize{
//'     \item \code{"unit"} (default): divides by \eqn{s' \Sigma^{-1} s} and
//'       scales by \code{shockSize}. Use when \code{s} is unit-normalised
//'       (first element = 1).
//'     \item \code{"sd"}: returns \eqn{U \Sigma^{-1} s} directly. Use when
//'       \code{s} is already sd-normalised.
//'   }
//'
//' @return A T x 1 numeric vector of recovered structural shocks.
//'
//' @export
// [[Rcpp::export]]
arma::vec fGetShock(const arma::mat& residuals, const arma::mat& sigma_full,
                    const arma::vec& s, double shockSize = 1.0,
                    std::string normalize = "unit") {
    return fGetShock_cpp(residuals, sigma_full, s, shockSize, normalize);
}
