#include <RcppArmadillo.h>
#include "fCheckRestrictions.h"
// [[Rcpp::depends(RcppArmadillo)]]

int fCheckRestrictions_cpp(const arma::cube&   irf,
                           int                 shock_idx,
                           const arma::rowvec& restr,
                           const arma::uvec&   hor_vec) {

    const arma::uword N = irf.n_rows;

    if (irf.n_cols != N) {
        Rcpp::stop("irf must be an N x N x H cube.");
    }
    if (shock_idx < 0 || static_cast<arma::uword>(shock_idx) >= N) {
        Rcpp::stop("shock index is out of range.");
    }
    if (restr.n_elem != N) {
        Rcpp::stop("restr must have length N.");
    }
    if (hor_vec.n_elem == 0u) {
        return 1;
    }

    arma::uword n_restricted = 0u;
    for (arma::uword i = 0; i < N; ++i) {
        n_restricted += (restr(i) != 0.0);
    }
    if (n_restricted == 0u) {
        return 1;
    }

    bool all_horizons_pass = true;
    bool all_horizons_fail = true;

    for (arma::uword r = 0; r < hor_vec.n_elem; ++r) {
        const arma::uword h = hor_vec(r);
        if (h >= irf.n_slices) {
            Rcpp::stop("hor_vec contains a horizon outside the IRF cube.");
        }

        const double* col = irf.slice(h).colptr(static_cast<arma::uword>(shock_idx));
        arma::uword matches = 0u;

        for (arma::uword i = 0; i < N; ++i) {
            const double sign_restr = restr(i);
            if (sign_restr == 0.0) {
                continue;
            }

            const double val = col[i];
            matches += (sign_restr > 0.0) ? (val > 0.0) : (val < 0.0);
        }

        all_horizons_pass = all_horizons_pass && (matches == n_restricted);
        all_horizons_fail = all_horizons_fail && (matches == 0u);

        if (!all_horizons_pass && !all_horizons_fail) {
            return -1;
        }
    }

    if (all_horizons_pass) {
        return 1;
    }
    if (all_horizons_fail) {
        return 0;
    }
    return -1;
}

//' Check Sign Restrictions for One Shock
//'
//' For a candidate orthonormal rotation Q applied to the Cholesky IRFs,
//' tests whether the impulse responses for a given shock satisfy the imposed
//' sign restrictions at all specified horizons.
//'
//' @param irf N x N x H cube of structural IRFs (e.g. \code{cholirf \%*\% Q}
//'   slice-by-slice).
//' @param shock Integer (1-indexed). Which column of \code{irf} to examine.
//' @param restr Numeric row vector of length N. Sign restrictions per variable:
//'   \code{+1} (must be positive), \code{-1} (must be negative), \code{0} (unrestricted).
//' @param hor_vec Integer vector of horizons to check (1-indexed, where 1 = impact).
//'
//' @return Integer scalar:
//'   \itemize{
//'     \item \code{+1} all restrictions pass at every checked horizon
//'     \item \code{ 0} all restrictions fail at every horizon; flipping the sign of
//'       this shock column would fix them all
//'     \item \code{-1} mixed result; this draw of Q must be discarded
//'   }
//'
//' @details
//' Implements the restriction-checking step of the algorithm in Rubio-Ramirez,
//' Waggoner & Zha (2010). Restrictions with \code{restr[i] = 0} are treated as
//' "generic" (no restriction) and are ignored when counting satisfied constraints.
//'
//' @examples
//' \dontrun{
//' # After computing cholirf and drawing Q:
//' Q   <- fGenerateQ(3)
//' irf <- array(NA, dim = c(3,3,21))
//' for (h in 1:21) irf[,,h] <- cholirf[,,h] %*% Q
//'
//' # Check monetary policy shock (col 3): FFR up, GDP and PCE down at impact
//' restr   <- c(-1, -1, 1)
//' hor_vec <- c(1L)
//' fCheckRestrictions(irf, shock = 3, restr = restr, hor_vec = hor_vec)
//' }
//'
//' @seealso \code{\link{fGenerateQ}}, \code{\link{fSignRestrictions}}
//'
//' @export
// [[Rcpp::export]]
int fCheckRestrictions(const arma::cube&   irf,
                       int                 shock,
                       const arma::rowvec& restr,
                       const arma::uvec&   hor_vec) {
    if (shock < 1 || static_cast<arma::uword>(shock) > irf.n_cols) {
        Rcpp::stop("shock must be between 1 and N.");
    }
    if (arma::any(hor_vec == 0u)) {
        Rcpp::stop("hor_vec must use 1-indexed horizons.");
    }

    return fCheckRestrictions_cpp(irf, shock - 1, restr, hor_vec - 1u);
}
