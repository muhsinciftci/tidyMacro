#include <RcppArmadillo.h>
#include "fevd_chol.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Internal C++ function (called from other C++ code)
FEVDCholResult fevd_chol_cpp(const arma::cube& chol_irf) {
    
    const int N = chol_irf.n_rows;
    const int H = chol_irf.n_slices;
    
    if (chol_irf.n_cols != static_cast<arma::uword>(N)) {
        Rcpp::stop("chol_irf must be N x N x (horizon+1)");
    }
    
    // Pre-allocate FEVD cube: fevd(i, j, h) = share of variance of variable i 
    // at horizon h explained by shock j
    arma::cube fevd(N, N, H, arma::fill::zeros);
    arma::cube cumsum_sq(N, N, H, arma::fill::zeros);
    
    for (int h = 0; h < H; ++h) {
        if (h == 0) {
            cumsum_sq.slice(h) = arma::square(chol_irf.slice(h));
        } else {
            cumsum_sq.slice(h) = cumsum_sq.slice(h - 1) + arma::square(chol_irf.slice(h));
        }
        
        for (int i = 0; i < N; ++i) {
            double total_var = arma::accu(cumsum_sq.slice(h).row(i));
            
            if (total_var > 1e-15) {
                for (int j = 0; j < N; ++j) {
                    fevd(i, j, h) = cumsum_sq(i, j, h) / total_var;
                }
            } else {
                for (int j = 0; j < N; ++j) {
                    fevd(i, j, h) = 1.0 / N;
                }
            }
        }
    }
    
    FEVDCholResult result;
    result.fevd = fevd;
    return result;
}

//' Forecast Error Variance Decomposition for Cholesky-Identified VARs
//'
//' Computes the forecast error variance decomposition (FEVD) showing the
//' proportion of forecast error variance of each variable attributable to
//' each structural shock identified via Cholesky decomposition.
//'
//' @param chol_irf N x N x (horizon+1) cube of Cholesky impulse response functions
//'   from \code{fcholeskyIRF}. Each slice chol_irf[,,h] contains the responses at
//'   horizon h. The horizon is automatically inferred from the cube dimensions.
//' @param shock Integer (1-indexed). If provided, returns only the N x horizon
//'   matrix for that specific shock (matching MATLAB \code{variance_decomp(irf, shock)}).
//'   If 0 (default), returns the full N x N x horizon cube for all shocks.
//'
//' @return A list containing:
//'   \item{fevd}{If \code{shock = 0}: N x N x horizon cube where \code{fevd[i,j,h]}
//'     is the share of forecast error variance of variable i at horizon h explained
//'     by shock j. If \code{shock > 0}: N x horizon matrix for the specified shock.}
//'
//' @details
//' The forecast error variance decomposition measures the proportion of the
//' h-step ahead forecast error variance of variable i that is attributable to
//' each structural shock j. For Cholesky-identified VARs:
//'
//' \deqn{FEVD(i,j,h) = \frac{\sum_{s=0}^{h} [\Psi_s P]_{ij}^2}{\sum_{k=1}^{N} \sum_{s=0}^{h} [\Psi_s P]_{ik}^2}}
//'
//' @examples
//' \dontrun{
//' # Estimate VAR
//' var_result <- fVAR(y, p = 2, c = 1)
//' wold <- fwoldIRF(var_result, horizon = 20)
//' S <- t(chol(var_result$sigma_full))
//' chol_irf <- fcholeskyIRF(wold, S)
//'
//' # Full decomposition (all shocks) — returns 3D cube
//' vardec <- fevd_chol(chol_irf)
//'
//' # Single shock decomposition — returns N x horizon matrix
//' shock <- match("UNCERT", colnames(y))
//' vardec <- fevd_chol(chol_irf, shock = shock)
//' }
//'
//' @seealso \code{\link{fcholeskyIRF}}, \code{\link{fwoldIRF}},
//'   \code{\link{fplot_vardec}}, \code{\link{fevd_iv}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fevd_chol(const arma::cube& chol_irf, int shock = 0) {
    
    FEVDCholResult result = fevd_chol_cpp(chol_irf);
    
    if (shock == 0) {
        // Return full 3D cube (all shocks)
        return Rcpp::List::create(
            Rcpp::Named("fevd") = result.fevd
        );
    } else {
        // Validate shock index
        const int N = chol_irf.n_rows;
        if (shock < 1 || shock > N) {
            Rcpp::stop("shock must be between 1 and N (%d)", N);
        }
        // Extract N x horizon matrix for the specified shock (0-indexed: shock-1)
        arma::mat fevd_shock = result.fevd.col(shock - 1);  // N x H matrix
        return Rcpp::List::create(
            Rcpp::Named("fevd") = fevd_shock
        );
    }
}
