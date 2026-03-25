#include "fbqIRF.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Compute Blanchard-Quah (BQ) Impulse Response Functions
//'
//' @param wold Wold representation cube (N x N x horizon+1), where
//'   \code{wold[,,h]} contains the Wold IRF at horizon h
//' @param K N x N lower triangular Cholesky factor of the long-run
//'   covariance matrix used for BQ identification
//' @param scaling Optional numeric vector of length 2. The first element
//'   specifies the variable index (1-based) used for normalisation, and the
//'   second element specifies the shock size. When omitted no normalisation
//'   is applied.
//'
//' @return A cube (N x N x horizon+1) of long-run identified impulse response
//'   functions.
//'
//' @details
//' Computes structural impulse response functions under the Blanchard-Quah
//' (1989) long-run identification scheme. For each horizon h the structural
//' IRF is computed as:
//' \deqn{IRF_h = \Psi_h \cdot K}
//' where \eqn{\Psi_h} is the Wold representation at horizon h and K is the
//' lower triangular Cholesky factor of the long-run covariance matrix.
//'
//' When \code{scaling} is supplied the entire cube is divided by
//' \eqn{IRF_0(\text{scaling}[1],\, \text{scaling}[1]) \;/\; \text{scaling}[2]},
//' normalising the impact response of the selected variable to
//' \code{scaling[2]}.
//'
//' @references
//' Blanchard, O. J., & Quah, D. (1989). The dynamic effects of aggregate
//' demand and supply disturbances. \emph{American Economic Review}, 79(4),
//' 655--673.
//'
//' @examples
//' \dontrun{
//' # Estimate a VAR and compute Wold IRFs
//' VAR  <- fVAR(y, p = 2, c = 1)
//' wold <- fwoldIRF(VAR, horizon = 20)
//'
//' # Obtain BQ long-run Cholesky factor K (from e.g. fBQ())
//' bqirf <- fbqIRF(wold, K)
//'
//' # With normalisation: unit shock to variable 1
//' bqirf_norm <- fbqIRF(wold, K, scaling = c(1, 1))
//' }
//'
//' @export
// [[Rcpp::export]]
arma::cube fbqIRF(const arma::cube& wold, const arma::mat& K,
                  Rcpp::Nullable<arma::vec> scaling = R_NilValue) {

  int N       = wold.n_rows;
  int horizon = wold.n_slices;

  // Pre-allocate output cube
  arma::cube bqirf(N, N, horizon, arma::fill::none);

  // BQ IRF: for each horizon h, IRF_h = Wold_h * K
  for (int h = 0; h < horizon; h++) {
    bqirf.slice(h) = wold.slice(h) * K;
  }

  // Optional normalisation: divide entire cube by
  //   bqirf(idx, idx, 0) / scaling(1)
  // so that the impact response of the chosen variable equals scaling(1).
  // scaling(0) is 1-based (R convention) -> subtract 1 for C++ indexing.
  if (scaling.isNotNull()) {
    arma::vec s    = Rcpp::as<arma::vec>(scaling);
    int    idx     = static_cast<int>(s(0)) - 1;
    double scale_val = bqirf(idx, idx, 0) / s(1);
    bqirf /= scale_val;
  }

  return bqirf;
}
