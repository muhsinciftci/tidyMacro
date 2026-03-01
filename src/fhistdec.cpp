#include <RcppArmadillo.h>
#include "fhistdec.h"
#include "fVAR.h"
#include "fwoldIRF.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Internal C++ function
HistDecResult fhistdec_cpp(const arma::mat& y,
                           const VARResult& var_result,
                           const arma::mat& K,
                           int series) {

    const int T       = y.n_rows;
    const int N       = y.n_cols;
    const int p       = var_result.p;
    const int T_eff   = T - p;          // effective sample length after lags

    // 0-indexed series position
    const int s_idx = series - 1;

    // ------------------------------------------------------------------
    // 1. Structural shocks: e_struct = residuals * K^{-T}
    //    MATLAB: struct_shock = (K \ residuals')' 
    //    i.e. solve K * X = residuals'  =>  X = K^{-1} * residuals'
    //    transposed: struct_shock = residuals * K^{-T}
    // ------------------------------------------------------------------
    const arma::mat& residuals = var_result.residuals;   // (T-p) x N
    arma::mat struct_shock = arma::solve(arma::trimatl(K), residuals.t()).t();
    // result: (T-p) x N

    // ------------------------------------------------------------------
    // 2. Demeaned series for the variable of interest
    //    MATLAB: ystar = y(p+1:end, series) - mean(y(p+1:end, series))
    // ------------------------------------------------------------------
    arma::vec y_eff = y.submat(p, s_idx, T - 1, s_idx);  // (T-p) x 1
    arma::vec ystar = y_eff - arma::mean(y_eff);

    // ------------------------------------------------------------------
    // 3. Structural MA coefficients: ma_coeff(:,:,h) = wold(:,:,h) * K
    //    Compute Wold IRFs over the full effective sample length
    // ------------------------------------------------------------------
    WoldIRFResult wold_result = fwoldIRF_cpp(var_result, T_eff - 1);
    // wold_result.irfwold: N x N x T_eff cube (horizons 0 ... T_eff-1)

    // Pre-multiply each slice by K to get structural MA coefficients
    // ma_coeff: N x N x T_eff
    arma::cube ma_coeff(N, N, T_eff, arma::fill::none);
    for (int h = 0; h < T_eff; ++h) {
        ma_coeff.slice(h) = wold_result.irfwold.slice(h) * K;
    }

    // ------------------------------------------------------------------
    // 4. Historical decomposition
    //    histdec(t, j) = sum_{h=0}^{t} ma_coeff(s_idx, j, h) * struct_shock(t-h, j)
    //    MATLAB: squeeze(ma_coeff(series, j, 1:t))' * flipud(struct_shock(1:t, j))
    // ------------------------------------------------------------------
    arma::mat histdec(T_eff, N, arma::fill::zeros);

    for (int t = 0; t < T_eff; ++t) {
        for (int j = 0; j < N; ++j) {
            double contrib = 0.0;
            for (int h = 0; h <= t; ++h) {
                // ma_coeff(s_idx, j, h) * struct_shock(t - h, j)
                contrib += ma_coeff(s_idx, j, h) * struct_shock(t - h, j);
            }
            histdec(t, j) = contrib;
        }
    }

    HistDecResult result;
    result.histdec = histdec;
    result.ystar   = ystar;
    return result;
}

//' Historical Decomposition of a VAR Variable
//'
//' Decomposes a chosen variable's realisation into contributions from each
//' structural shock, identified via a lower-triangular Cholesky factor K.
//' Mirrors the MATLAB \code{hist_decmp(y, beta, residuals, c, p, K, series)}.
//'
//' @param y TxN numeric matrix of original (undemeaned) data.
//' @param fVAR List returned by \code{fVAR()}, containing at minimum
//'   \code{beta}, \code{residuals}, \code{sigma_full}, \code{p}, \code{c},
//'   and \code{n_exog}.
//' @param K NxN lower-triangular Cholesky factor of the residual covariance
//'   matrix (i.e. \code{t(chol(sigma_full))}).
//' @param series Integer (1-indexed) selecting which variable to decompose.
//'
//' @return A list with two elements:
//'   \item{histdec}{(T-p) x N numeric matrix. Column \code{j} is the
//'     cumulative contribution of structural shock \code{j} to the chosen
//'     variable at each point in time.}
//'   \item{ystar}{(T-p) numeric vector of the demeaned realisation of the
//'     chosen variable (benchmark series for the plot).}
//'
//' @details
//' The structural shocks are recovered as
//' \deqn{\varepsilon_t = K^{-1} u_t}
//' where \eqn{u_t} are the reduced-form residuals. The structural MA
//' representation is built by multiplying each Wold matrix by K:
//' \deqn{\Theta_h = \Psi_h K}
//' The contribution of shock \eqn{j} to variable \eqn{i} at time \eqn{t} is
//' then the inner product of \eqn{\Theta_{0:t}[i,j]} with the time-reversed
//' structural shocks \eqn{\varepsilon_{t:-1:0,j}}:
//' \deqn{HD(t,j) = \sum_{h=0}^{t} \Theta_h[i,j]\, \varepsilon_{t-h,j}}
//'
//' @examples
//' \dontrun{
//' var_result <- fVAR(y, p = 12, c = 1)
//' K <- t(chol(var_result$sigma_full))
//' series <- match("UNCERT", colnames(y))
//'
//' hd <- fhistdec(y, var_result, K, series)
//'
//' # hd$histdec is (T-p) x N; hd$ystar is the demeaned series
//' matplot(hd$histdec, type = "l")
//' lines(hd$ystar, lwd = 2)
//' }
//'
//' @seealso \code{\link{fVAR}}, \code{\link{fcholeskyIRF}},
//'   \code{\link{fwoldIRF}}, \code{\link{plothistdec}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fhistdec(const arma::mat& y,
                    const Rcpp::List& fVAR,
                    const arma::mat& K,
                    int series) {

    // Reconstruct VARResult struct from R list
    VARResult var_result;
    var_result.beta       = Rcpp::as<arma::mat>(fVAR["beta"]);
    var_result.residuals  = Rcpp::as<arma::mat>(fVAR["residuals"]);
    var_result.sigma_full = Rcpp::as<arma::mat>(fVAR["sigma_full"]);
    var_result.p          = Rcpp::as<int>(fVAR["p"]);
    var_result.c          = Rcpp::as<int>(fVAR["c"]);
    var_result.n_exog     = Rcpp::as<int>(fVAR["n_exog"]);

    const int N = y.n_cols;

    // Input validation
    if (K.n_rows != static_cast<arma::uword>(N) ||
        K.n_cols != static_cast<arma::uword>(N)) {
        Rcpp::stop("K must be an N x N matrix (N = number of variables)");
    }
    if (series < 1 || series > N) {
        Rcpp::stop("series must be between 1 and N (%d)", N);
    }

    HistDecResult result = fhistdec_cpp(y, var_result, K, series);

    return Rcpp::List::create(
        Rcpp::Named("histdec") = result.histdec,
        Rcpp::Named("ystar")   = result.ystar
    );
}
