#include <RcppArmadillo.h>
#include "fevd_iv.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Internal C++ function (called from other C++ code)
FEVDIVResult fevd_iv_cpp(arma::vec s, arma::mat S, arma::cube wold, int N, int hor,
                         arma::mat sigma, arma::mat u, int T1, int p) {
    arma::vec scaler = s / arma::norm(arma::solve(S, s), 2);
    
    arma::mat ivirf_scaled = arma::zeros(N, hor + 1);
    for (int h = 0; h < hor + 1; h++) {
        ivirf_scaled.col(h) = wold.slice(h) * scaler;
    }
    
    arma::rowvec oil_news_unit_var = scaler.t() * arma::solve(sigma, u.t());
    
    arma::mat check_uv = oil_news_unit_var * oil_news_unit_var.t() / (T1 - 1 - p - N * p);
    
    // Get the total variation of the forecast errors
    arma::mat temp = arma::zeros(N, N);
    arma::mat denom = arma::zeros(N, hor + 1);
    for (int h = 0; h < hor + 1; h++) {
        temp = temp + wold.slice(h) * sigma * wold.slice(h).t();
        denom.col(h) = arma::diagvec(temp);  // we only want the variances, not the covariances
    }
    
    // Get the variation in forecast errors due to a specific shock
    arma::mat irf_sq = arma::square(ivirf_scaled);
    arma::mat num = arma::cumsum(irf_sq, 1);  // cumsum along rows (dimension 1)
    
    // Compute the share of the total
    arma::mat fevd_iv_result = num / denom;
    
    // Return results as struct
    FEVDIVResult result;
    result.scaler = scaler;
    result.ivirf_scaled = ivirf_scaled;
    result.oil_news_unit_var = oil_news_unit_var;
    result.check_uv = check_uv;
    result.denom = denom;
    result.fevd_iv = fevd_iv_result;
    return result;
}

//' Forecast Error Variance Decomposition for IV-Identified Shocks
//'
//' Computes the forecast error variance decomposition (FEVD) for instrumental
//' variable identified structural shocks. The function scales shocks to unit
//' variance and calculates the contribution of a specific shock to the forecast
//' error variance of each variable over time.
//'
//' @param s N x 1 vector of structural impact coefficients
//' @param S N x N lower triangular Cholesky factor of the residual covariance matrix
//' @param wold N x N x (hor+1) cube of Wold impulse response functions
//' @param N Integer number of variables in the VAR system
//' @param hor Integer maximum forecast horizon
//' @param sigma N x N residual covariance matrix
//' @param u (T-p) x N matrix of VAR residuals
//' @param T1 Integer effective sample size (T-p)
//' @param p Integer lag order of the VAR model
//'
//' @return A list containing:
//'   - scaler: N x 1 vector of scaling factors for unit variance normalization
//'   - ivirf_scaled: N x (hor+1) matrix of unit variance impulse responses
//'   - oil_news_unit_var: 1 x T1 vector of unit variance shock series
//'   - check_uv: Scalar variance check (should be close to 1)
//'   - denom: N x (hor+1) matrix of total forecast error variances
//'   - fevd_iv: N x (hor+1) matrix of FEVD shares (between 0 and 1)
//'
//' @details
//' The FEVD measures the proportion of the h-step ahead forecast error variance
//' of each variable that is attributable to the identified structural shock.
//'
//' The function performs the following steps:
//' 1. Normalizes the structural shock to have unit variance
//' 2. Computes unit variance impulse responses
//' 3. Recovers the unit variance shock series
//' 4. Computes total forecast error variance at each horizon
//' 5. Computes variance due to identified shock
//' 6. Calculates FEVD as the ratio
//'
//' The FEVD values are proportions between 0 and 1, where higher values indicate
//' that the identified shock explains a larger share of the forecast error variance.
//'
//' @examples
//' \dontrun{
//' # Estimate VAR and compute Wold IRF
//' var_result <- fVAR(y, p = 2, c = 1)
//' wold <- fwoldIRF(var_result, horizon = 20)
//' 
//' # IV identification
//' s <- c(1, 0.5, 0.3)
//' sigma <- var_result$sigma_full
//' S <- t(chol(sigma))
//' 
//' # Compute FEVD
//' result <- fevd_iv(s, S, wold, N = 3, hor = 20,
//'                   sigma = sigma, u = var_result$residuals,
//'                   T1 = nrow(var_result$residuals), p = 2)
//' 
//' # Plot FEVD for first variable
//' plot(result$fevd_iv[1, ], type = "l", ylim = c(0, 1),
//'      main = "FEVD: Share of Variance Explained by Shock",
//'      xlab = "Horizon", ylab = "Share")
//' }
//'
//' @seealso \code{\link{fwoldIRF}} for computing Wold impulse responses,
//'   \code{\link{fVAR}} for VAR estimation
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fevd_iv(arma::vec s, arma::mat S, arma::cube wold, int N, int hor,
                   arma::mat sigma, arma::mat u, int T1, int p) {
    // Call the C++ function
    FEVDIVResult result = fevd_iv_cpp(s, S, wold, N, hor, sigma, u, T1, p);
    
    // Return results as a list for R
    return Rcpp::List::create(
        Rcpp::Named("scaler") = result.scaler,
        Rcpp::Named("ivirf_scaled") = result.ivirf_scaled,
        Rcpp::Named("oil_news_unit_var") = result.oil_news_unit_var,
        Rcpp::Named("check_uv") = result.check_uv,
        Rcpp::Named("denom") = result.denom,
        Rcpp::Named("fevd_iv") = result.fevd_iv
    );
}
