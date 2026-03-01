#include <RcppArmadillo.h>
#include "flagmakerMatrix.h"
#include "fVARX.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// Internal C++ function (called from other C++ code)
VARXResult fVARX_cpp(const arma::mat& y, const arma::mat& ex, int p, int c) {

    const int T = y.n_rows;
    const arma::mat yfinal = y.rows(p, T - 1);

    arma::mat X;
    if (c == 1) {
        X = arma::ones<arma::mat>(T - p, 1);
        X = arma::join_rows(X, flagmakerMatrix(y, p));
        X = arma::join_rows(X, ex.rows(p, T - 1));
        if (p > 0) {
            X = arma::join_rows(X, flagmakerMatrix(ex, p));
        }
    } else {
        X = flagmakerMatrix(y, p);
        X = arma::join_rows(X, ex.rows(p, T - 1));
        if (p > 0) {
            X = arma::join_rows(X, flagmakerMatrix(ex, p));
        }
    }

    const arma::mat beta = arma::solve(X.t() * X, X.t() * yfinal,
                                       arma::solve_opts::fast);

    // Return results as struct
    VARXResult result;
    result.beta = beta;
    result.residuals = yfinal - X * beta;
    return result;
}

//' Vector Autoregression with Exogenous Variables (VARX) Model Estimation
//'
//' Estimates a VARX(p) model where both endogenous and exogenous variables
//' are lagged using ordinary least squares (OLS). Optimised for maximum
//' speed and memory efficiency.
//'
//' @param y A numeric matrix of endogenous time series data with T rows
//'   (time periods) and N columns (variables).
//' @param ex A numeric matrix of exogenous time series data with T rows
//'   (time periods) and M columns (variables).
//' @param p An integer specifying the number of lags for both endogenous
//'   and exogenous variables.
//' @param c An integer indicator for including a constant term (1 = include
//'   intercept, 0 = no intercept).
//'
//' @return A list with two elements: beta — coefficient matrix of dimensions
//'   (Np + M + Mp + c) x N; residuals — (T-p) x N residual matrix.
//'
//' @details
//' The VARX(p) model is:
//' \deqn{Y_t = c + A_1 Y_{t-1} + \cdots + A_p Y_{t-p}
//'            + B_0 X_t + B_1 X_{t-1} + \cdots + B_p X_{t-p} + e_t}
//' Both endogenous and exogenous variables are lagged p periods.
//' Exogenous variables also enter contemporaneously at time t.
//'
//' @examples
//' \dontrun{
//' y  <- matrix(rnorm(200), ncol = 2)
//' ex <- matrix(rnorm(100), ncol = 1)
//' result    <- fVARX(y, ex, p = 2, c = 1)
//' beta      <- result$beta
//' residuals <- result$residuals
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fVARX(const arma::mat& y, const arma::mat& ex, int p, int c) {
    // Call the C++ function
    VARXResult result = fVARX_cpp(y, ex, p, c);
    
    // Return results as a list for R
    return Rcpp::List::create(
        Rcpp::Named("beta")      = result.beta,
        Rcpp::Named("residuals") = result.residuals
    );
}
