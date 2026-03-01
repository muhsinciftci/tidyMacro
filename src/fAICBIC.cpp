#include <RcppArmadillo.h>
#include "fVAR.h"
#include "fAICBIC.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Internal C++ function (called from other C++ code)
// Calls fVAR_cpp directly - no list overhead in the loop!
AICBICResult fAICBIC_cpp(const arma::mat& y, int pmax, int c,
                         Rcpp::Nullable<arma::mat> exog) {

    const int T = y.n_rows;
    const int N = y.n_cols;

    int M = 0;
    if (exog.isNotNull()) {
        arma::mat exog_mat = Rcpp::as<arma::mat>(exog);
        M = exog_mat.n_cols;
    }

    arma::vec aic(pmax);
    arma::vec bic(pmax);
    arma::vec hq(pmax);

    for (int p = 1; p <= pmax; ++p) {
        // Call fVAR_cpp directly - much faster than calling the R wrapper!
        VARResult var_result = fVAR_cpp(y, p, c, exog);
        const arma::mat& e = var_result.residuals;

        double a = std::log(arma::det((e.t() * e) / T));
        double b = (p * N * N + N + M) / static_cast<double>(T);

        aic(p - 1) = a + 2.0 * b;
        bic(p - 1) = a + std::log(static_cast<double>(T)) * b;
        hq(p - 1)  = a + 2.0 * b * std::log(std::log(static_cast<double>(T)));
    }

    // Return results as struct
    AICBICResult result;
    result.aic = aic.index_min() + 1;  // Implicit conversion from arma::uword to int
    result.bic = bic.index_min() + 1;  // Implicit conversion from arma::uword to int
    result.hq = hq.index_min() + 1;    // Implicit conversion from arma::uword to int
    return result;
}

//' Determine Optimal Lag Order for VAR Model Using Information Criteria
//'
//' Fits VAR models for lag orders 1 through \code{pmax} and returns the lag
//' length that minimises each of three information criteria.
//'
//' @param y A T x N matrix of endogenous variables (T observations, N series).
//' @param pmax An integer specifying the maximum VAR lag order to consider.
//' @param c An integer (0 or 1) indicating whether to include a constant term
//'   (1 = include, 0 = exclude).
//' @param exog An optional T x M matrix of exogenous variables. Default is
//'   NULL.
//'
//' @return A list with three integer elements: aic — optimal lag length by
//'   Akaike IC; bic — optimal lag length by Bayesian IC; hq — optimal lag
//'   length by Hannan-Quinn IC.
//'
//' @details
//' For each candidate lag order p the residual covariance Sigma = (e'e)/T is
//' computed and the three criteria are evaluated as:
//'
//' AIC = log(det(Sigma)) + 2 * (p*N^2 + N + M) / T
//'
//' BIC = log(det(Sigma)) + log(T) * (p*N^2 + N + M) / T
//'
//' HQ  = log(det(Sigma)) + 2 * log(log(T)) * (p*N^2 + N + M) / T
//'
//' where M is the number of exogenous variables.
//'
//' @examples
//' \dontrun{
//' set.seed(123)
//' y      <- matrix(rnorm(200), ncol = 2)
//' result <- fAICBIC(y, pmax = 10, c = 1)
//' result$aic
//' result$bic
//' }
//'
//' @seealso \code{\link{fVAR}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fAICBIC(const arma::mat& y, int pmax, int c,
                   Rcpp::Nullable<arma::mat> exog = R_NilValue) {
    // Call the C++ function
    AICBICResult result = fAICBIC_cpp(y, pmax, c, exog);
    
    // Return results as a list for R
    return Rcpp::List::create(
        Rcpp::Named("aic") = result.aic,
        Rcpp::Named("bic") = result.bic,
        Rcpp::Named("hq")  = result.hq
    );
}
