#ifndef FMSW_H
#define FMSW_H

#include <RcppArmadillo.h>

// Results of Montiel-Olea, Stock, Watson (2021) weak-IV robust inference
struct MSWResult {
    arma::mat IRF;        // n x (hor+1) plug-in (delta method) IRF
    arma::mat IRFstderr;  // n x (hor+1) delta-method standard errors
    arma::mat Dmlbound;   // n x (hor+1) delta-method lower bound
    arma::mat Dmubound;   // n x (hor+1) delta-method upper bound
    arma::mat MSWlbound;  // n x (hor+1) Anderson-Rubin lower bound
    arma::mat MSWubound;  // n x (hor+1) Anderson-Rubin upper bound
    double Waldstat;      // Wald stat for instrument relevance (= robust F-stat)
    double Fstat;         // Standard (non-robust) F-stat (informational)
};

// Internal C++ function
// AL       : n x n*p VAR coefficients (no constant), as [A1 | A2 | ... | Ap]
// Sigma_sub: n x n sample covariance of proxy-sample residuals (denom T1-1)
// eta      : n x T1 residuals (transposed) for the proxy sample
// X        : T1 x (m + n*p) design matrix (intercept + lags), proxy sample
// Z        : T1 x k instrument matrix (proxy sample)
// nvar     : 1-based normalization variable index
MSWResult fMSW_cpp(const arma::mat& AL,
                   const arma::mat& Sigma_sub,
                   const arma::mat& eta,
                   const arma::mat& X,
                   const arma::mat& Z,
                   int p, int hor,
                   int nvar,
                   double scale,
                   double confidence,
                   int NWlags);

// R wrapper
Rcpp::List fMSW(const Rcpp::List& var_result,
                const arma::mat& Z,
                const arma::mat& finaldata,
                const arma::ivec& adjustu,
                int hor,
                int nvar,
                double scale,
                double confidence,
                int NWlags);

#endif // FMSW_H
