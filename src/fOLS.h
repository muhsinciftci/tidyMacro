#ifndef FOLS_H
#define FOLS_H

#include <RcppArmadillo.h>

// Struct to hold OLS estimation results
struct OLSResult {
  arma::mat beta;           // Coefficient estimates
  arma::mat fitted;         // Fitted values
  arma::mat err;            // Residuals
  double r2;                // R-squared
  double r2adj;             // Adjusted R-squared
  double F;                 // F-statistic (overall significance)
  double Frobust;           // Robust F-statistic (White or HAC)
  arma::mat fitted_partial; // Fitted values excluding intercept
};

// Internal C++ function (for use in other C++ code)
// robust: 0 = standard, 1 = White, 2 = Newey-West HAC
// lag: number of lags for HAC (0 = use Newey-West rule of thumb when robust==2)
// mode: 0 = full output, 1 = beta + err only (skips r2, varbhat, F-stat)
OLSResult fOLS_cpp(arma::mat y, arma::mat X, int c, int robust, int lag, int mode);

// R wrapper function (for calling from R)
// lag=0: White robust F; lag>0: Newey-West HAC robust F with that lag
Rcpp::List fOLS(arma::mat y, arma::mat X, int c, int lag);

#endif // FOLS_H
