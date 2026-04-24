#include <RcppArmadillo.h>
#include "fwoldIRF.h"
#include "fVAR.h"
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

// Internal C++ function (called from other C++ code)
// Uses block recursion Phi_h = sum_{l=1..min(h,p)} Phi_{h-l} * A_l
// O(N^3 * p) per step vs O((N*p)^3) for companion matrix powers.
WoldIRFResult fwoldIRF_cpp(const VARResult& var_result, int horizon) {
  const arma::mat& beta = var_result.beta;
  const int c = var_result.c;
  const int p = var_result.p;
  const int N = static_cast<int>(beta.n_cols);

  // Pre-extract A_l coefficient matrices (each N x N, 0-indexed: A[0]=A_1...)
  std::vector<arma::mat> A(p);
  for (int l = 0; l < p; ++l)
    A[l] = beta.rows(c + l * N, c + (l + 1) * N - 1).t();

  arma::cube irfwold(N, N, horizon + 1, arma::fill::none);
  irfwold.slice(0) = arma::eye<arma::mat>(N, N);

  for (int h = 1; h <= horizon; ++h) {
    arma::mat Phi_h(N, N, arma::fill::zeros);
    const int lmax = std::min(h, p);
    for (int l = 1; l <= lmax; ++l)
      Phi_h += irfwold.slice(h - l) * A[l - 1];
    irfwold.slice(h) = Phi_h;
  }

  WoldIRFResult result;
  result.irfwold = irfwold;
  return result;
}

//' Compute Wold Impulse Response Functions for VAR Model
//'
//' @param fVAR A list containing VAR estimation results with elements:
//'   \itemize{
//'     \item beta: Coefficient matrix
//'     \item c: Integer indicator for intercept (1 if intercept, 0 otherwise)
//'     \item p: Integer lag order
//'     \item n_exog: Number of exogenous variables (optional)
//'   }
//' @param horizon Integer, the IRF horizon (number of periods ahead)
//'
//' @return A 3D array (cube) of Wold impulse response functions with dimensions
//'   N x N x (horizon+1), where irfwold[,,h] contains the response at horizon h
//'
//' @details
//' This function computes the Wold (reduced-form) impulse response functions
//' for a VAR(p) model with optional exogenous variables by calculating powers 
//' of the companion matrix. The Wold IRF shows the effect of a one-unit shock 
//' to each variable on all variables in the system over time, without imposing 
//' any identifying restrictions.
//'
//' Exogenous variables do not contribute to the dynamic propagation of shocks
//' and are therefore excluded from the companion matrix representation.
//'
//' The IRFs are computed as:
//' \deqn{\Phi_h = A^h[1:N, 1:N]}
//' where A is the companion matrix and h is the horizon.
//'
//' @examples
//' \dontrun{
//' # VAR(2) model with 3 variables and 1 exogenous variable
//' y <- matrix(rnorm(200), ncol = 2)
//' exog <- matrix(rnorm(100), ncol = 1)
//' VAR <- fVAR(y, p = 2, c = 1, exog = exog)
//' irf <- fwoldIRF(VAR, horizon = 10)
//' # irf is a 2x2x11 array
//' 
//' # Plot impulse response
//' plot(irf[1, 2, ], type = "l", 
//'      main = "Response of Variable 2 to Shock in Variable 1")
//' }
//'
//' @export
// [[Rcpp::export]]
arma::cube fwoldIRF(const Rcpp::List& fVAR, int horizon) {
  // Extract elements from R list and construct VARResult struct
  VARResult var_result;
  var_result.beta = Rcpp::as<arma::mat>(fVAR["beta"]);
  var_result.residuals = Rcpp::as<arma::mat>(fVAR["residuals"]);
  var_result.sigma = Rcpp::as<arma::mat>(fVAR["sigma"]);
  var_result.p = Rcpp::as<int>(fVAR["p"]);
  var_result.c = Rcpp::as<int>(fVAR["c"]);
  var_result.n_exog = Rcpp::as<int>(fVAR["n_exog"]);
  
  // Call the C++ function with the struct
  WoldIRFResult result = fwoldIRF_cpp(var_result, horizon);
  
  // Return the IRF cube for R
  return result.irfwold;
}
