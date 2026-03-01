#include <RcppArmadillo.h>
#include "fcompanionMatrix.h"
#include "fwoldIRF.h"
#include "fVAR.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Internal C++ function (called from other C++ code)
// Accepts VARResult struct directly - no list overhead!
WoldIRFResult fwoldIRF_cpp(const VARResult& var_result, int horizon) {
         
  // Extract VAR elements from struct (much faster than from Rcpp::List)
  const arma::mat& beta = var_result.beta;
  const int c = var_result.c;
  const int p = var_result.p;
  
  // Get companion matrix using _cpp version - no list overhead!
  CompanionMatrixResult comp_result = fcompanionMatrix_cpp(beta, c, p);
  const arma::mat& BigA = comp_result.comp;  // Direct struct access, no conversion needed
  const int N = comp_result.N;
  const int comp_size = BigA.n_rows;  // Implicit conversion from arma::uword to int
  
  // Pre-allocate IRF cube
  arma::cube irfwold(N, N, horizon + 1, arma::fill::none);
  
  // Initialize with identity matrix for h=0
  arma::mat A_power = arma::eye<arma::mat>(comp_size, comp_size);
  
  // Extract first N x N block directly into cube (avoid copy)
  irfwold.slice(0) = A_power.submat(0, 0, N - 1, N - 1);
  
  // Iterative multiplication: A^h = A^(h-1) * A
  // Use in-place operations where possible
  arma::mat A_temp(comp_size, comp_size, arma::fill::none);
  
  for (int h = 1; h <= horizon; ++h) {
    A_temp = A_power * BigA;  // Compute A^h
    A_power.swap(A_temp);      // Swap instead of copy (O(1) operation)
    
    // Direct memory access to cube slice
    irfwold.slice(h) = A_power.submat(0, 0, N - 1, N - 1);
  }
  
  // Return results as struct
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
  var_result.sigma_full = Rcpp::as<arma::mat>(fVAR["sigma_full"]);
  var_result.p = Rcpp::as<int>(fVAR["p"]);
  var_result.c = Rcpp::as<int>(fVAR["c"]);
  var_result.n_exog = Rcpp::as<int>(fVAR["n_exog"]);
  
  // Call the C++ function with the struct
  WoldIRFResult result = fwoldIRF_cpp(var_result, horizon);
  
  // Return the IRF cube for R
  return result.irfwold;
}
