#include "fCheckStability.h"
#include "fCompanionMatrix.h"
#include "fVAR.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Internal C++ function (called from other C++ code)
StabilityResult fCheckStability_cpp(const arma::mat& beta, int c, int p) {

  // Stability of a VAR(p) is a property of the (Np x Np) companion matrix F,
  // not of A_1 alone: for p > 1 the eigenvalues of A_1 are neither necessary
  // nor sufficient for stationarity.
  CompanionMatrixResult comp_res = fCompanionMatrix_cpp(beta, c, p);

  // Eigenvalues of the (generally non-symmetric) companion matrix
  arma::cx_vec eigval = arma::eig_gen(comp_res.comp);

  StabilityResult result;
  // Moduli |lambda_i| (companion eigenvalues can be complex)
  result.max_eig = arma::abs(eigval).max();
  // A single eigenvalue with |lambda_i| >= 1 is enough for non-stationarity
  result.stable  = (result.max_eig < 1.0);
  return result;
}

// Internal C++ function (VAR result overload)
StabilityResult fCheckStability_cpp(const VARResult& var_result) {
  return fCheckStability_cpp(var_result.beta, var_result.c, var_result.p);
}

//' Check VAR Stability
//'
//' @param var_model A list containing VAR estimation results with elements:
//'   \itemize{
//'     \item beta: Coefficient matrix
//'     \item c: Integer indicator for intercept (1 if intercept, 0 otherwise)
//'     \item p: Integer lag order
//'   }
//'
//' @return Called for its side-effect of printing to the console.
//'   Returns \code{NULL} invisibly.
//'
//' @export
// [[Rcpp::export]]
void fCheckStability(const Rcpp::List& var_model) {

  // Extract only the elements needed for the companion matrix
  VARResult var_result;
  var_result.beta = Rcpp::as<arma::mat>(var_model["beta"]);
  var_result.p    = Rcpp::as<int>(var_model["p"]);
  var_result.c    = Rcpp::as<int>(var_model["c"]);

  StabilityResult result = fCheckStability_cpp(var_result);

  Rcpp::Rcout << "Maximum eigenvalue (modulus): " << result.max_eig
              << std::endl;
  if (result.stable) {
    Rcpp::Rcout << "The VAR is stable (all |lambda_i| < 1)." << std::endl;
  } else {
    Rcpp::Rcout << "The VAR is unstable (at least one |lambda_i| >= 1)."
                << std::endl;
  }
}
