#include "fcompanionMatrix.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Internal C++ function (called from other C++ code)
CompanionMatrixResult fcompanionMatrix_cpp(arma::mat beta, int c, int p) {

  // N is the number of variables (columns of beta)
  int N = beta.n_cols;

  // Remove intercept row if c == 1
  arma::mat beta_temp;
  if (c == 1) {
    beta_temp = beta.rows(1, beta.n_rows - 1);
  } else {
    beta_temp = beta;
  }

  // Extract only the VAR coefficient rows (exclude exogenous variable
  // coefficients) VAR coefficients are in rows 0 to (N*p - 1), exogenous are in
  // rows (N*p) to end
  int var_coeff_rows = N * p;
  arma::mat beta_coeffs = beta_temp.rows(0, var_coeff_rows - 1);

  // Initialize companion matrix
  int comp_size = N * p;
  arma::mat comp = arma::zeros<arma::mat>(comp_size, comp_size);

  // Fill first N rows with transposed coefficients
  comp.rows(0, N - 1) = beta_coeffs.t();

  // Fill identity block for lagged variables (if p > 1)
  if (p > 1) {
    comp.submat(N, 0, comp_size - 1, N * (p - 1) - 1) =
        arma::eye<arma::mat>(N * (p - 1), N * (p - 1));
  }

  // Return results as struct
  CompanionMatrixResult result;
  result.comp = comp;
  result.N = N;
  return result;
}

//' Compute VAR Companion Matrix
//'
//' @param beta Coefficient matrix. Dimensions: (Np+c+M) x N,
//'   where N is the number of variables, p is the lag order,
//'   c is 1 if intercept included (0 otherwise), and M is the number of
//'   exogenous variables
//' @param c Integer indicator for intercept (1 if intercept included, 0 otherwise)
//' @param p Integer lag order of the VAR model
//'
//' @return A list containing:
//'   \itemize{
//'     \item comp: Companion matrix (Np x Np)
//'     \item N: Number of variables in the VAR system
//'   }
//'
//' @details
//' This function constructs the companion form matrix of a VAR(p) model.
//' The companion matrix representation allows the VAR(p) to be written as a
//' VAR(1) system, which is useful for various computations including impulse
//' response functions and forecasting.
//'
//' When exogenous variables are present, they are excluded from the companion
//' matrix as they do not have a dynamic feedback structure in the VAR system.
//'
//' The companion matrix has the structure:
//' \deqn{
//'   \begin{bmatrix}
//'   A_1 & A_2 & \cdots & A_{p-1} & A_p \\
//'   I_N & 0 & \cdots & 0 & 0 \\
//'   0 & I_N & \cdots & 0 & 0 \\
//'   \vdots & \vdots & \ddots & \vdots & \vdots \\
//'   0 & 0 & \cdots & I_N & 0
//'   \end{bmatrix}
//' }
//'
//' @examples
//' \dontrun{
//' # VAR(2) model with 3 variables, no intercept
//' beta <- matrix(rnorm(18), 6, 3)  # 6 rows (2*3), 3 columns
//' result <- fcompanionMatrix(beta, c = 0, p = 2)
//'
//' # VAR(2) model with 3 variables, with intercept and 2 exogenous variables
//' beta_int <- matrix(rnorm(27), 9, 3)  # 9 rows (2*3+1+2), 3 columns
//' result <- fcompanionMatrix(beta_int, c = 1, p = 2)
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fcompanionMatrix(arma::mat beta, int c, int p) {
  // Call the C++ function
  CompanionMatrixResult result = fcompanionMatrix_cpp(beta, c, p);
  
  // Return results as a list for R
  return Rcpp::List::create(
    Rcpp::Named("comp") = result.comp, 
    Rcpp::Named("N") = result.N
  );
}
