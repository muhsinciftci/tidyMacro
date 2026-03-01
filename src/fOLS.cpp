#include <RcppArmadillo.h>
#include "fOLS.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Internal C++ function (called from other C++ code)
OLSResult fOLS_cpp(arma::mat y, arma::mat X, int c) {
         
  // Get dimensions
  int T = X.n_rows;
  int N = X.n_cols;
  
  // Add intercept if c == 1
  arma::mat X_reg;
  if (c == 1) {
    X_reg = arma::join_rows(arma::ones<arma::vec>(T), X);
  } else {
    X_reg = X;
  }
  
  // Compute OLS estimates: beta = (X'X)^(-1) X'y
  // beta will be (N+c) x K matrix
  arma::mat beta = arma::solve(X_reg.t() * X_reg, X_reg.t() * y);
  
  // Compute fitted values
  arma::mat fitted = X_reg * beta;
  
  // Compute fitted values excluding intercept (if applicable)
  arma::mat fitted_partial;
  if (c == 1) {
    fitted_partial = X_reg.cols(1, N) * beta.rows(1, N);
  } else {
    fitted_partial = fitted;
  }
  
  // Compute residuals
  arma::mat err = y - fitted;
  
  // Compute overall R-squared for the system
  arma::vec err_vec = arma::vectorise(err);
  arma::vec y_vec = arma::vectorise(y);
  
  double rss = arma::dot(err_vec, err_vec);
  double y_mean = arma::mean(y_vec);
  arma::vec y_demeaned = y_vec - y_mean;
  double tss = arma::dot(y_demeaned, y_demeaned);
  double r2 = 1.0 - rss / tss;
  
  // Return results as struct
  OLSResult result;
  result.beta = beta;
  result.fitted = fitted;
  result.err = err;
  result.r2 = r2;
  result.fitted_partial = fitted_partial;
  return result;
}

//' Ordinary Least Squares Regression
//'
//' @param y Dependent variable matrix! not vector
//' @param X Independent variables matrix (T x N)
//' @param c Integer indicator for intercept (1 if intercept included, 0 otherwise)
//'
//' @return A list containing:
//'   \itemize{
//'     \item beta: Coefficient estimates ((N+1) x 1 if c=1, N x 1 if c=0)
//'     \item fitted: Fitted values (T x 1)
//'     \item err: Residuals (T x 1)
//'     \item r2: R-squared statistic (scalar)
//'     \item fitted_partial: Fitted values excluding intercept (T x 1)
//'   }
//'
//' @details
//' This function performs ordinary least squares (OLS) regression. The coefficient
//' estimates are computed using the normal equations:
//' \deqn{\hat{\beta} = (X'X)^{-1}X'y}
//'
//' The R-squared statistic measures the proportion of variance explained:
//' \deqn{R^2 = 1 - \frac{RSS}{TSS}}
//' where RSS is the residual sum of squares and TSS is the total sum of squares.
//'
//' If an intercept is included (c=1), fitted_partial contains the fitted values
//' excluding the intercept contribution, useful for assessing the explanatory
//' power of the regressors alone.
//'
//' @examples
//' \dontrun{
//' # Generate sample data
//' set.seed(123)
//' y <- rnorm(100)
//' X <- matrix(rnorm(200), 100, 2)
//' 
//' # OLS with intercept
//' result <- fOLS(y, X, c = 1)
//' print(result$beta)
//' print(result$r2)
//' 
//' # OLS without intercept
//' result_no_int <- fOLS(y, X, c = 0)
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fOLS(arma::mat y, arma::mat X, int c = 1) {
  // Call the C++ function
  OLSResult result = fOLS_cpp(y, X, c);
  
  // Return results as a list for R
  return Rcpp::List::create(
    Rcpp::Named("beta") = result.beta,
    Rcpp::Named("fitted") = result.fitted,
    Rcpp::Named("err") = result.err,
    Rcpp::Named("r2") = result.r2,
    Rcpp::Named("fitted_partial") = result.fitted_partial
  );
}
