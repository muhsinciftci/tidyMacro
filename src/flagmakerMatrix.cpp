#include "flagmakerMatrix.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

//' Create Lagged Matrix
//'
//' Creates a matrix of lagged values from a time series matrix. For each
//' variable and each lag from 1 to p, this function creates a column containing
//' the lagged values, removing the first p rows to align the lags properly.
//'
//' @param y A numeric matrix, data frame, or tibble of time series data with T
//'   rows (time periods) and N columns (variables)
//' @param p An integer specifying the number of lags to create
//'
//' @return A numeric matrix with (T-p) rows and (N*p) columns. The columns are
//'   organized as: all N variables at lag 1, then all N variables at lag 2, etc.
//'   Column names are generated as "var_name_lag_k" where applicable.
//'
//' @details The function creates a lagged matrix where each column represents a
//'   lagged version of the original variables. The lags are organized by lag
//'   order first, then by variable. For example, with 2 variables and 2 lags,
//'   the column order would be: var1_lag1, var2_lag1, var1_lag2, var2_lag2.
//'
//' @examples
//' \dontrun{
//' # Create sample data as data frame
//' y <- data.frame(var1 = 1:10, var2 = 11:20)
//'
//' # Create 2 lags
//' x <- lagmakerMatrix(y, 2)
//'
//' # Works with tibbles too
//' library(tibble)
//' y_tbl <- tibble(var1 = 1:10, var2 = 11:20)
//' x <- lagmakerMatrix(y_tbl, 2)
//' }
//'
//' @export
// [[Rcpp::export]]
arma::mat flagmakerMatrix(arma::mat y, int p) {

  int T = y.n_rows;
  int N = y.n_cols;

  // Create output matrix with (T-p) rows and (N*p) columns
  arma::mat x(T - p, N * p);

  int counter = 0;

  // Loop through lags
  for (int i = 0; i < p; i++) {
    // Loop through variables
    for (int j = 0; j < N; j++) {
      // Extract lagged values: y(p-i:T-1-i, j) in 0-indexed
      // In MATLAB: y(p+1-i:T-i,j) is 1-indexed
      // In C++: y.submat(p-i-1, j, T-i-2, j) for 0-indexed
      x.col(counter) = y.submat(p - i - 1, j, T - i - 2, j);
      counter++;
    }
  }

  return x;
}
