#include <RcppArmadillo.h>
#include "fOLS.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Internal C++ function (called from other C++ code)
OLSResult fOLS_cpp(arma::mat y, arma::mat X, int c, int robust, int lag, int mode) {

  // Get dimensions
  int T = X.n_rows;
  int N = X.n_cols;

  // Set default lag for HAC (Newey-West rule of thumb)
  if (robust == 2 && lag == 0) {
    lag = (int)std::round(4.0 * std::pow((double)T / 100.0, 2.0 / 9.0));
  }

  // Add intercept if c == 1 (prepended as first column)
  arma::mat X_reg;
  if (c == 1) {
    X_reg = arma::join_rows(arma::ones<arma::vec>(T), X);
  } else {
    X_reg = X;
  }

  int k_reg = X_reg.n_cols;  // Total regressors (including constant)

  // Compute invXpX and OLS estimates
  arma::mat invXpX = arma::inv_sympd(X_reg.t() * X_reg);
  arma::mat beta   = invXpX * (X_reg.t() * y);

  // Lightweight mode: skip r2, varbhat, F-stat — used in hot inner loops
  if (mode == 1) {
    OLSResult result;
    result.beta = beta;
    result.err  = y - X_reg * beta;
    return result;
  }

  // Fitted values
  arma::mat fitted = X_reg * beta;

  // Fitted values excluding intercept
  arma::mat fitted_partial;
  if (c == 1) {
    fitted_partial = X_reg.cols(1, N) * beta.rows(1, N);
  } else {
    fitted_partial = fitted;
  }

  // Residuals and sums of squares
  arma::mat err      = y - fitted;
  arma::vec err_vec  = arma::vectorise(err);
  arma::vec y_vec    = arma::vectorise(y);
  double SSR         = arma::dot(err_vec, err_vec);
  double y_mean      = arma::mean(y_vec);
  arma::vec y_dem    = y_vec - y_mean;
  double tss         = arma::dot(y_dem, y_dem);

  // R2, adjusted R2, RMSE
  double r2    = 1.0 - SSR / tss;
  double r2adj = 1.0 - ((double)(T - 1) / (double)(T - k_reg)) * (SSR / tss);
  double rmse  = std::sqrt(SSR / (double)(T - k_reg));

  // Standard variance-covariance
  arma::mat varbhat = rmse * rmse * invXpX;

  // Robust variance-covariance (White or HAC/Newey-West)
  arma::mat Shat        = arma::zeros(k_reg, k_reg);
  arma::mat varbhatrobust;
  if (robust >= 1) {
    for (int i = 0; i < T; i++) {
      arma::vec xi = X_reg.row(i).t();
      Shat += err_vec(i) * err_vec(i) * xi * xi.t();
    }
    if (robust == 2) {
      for (int l = 1; l <= lag; l++) {
        double w = 1.0 - (double)l / ((double)lag + 1.0);  // Bartlett weight
        arma::mat Shatl = arma::zeros(k_reg, k_reg);
        for (int i = l; i < T; i++) {
          arma::vec xi  = X_reg.row(i).t();
          arma::vec xil = X_reg.row(i - l).t();
          Shatl += w * err_vec(i) * err_vec(i - l) * (xi * xil.t() + xil * xi.t());
        }
        Shat += Shatl;
      }
    }
    varbhatrobust = ((double)T / (double)(T - k_reg)) * invXpX * Shat * invXpX;
  }

  // F-test for overall significance
  // Constant is column 0 (if c==1); test all slope coefficients == 0
  int num_restr = k_reg - c;  // N slope coefficients (or k_reg if no constant)
  double F_val = 0.0, Frobust_val = 0.0;

  if (num_restr > 0 && y.n_cols == 1) {
    // Selection matrix R: picks non-constant columns of beta
    arma::mat R = arma::zeros(num_restr, k_reg);
    if (c == 1) {
      R.cols(1, k_reg - 1) = arma::eye(num_restr, num_restr);
    } else {
      R = arma::eye(k_reg, k_reg);
    }

    arma::mat Rbeta = R * beta;
    arma::mat RVR   = R * varbhat * R.t();
    F_val = arma::as_scalar((1.0 / num_restr) * Rbeta.t() * arma::inv(RVR) * Rbeta);

    if (robust >= 1) {
      arma::mat RVRr = R * varbhatrobust * R.t();
      Frobust_val    = arma::as_scalar((1.0 / num_restr) * Rbeta.t() * arma::inv(RVRr) * Rbeta);
    }
  }

  // Return results as struct
  OLSResult result;
  result.beta          = beta;
  result.fitted        = fitted;
  result.err           = err;
  result.r2             = r2;
  result.r2adj          = r2adj;
  result.F              = F_val;
  result.Frobust        = Frobust_val;
  result.fitted_partial = fitted_partial;
  return result;
}

//' Ordinary Least Squares Regression
//'
//' @param y Dependent variable matrix (T x 1)
//' @param X Independent variables matrix (T x N)
//' @param c Integer indicator for intercept (1 if included, 0 otherwise)
//' @param robust SE type: 0 = standard, 1 = White (heteroskedasticity-robust),
//'   2 = Newey-West HAC
//' @param lag Number of lags for HAC (0 = Newey-West rule of thumb)
//'
//' @return A list containing:
//'   \itemize{
//'     \item beta: Coefficient estimates
//'     \item fitted: Fitted values
//'     \item err: Residuals
//'     \item r2: R-squared
//'     \item r2adj: Adjusted R-squared
//'     \item F: F-statistic for overall significance
//'     \item Frobust: Robust F-statistic (only if robust > 0)
//'     \item fitted_partial: Fitted values excluding intercept
//'   }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fOLS(arma::mat y, arma::mat X, int c = 1, int lag = 0) {
  // lag=0: White robust F; lag>0: Newey-West HAC with that lag
  int robust = (lag > 0) ? 2 : 1;
  OLSResult result = fOLS_cpp(y, X, c, robust, lag, 0);

  return Rcpp::List::create(
    Rcpp::Named("beta")           = result.beta,
    Rcpp::Named("fitted")         = result.fitted,
    Rcpp::Named("err")            = result.err,
    Rcpp::Named("r2")             = result.r2,
    Rcpp::Named("r2adj")          = result.r2adj,
    Rcpp::Named("F")              = result.F,
    Rcpp::Named("Frobust")        = result.Frobust,
    Rcpp::Named("fitted_partial") = result.fitted_partial
  );
}
