#include "fVAR.h"
#include "flagmakerMatrix.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// Shared implementation: takes either an empty mat (no exog) or a T x M exog mat
// Kept static so the public _cpp overloads stay as the only external entry points.
static VARResult fVAR_cpp_impl(const arma::mat &y, int p, int c,
                               const arma::mat *exog_ptr) {
  int T = y.n_rows;

  arma::mat yfinal = y.rows(p, T - 1);
  arma::mat lags = flagmakerMatrix(y, p);

  arma::mat X;
  if (c == 1) {
    X = arma::join_rows(arma::ones<arma::mat>(T - p, 1), lags);
  } else {
    X = lags;
  }

  int n_exog = 0;
  if (exog_ptr != nullptr) {
    const arma::mat &exog_mat = *exog_ptr;
    arma::mat exog_eff = exog_mat.rows(p, T - 1);
    X = arma::join_rows(X, exog_eff);
    n_exog = exog_mat.n_cols;
  }

  arma::mat XtX = X.t() * X;
  arma::mat Xty = X.t() * yfinal;
  arma::mat beta = arma::solve(XtX, Xty, arma::solve_opts::fast);

  arma::mat residuals = yfinal - X * beta;
  int n_obs = residuals.n_rows;   // T_full - p (residual sample)
  int N     = y.n_cols;
  arma::mat sigma = (residuals.t() * residuals) / (n_obs - c - N * p - n_exog);

  VARResult result;
  result.beta = beta;
  result.residuals = residuals;
  result.sigma = sigma;
  result.p = p;
  result.c = c;
  result.n_exog = n_exog;
  return result;
}

// Internal C++ function (called from other C++ code) — Nullable overload
VARResult fVAR_cpp(const arma::mat &y, int p, int c,
                   Rcpp::Nullable<arma::mat> exog) {
  if (exog.isNotNull()) {
    arma::mat exog_mat = Rcpp::as<arma::mat>(exog);
    return fVAR_cpp_impl(y, p, c, &exog_mat);
  }
  return fVAR_cpp_impl(y, p, c, nullptr);
}

// Pre-converted exog overload: avoids Rcpp::as in hot loops
VARResult fVAR_cpp(const arma::mat &y, int p, int c,
                   const arma::mat &exog) {
  return fVAR_cpp_impl(y, p, c, &exog);
}


//' Vector Autoregression (VAR) Model Estimation
//'
//' Estimates a Vector Autoregression model with optional exogenous variables
//' using ordinary least squares (OLS). Optimised for maximum speed and memory
//' efficiency via RcppArmadillo.
//'
//' @param y A numeric matrix of time series data with T rows (time periods)
//'   and N columns (endogenous variables).
//' @param p An integer specifying the number of lags to include in the VAR.
//' @param c An integer indicator for including a constant term (1 = include
//'   intercept, 0 = no intercept).
//' @param exog An optional T x M matrix of exogenous variables. If provided,
//'   these variables enter contemporaneously (not lagged). Default is NULL.
//'
//' @return A list with elements: beta — coefficient matrix of dimensions
//'   (Np + c + M) x N, where the first row is the intercept (if c = 1),
//'   followed by N*p lag-coefficient rows, then M exogenous-coefficient rows;
//'   residuals — (T-p) x N matrix of OLS residuals;
//'   sigma — N x N residual covariance matrix with degrees-of-freedom
//'   correction: denominator is (T - p) - c - N*p - n_exog, matching
//'   Bloom (2009) and Kaenzig (2021);
//'   p — lag order (echoed from input);
//'   c — intercept indicator (echoed from input);
//'   n_exog — number of exogenous variables (0 if none provided).
//'
//' @details
//' The model is \deqn{Y_t = c + A_1 Y_{t-1} + \cdots + A_p Y_{t-p} + B X_t + e_t}
//' where \eqn{X_t} are contemporaneous exogenous variables. Without exogenous
//' variables the model reduces to a standard VAR(p). Estimation uses
//' \code{arma::solve} for numerical stability.
//'
//' @examples
//' \dontrun{
//' # VAR(2) without exogenous variables
//' y      <- matrix(rnorm(200), ncol = 2)
//' result <- fVAR(y, p = 2, c = 1)
//'
//' # VAR(2) with one exogenous variable
//' exog   <- matrix(rnorm(100), ncol = 1)
//' result <- fVAR(y, p = 2, c = 1, exog = exog)
//' }
//'
//' @seealso \code{\link{fwoldIRF}}, \code{\link{fbootstrapChol}},
//'   \code{\link{fAICBIC}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fVAR(const arma::mat &y, int p, int c,
                Rcpp::Nullable<arma::mat> exog = R_NilValue) {
  VARResult result = fVAR_cpp(y, p, c, exog);

  return Rcpp::List::create(Rcpp::Named("beta") = result.beta,
                            Rcpp::Named("residuals") = result.residuals,
                            Rcpp::Named("sigma") = result.sigma,
                            Rcpp::Named("p") = result.p,
                            Rcpp::Named("c") = result.c,
                            Rcpp::Named("n_exog") = result.n_exog);
}
