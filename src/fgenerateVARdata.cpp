#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Generate VAR Data
//'
//' Recursively generates VAR data from initial values, coefficients, and
//' residuals. Optimised for maximum performance and memory efficiency.
//'
//' @param y T x N matrix of initial/historical endogenous variables.
//' @param p Integer lag order of the VAR model.
//' @param c Integer indicator for intercept (1 if included, 0 otherwise).
//' @param beta Coefficient matrix: (Np+c) x N if c = 1, or Np x N if c = 0.
//' @param residuals (T-p) x N matrix of residuals/shocks to add.
//'
//' @return T x N matrix of generated VAR data. The first p rows are copied
//'   from \code{y}; the remaining T-p rows are generated recursively.
//'
//' @details
//' The recursion follows:
//' \deqn{y_t = c + A_1 y_{t-1} + \cdots + A_p y_{t-p} + e_t}
//'
//' @examples
//' \dontrun{
//' y_init    <- matrix(rnorm(200), ncol = 2)
//' beta      <- matrix(rnorm(10), 5, 2)
//' residuals <- matrix(rnorm(190), 95, 2)
//' y_new     <- fgenerateVARdata(y_init, p = 2, c = 1, beta, residuals)
//' }
//'
//' @export
// [[Rcpp::export]]
arma::mat fgenerateVARdata(const arma::mat& y,
                           int p,
                           int c,
                           const arma::mat& beta,
                           const arma::mat& residuals) {

    const int T  = y.n_rows;
    const int N  = y.n_cols;
    const int Np = N * p;

    arma::mat ynext(T, N, arma::fill::none);
    ynext.rows(0, p - 1) = y.rows(0, p - 1);

    arma::rowvec const_term;
    arma::mat    pi;

    if (c == 1) {
        const_term = beta.row(0);
        pi         = beta.rows(1, beta.n_rows - 1);
    } else {
        const_term = arma::zeros<arma::rowvec>(N);
        pi         = beta;
    }

    arma::vec ylag_vec(Np);
    for (int lag = 0; lag < p; ++lag) {
        for (int j = 0; j < N; ++j) {
            ylag_vec(lag * N + j) = y(p - 1 - lag, j);
        }
    }

    const int T_iter = T - p;
    for (int i = 0; i < T_iter; ++i) {
        for (int j = 0; j < N; ++j) {
            double sum = const_term(j);
            for (int k = 0; k < Np; ++k) {
                sum += pi(k, j) * ylag_vec(k);
            }
            sum += residuals(i, j);
            ynext(p + i, j) = sum;
        }
        for (int lag = 0; lag < p; ++lag) {
            for (int j = 0; j < N; ++j) {
                ylag_vec(lag * N + j) = ynext(p + i - lag, j);
            }
        }
    }

    return ynext;
}
