#include <RcppArmadillo.h>
#include "fbootstrapVAR.h"
#include "fVAR.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Internal C++ function (called from other C++ code)
// Accepts VARResult struct directly - no list overhead!
BootstrapVARResult fbootstrapVAR_cpp(const arma::mat&    y,
                                     const VARResult&    var_result,
                                     const std::string&  bootscheme) {

    // ------------------------------------------------------------------ //
    //  Unpack inputs from struct (much faster than from Rcpp::List)       //
    // ------------------------------------------------------------------ //
    const arma::mat& beta      = var_result.beta;
    const arma::mat& residuals = var_result.residuals;
    const int        p         = var_result.p;
    const int        c         = var_result.c;

    const int T = y.n_rows;  // Implicit conversion from arma::uword to int
    const int N = y.n_cols;  // Implicit conversion from arma::uword to int
    const int T_resid = residuals.n_rows;  // Implicit conversion from arma::uword to int
    const int T_iter = T - p;

    // ------------------------------------------------------------------ //
    //  Separate intercept and lag-coefficient blocks                       //
    // ------------------------------------------------------------------ //
    // beta layout (rows):  [intercept (if c=1)] | [A_1 ; A_2 ; ... ; A_p]
    // Each A_i is N x N, stored as N rows of beta.
    arma::rowvec intercept;
    arma::mat    Pi;             // (N*p) x N

    if (c == 1) {
        intercept = beta.row(0);                      // 1 x N
        Pi        = beta.rows(1, beta.n_rows - 1);    // (N*p) x N
    } else {
        intercept = arma::zeros<arma::rowvec>(N);
        Pi        = beta;                              // (N*p) x N
    }

    // ------------------------------------------------------------------ //
    //  Build bootstrap residual matrix                                     //
    // ------------------------------------------------------------------ //
    arma::mat boot_resid(T_iter, N, arma::fill::none);
    arma::vec rademacher(T_iter, arma::fill::zeros);

    if (bootscheme == "residual") {
        // iid resample rows of the residual matrix
        arma::ivec idx = arma::randi<arma::ivec>(
                             T_iter, arma::distr_param(0, T_resid - 1));
        for (int i = 0; i < T_iter; ++i) {
            boot_resid.row(i) = residuals.row(idx(i));
        }

    } else if (bootscheme == "wild") {
        // Rademacher signs: +1 or -1 with equal probability
        arma::vec u = arma::randu<arma::vec>(T_resid);
        rademacher  = arma::sign(u - 0.5);
        // Multiply each residual row by its sign
        for (int i = 0; i < T_iter; ++i) {
            const int ri  = i % T_resid;
            boot_resid.row(i) = residuals.row(ri) * rademacher(ri);
        }
    }
    // (any other scheme: boot_resid is zero — caller's responsibility)

    // ------------------------------------------------------------------ //
    //  Recursive simulation                                                //
    // ------------------------------------------------------------------ //
    arma::mat ynext(T, N, arma::fill::none);
    ynext.rows(0, p - 1) = y.rows(0, p - 1);   // seed with original data

    // lag vector: [y_{t-1}' , y_{t-2}' , ... , y_{t-p}']  (column, length N*p)
    arma::vec ylag(N * p, arma::fill::none);
    for (int lag = 0; lag < p; ++lag) {
        ylag.rows(lag * N, lag * N + N - 1) = y.row(p - 1 - lag).t();
    }

    for (int i = 0; i < T_iter; ++i) {
        // y_t = intercept + Pi' * ylag + e_t*
        // Pi is (N*p) x N, so Pi.t() * ylag gives N x 1
        arma::rowvec yt = intercept + (Pi.t() * ylag).t() + boot_resid.row(i);
        ynext.row(p + i) = yt;

        // Update lag vector: shift down and insert new observation at top
        if (p > 1) {
            ylag.rows(N, N * p - 1) = ylag.rows(0, N * (p - 1) - 1);
        }
        ylag.rows(0, N - 1) = yt.t();
    }

    // Return results as struct
    BootstrapVARResult result;
    result.ynext = ynext;
    result.rademacher = rademacher;
    return result;
}

//' Bootstrap VAR Model
//'
//' Generates a bootstrap pseudo-sample from a fitted VAR model using either
//' residual (iid) resampling or a wild (Rademacher) bootstrap scheme.
//'
//' @param y T x N matrix of original endogenous variables.
//' @param fVAR_result List from \code{fVAR} containing beta, residuals, p,
//'   and c.
//' @param bootscheme String indicating the bootstrap method: "residual" for
//'   iid resampling of residuals, or "wild" for a Rademacher wild bootstrap.
//'
//' @return A list with two elements: ynext — T x N bootstrapped data matrix
//'   (first p rows copied from \code{y}); rademacher — T_eff-length vector of
//'   Rademacher signs used in the wild scheme (all zeros for residual scheme).
//'
//' @details
//' The recursion is \deqn{y_t = c + A_1 y_{t-1} + \cdots + A_p y_{t-p} + e_t^*}
//' where \eqn{e_t^*} is either a randomly resampled residual (residual scheme)
//' or the t-th residual multiplied by a Rademacher sign (wild scheme).
//'
//' @seealso \code{\link{fbootstrapChol}}, \code{\link{fVAR}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fbootstrapVAR(const arma::mat&    y,
                         const Rcpp::List&   fVAR_result,
                         const std::string&  bootscheme = "residual") {

    // Extract elements from R list and construct VARResult struct
    VARResult var_result;
    var_result.beta = Rcpp::as<arma::mat>(fVAR_result["beta"]);
    var_result.residuals = Rcpp::as<arma::mat>(fVAR_result["residuals"]);
    var_result.sigma_full = Rcpp::as<arma::mat>(fVAR_result["sigma_full"]);
    var_result.p = Rcpp::as<int>(fVAR_result["p"]);
    var_result.c = Rcpp::as<int>(fVAR_result["c"]);
    var_result.n_exog = Rcpp::as<int>(fVAR_result["n_exog"]);
    
    // Call the C++ function with the struct
    BootstrapVARResult result = fbootstrapVAR_cpp(y, var_result, bootscheme);
    
    // Return results as a list for R
    return Rcpp::List::create(
        Rcpp::Named("ynext")      = result.ynext,
        Rcpp::Named("rademacher") = result.rademacher
    );
}
