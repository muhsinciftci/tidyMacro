// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_DONT_USE_OPENMP

#include "fVARPosterior.h"
#include "fVAR.h"
#include "fLagMakerMatrix.h"
#include <RcppArmadillo.h>

void fVARDesign_cpp(const arma::mat& y, int p, int c, const arma::mat& exog,
                    arma::mat& Y, arma::mat& X) {
    const arma::uword T = y.n_rows;
    Y = y.rows(p, T - 1);

    arma::mat lags = fLagMakerMatrix(y, p);
    if (c == 1) {
        X = arma::join_rows(arma::ones<arma::mat>(T - p, 1), lags);
    } else {
        X = lags;
    }
    if (!exog.is_empty()) {
        X = arma::join_rows(X, exog.rows(p, T - 1));
    }
}

NIWPosterior fNIWPosteriorPrep_cpp(const arma::mat& beta_hat,
                                   const arma::mat& sigma_hat,
                                   const arma::mat& XtX,
                                   int nobs) {

    NIWPosterior post;
    post.N    = static_cast<int>(sigma_hat.n_cols);
    post.k    = static_cast<int>(beta_hat.n_rows);
    post.nobs = nobs;

    if (nobs <= post.N) {
        Rcpp::stop("Posterior draws require nobs > nvar; the inverse-Wishart "
                   "degrees of freedom (%d) are too few for %d variables.",
                   nobs, post.N);
    }
    if (static_cast<int>(XtX.n_rows) != post.k) {
        Rcpp::stop("XtX and beta_hat have incompatible dimensions.");
    }

    post.beta_hat = beta_hat;

    // Psi = nobs * sigma_hat = Lsig * Lsig'
    arma::mat Lsig;
    if (!arma::chol(Lsig, static_cast<double>(nobs) * sigma_hat, "lower")) {
        Rcpp::stop("The residual covariance matrix is not positive definite.");
    }
    post.Lsig = Lsig;

    // (X'X)^{-1} = Rx^{-1} Rx^{-T} with X'X = Rx' Rx (upper Cholesky).
    arma::mat Rx;
    if (!arma::chol(Rx, XtX, "upper")) {
        Rcpp::stop("X'X is not positive definite; the VAR design is rank deficient.");
    }
    post.Lx = arma::inv(arma::trimatu(Rx));

    return post;
}

void fNIWPosteriorDraw_cpp(const NIWPosterior& post,
                           tidymacro::RNG&     rng,
                           NIWScratch&         scratch,
                           arma::mat&          G,
                           arma::mat&          sigma_draw,
                           arma::mat&          beta_draw) {

    const arma::uword N  = static_cast<arma::uword>(post.N);
    const arma::uword k  = static_cast<arma::uword>(post.k);
    const double      df = static_cast<double>(post.nobs);

    if (scratch.A.n_rows != N)  scratch.A.set_size(N, N);
    if (scratch.Z.n_rows != k || scratch.Z.n_cols != N) scratch.Z.set_size(k, N);

    // ---- 1. sigma ~ IW(nobs * sigma_hat, nobs) -------------------------
    // Bartlett factor A (lower triangular) of a Wishart(I, df) draw, then
    // sigma = (Lsig A^{-T}) (Lsig A^{-T})'.  Sampling the triangular factor
    // costs O(N^2) normals instead of the O(N * df) of a naive Wishart.
    scratch.A.zeros();
    for (arma::uword i = 0; i < N; ++i) {
        scratch.A(i, i) = std::sqrt(rng.chisq(df - static_cast<double>(i)));
        for (arma::uword j = 0; j < i; ++j) {
            scratch.A(i, j) = rng.norm();
        }
    }
    scratch.Ainv = arma::inv(arma::trimatl(scratch.A));

    G          = post.Lsig * scratch.Ainv.t();
    sigma_draw = G * G.t();
    sigma_draw = 0.5 * (sigma_draw + sigma_draw.t());   // enforce exact symmetry

    // ---- 2. beta | sigma ~ N(beta_hat, kron(sigma, (X'X)^{-1})) --------
    for (arma::uword j = 0; j < N; ++j) {
        double* zc = scratch.Z.colptr(j);
        for (arma::uword i = 0; i < k; ++i) zc[i] = rng.norm();
    }
    scratch.tmp = post.Lx * scratch.Z;      // k x N
    beta_draw   = post.beta_hat + scratch.tmp * G.t();
}

//' Draw from the Normal-Inverse-Wishart Posterior of a VAR
//'
//' Samples reduced-form VAR parameters from their exact posterior under a
//' flat (diffuse) prior, the sampler underlying Bayesian sign-restriction
//' inference.
//'
//' @param y A T x N numeric matrix of endogenous variables.
//' @param p Integer lag order.
//' @param c Integer intercept indicator (1 = include, 0 = exclude).
//' @param ndraws Integer number of posterior draws.
//' @param seed Integer base seed. Draw \code{d} uses \code{seed + d}, so the
//'   output does not depend on how the work is scheduled.
//' @param exog Optional T x M matrix of exogenous regressors (default NULL).
//'
//' @return A list with elements \code{beta_draws}, a (Np + c + M) x N x ndraws
//'   cube of coefficient draws, and \code{sigma_draws}, an N x N x ndraws cube
//'   of covariance draws.
//'
//' @details
//' The posterior is
//' \deqn{\Sigma \mid y \sim IW(T_{eff}\hat{\Sigma},\ T_{eff}), \qquad
//'       vec(B) \mid \Sigma, y \sim N(vec(\hat{B}),\ \Sigma \otimes (X'X)^{-1})}
//' with \eqn{T_{eff} = T - p}. The Kronecker covariance is never formed: the
//' coefficient draw is built from the Cholesky factors of \eqn{(X'X)^{-1}} and
//' of the drawn \eqn{\Sigma}.
//'
//' @references
//' Uhlig, H. (2005). What are the effects of monetary policy on output?
//' \emph{Journal of Monetary Economics}, 52(2), 381--419.
//'
//' @seealso \code{\link{fSR_cpp}}, \code{\link{fVAR}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fVARPosterior_cpp(const arma::mat& y, int p, int c, int ndraws,
                             int seed = 42,
                             Rcpp::Nullable<arma::mat> exog = R_NilValue) {

    if (ndraws <= 0) Rcpp::stop("ndraws must be positive.");

    arma::mat exog_mat;
    if (exog.isNotNull()) exog_mat = Rcpp::as<arma::mat>(exog);

    VARResult var_result = exog.isNotNull()
                             ? fVAR_cpp(y, p, c, exog_mat)
                             : fVAR_cpp(y, p, c, R_NilValue);

    arma::mat Y, X;
    fVARDesign_cpp(y, p, c, exog_mat, Y, X);

    const int nobs = static_cast<int>(X.n_rows);
    NIWPosterior post = fNIWPosteriorPrep_cpp(var_result.beta, var_result.sigma,
                                              X.t() * X, nobs);

    arma::cube beta_draws(post.k, post.N, ndraws, arma::fill::none);
    arma::cube sigma_draws(post.N, post.N, ndraws, arma::fill::none);

    NIWScratch scratch;
    arma::mat G, sigma_draw, beta_draw;
    for (int d = 0; d < ndraws; ++d) {
        tidymacro::RNG rng(static_cast<std::uint64_t>(seed) +
                           static_cast<std::uint64_t>(d));
        fNIWPosteriorDraw_cpp(post, rng, scratch, G, sigma_draw, beta_draw);
        beta_draws.slice(d)  = beta_draw;
        sigma_draws.slice(d) = sigma_draw;
    }

    return Rcpp::List::create(Rcpp::Named("beta_draws")  = beta_draws,
                              Rcpp::Named("sigma_draws") = sigma_draws);
}
