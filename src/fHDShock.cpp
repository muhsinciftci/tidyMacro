// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_DONT_USE_OPENMP

#include "fHDShock.h"
#include "fVARPosterior.h"
#include <RcppArmadillo.h>

HDShockResult fHDCompute_cpp(const arma::mat& y,
                             const arma::mat& beta,
                             const arma::mat& B,
                             int p, int c,
                             const arma::mat& exog) {

    const arma::uword N    = y.n_cols;
    const arma::uword Np   = N * static_cast<arma::uword>(p);
    const arma::uword nex  = exog.is_empty() ? 0u : exog.n_cols;

    if (B.n_rows != N || B.n_cols != N) Rcpp::stop("B must be N x N.");
    if (beta.n_cols != N)               Rcpp::stop("beta must have one column per variable.");

    arma::mat Y, X;
    fVARDesign_cpp(y, p, c, exog, Y, X);
    const arma::uword nobs = Y.n_rows;
    if (beta.n_rows != X.n_cols) Rcpp::stop("beta and the VAR design are incompatible.");

    // Companion matrix of the endogenous lag block.
    arma::mat Fcomp(Np, Np, arma::fill::zeros);
    Fcomp.rows(0, N - 1) = beta.rows(c, c + Np - 1).t();
    if (Np > N) Fcomp.submat(N, 0, Np - 1, Np - N - 1) = arma::eye<arma::mat>(Np - N, Np - N);

    // Structural shocks e_t = B^{-1} u_t, stored as N x nobs.
    const arma::mat resid = Y - X * beta;
    arma::mat eps;
    if (!arma::solve(eps, B, resid.t())) Rcpp::stop("The impact matrix B is singular.");

    arma::mat B_big(Np, N, arma::fill::zeros);
    B_big.rows(0, N - 1) = B;

    // ---- contribution of each structural shock ------------------------
    // Companion recursion  s_t = B_big e_t + Fcomp s_{t-1}, run once per shock.
    arma::cube shock(nobs + p, N, N, arma::fill::zeros);
    arma::vec  state(Np);
    for (arma::uword j = 0; j < N; ++j) {
        state.zeros();
        for (arma::uword t = 0; t < nobs; ++t) {
            state = Fcomp * state + B_big.col(j) * eps(j, t);
            for (arma::uword i = 0; i < N; ++i) shock(t + p, i, j) = state(i);
        }
    }

    // ---- initial condition --------------------------------------------
    // X.row(0) holds y_{0}, ..., y_{-p+1}: the state vector at the last
    // pre-sample date, propagated forward with no forcing.
    arma::mat init(nobs + p, N, arma::fill::zeros);
    state = X.row(0).cols(c, c + Np - 1).t();
    for (arma::uword i = 0; i < N; ++i) init(p - 1, i) = state(i);
    for (arma::uword t = 0; t < nobs; ++t) {
        state = Fcomp * state;
        for (arma::uword i = 0; i < N; ++i) init(t + p, i) = state(i);
    }

    // ---- intercept -----------------------------------------------------
    arma::mat cons(nobs + p, N, arma::fill::zeros);
    if (c > 0) {
        arma::vec CC(Np, arma::fill::zeros);
        CC.head(N) = beta.row(0).t();
        state.zeros();
        for (arma::uword t = 0; t < nobs; ++t) {
            state = Fcomp * state + CC;
            for (arma::uword i = 0; i < N; ++i) cons(t + p, i) = state(i);
        }
    }

    // ---- exogenous regressors ------------------------------------------
    arma::cube exo(nobs + p, N, std::max<arma::uword>(nex, 1u), arma::fill::zeros);
    if (nex > 0) {
        const arma::mat Xex = exog.rows(p, exog.n_rows - 1);
        for (arma::uword e = 0; e < nex; ++e) {
            arma::vec EE(Np, arma::fill::zeros);
            EE.head(N) = beta.row(c + Np + e).t();
            state.zeros();
            for (arma::uword t = 0; t < nobs; ++t) {
                state = Fcomp * state + EE * Xex(t, e);
                for (arma::uword i = 0; i < N; ++i) exo(t + p, i, e) = state(i);
            }
        }
    }

    // ---- total ----------------------------------------------------------
    arma::mat endo = init + cons;
    for (arma::uword j = 0; j < N; ++j) endo += shock.slice(j);
    if (nex > 0) for (arma::uword e = 0; e < nex; ++e) endo += exo.slice(e);

    // The first p rows are pre-sample; NaN-pad them so plots align with y.
    // `init` carries the state at row p-1, which is observable.
    const double nan = arma::datum::nan;
    for (arma::uword t = 0; t + 1 < static_cast<arma::uword>(p); ++t) {
        for (arma::uword i = 0; i < N; ++i) init(t, i) = nan;
    }
    for (arma::uword t = 0; t < static_cast<arma::uword>(p); ++t) {
        for (arma::uword i = 0; i < N; ++i) {
            for (arma::uword j = 0; j < N; ++j) shock(t, i, j) = nan;
            cons(t, i) = nan;
            endo(t, i) = nan;
            if (nex > 0) for (arma::uword e = 0; e < nex; ++e) exo(t, i, e) = nan;
        }
    }

    HDShockResult out;
    out.shock = shock;
    out.init  = init;
    out.cons  = cons;
    out.exo   = exo;
    out.endo  = endo;
    return out;
}

//' Historical Decomposition for a Sign-Restricted SVAR
//'
//' Decomposes each observed series into the cumulated contribution of every
//' structural shock plus the initial condition, the intercept and any exogenous
//' regressors, given a structural impact matrix.
//'
//' @param y A T x N numeric matrix of endogenous variables.
//' @param beta A (Np + c + M) x N coefficient matrix, as returned by
//'   \code{fVAR} or drawn by \code{fVARPosterior_cpp}.
//' @param B An N x N structural impact matrix, for instance \code{Bfp} from
//'   \code{fSR_cpp}. It must be invertible.
//' @param p Integer lag order.
//' @param c Integer intercept indicator (1 = include, 0 = exclude).
//' @param exog Optional T x M matrix of exogenous regressors (default NULL).
//'
//' @return A list with \code{shock}, a T x N x N array indexed
//'   \code{[time, variable, shock]}; \code{init}, \code{const} and \code{endo},
//'   each T x N; and \code{exo}, a T x N x M array. The first \code{p} rows are
//'   \code{NA} because they are absorbed as initial conditions. The components
//'   sum to \code{y} over the estimation sample.
//'
//' @details
//' Contributions are accumulated through the companion-form recursion
//' \eqn{s_t = F s_{t-1} + \tilde{B} e_t} with \eqn{e_t = B^{-1} u_t}, run once
//' per structural shock. This is the decomposition of \code{compute_HD.m} in
//' the VAR Toolbox; the trend block is omitted because \code{fVAR} supports an
//' intercept only.
//'
//' @seealso \code{\link{fSR_cpp}}, \code{\link{fHistDec}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fHDShock_cpp(const arma::mat& y,
                        const arma::mat& beta,
                        const arma::mat& B,
                        int p, int c,
                        Rcpp::Nullable<arma::mat> exog = R_NilValue) {

    arma::mat exog_mat;
    if (exog.isNotNull()) exog_mat = Rcpp::as<arma::mat>(exog);

    HDShockResult r = fHDCompute_cpp(y, beta, B, p, c, exog_mat);

    return Rcpp::List::create(Rcpp::Named("shock") = r.shock,
                              Rcpp::Named("init")  = r.init,
                              Rcpp::Named("const") = r.cons,
                              Rcpp::Named("exo")   = r.exo,
                              Rcpp::Named("endo")  = r.endo);
}
