#ifndef FLP_H
#define FLP_H

#include <RcppArmadillo.h>
#include <vector>

// -----------------------------------------------------------------------
// Classical Local Projections — C++ engine
//
//   Standard:   y_{t+h}           = a_h + B^h X_t + u_{t+h}
//   Cumulative: y_{t+h} - y_{t-1} = a_h + B^h X_t + u_{t+h}
//
// Cumulative base row (Y_pre). The long difference needs y_{t-1}, which for
// the first estimable date lies one period BEFORE row 0 of Y. Two input
// conventions are supported:
//
//   Y_pre is 1 x n_y : the outcome at the date immediately preceding Y row 0.
//                      Every row of Y is then an estimable LHS date and no
//                      observation is lost. This reproduces the endo_lag1
//                      alignment of Cesa-Bianchi's LPmodel.m
//                      (endo_lag1 = ENDO(lag_trim:end-1)).
//   Y_pre is empty   : Y row 0 is consumed as the base and estimation starts
//                      at row 1, costing one observation. Correct only when
//                      no earlier outcome exists in the data.
//
// IRF at horizon h = coefficient on the shock variable in B^h.
//
// Newey-West HAC bandwidth: nwLag = min(max(nw_lags_base + h + nw_offset, 0), T_h - 1).
//
//   nw_offset = +1 (default): Miranda-Agrippino & Ricco rule, bandwidth
//               grows as h+1.
//   nw_offset =  0: classic Jorda (2005) rule, bandwidth = h (LP residuals
//               at horizon h are MA(h) under H0 of correct specification).
//               fLP's h is 0-indexed (h=0 is the impact horizon), which
//               makes this the direct equivalent of Cesa-Bianchi's VAR
//               Toolbox (LPmodel.m: OLSmodel(Y,X,0,hh-1) with 1-indexed
//               hh) used e.g. in the Jorda-Taylor (2025) replication.
//
// Variance: full Newey-West HAC sandwich (only supported form)
//   V_h = (X'X)^{-1} * ( sum_{a=0}^{nwL} w_a * (Gamma_a + Gamma_a', a>0) ) * (X'X)^{-1}
//   with Gamma_a = sum_t (X_t u_t)(X_{t-a} u_{t-a})' and Bartlett weight
//   w_a = (nwL + 1 - a) / (nwL + 1). Matches Cesa-Bianchi's OLSmodel.m.
//
// Horizon loop is parallelized via OpenMP (each horizon h is independent).
//
// Author: Dr. Muhsin Ciftci
// -----------------------------------------------------------------------

struct LPResult {
  arma::mat irfs;        // (H+1) x n_y — IRF point estimates
  arma::mat irfs_upper;  // (H+1) x n_y — upper confidence band (at conf_level)
  arma::mat irfs_lower;  // (H+1) x n_y — lower confidence band (at conf_level)
  arma::mat irfs_se;     // (H+1) x n_y — raw (unscaled) SE of the shock
                         //   coefficient. Independent of conf_level, so R
                         //   can build additional confidence bands at any
                         //   level post-hoc without re-running the engine:
                         //   band = irfs +/- qnorm(0.5*(1+level)) * irfs_se.
  // Rank diagnostic. inv_sympd() on X'X fails exactly when [1, X] is
  // numerically rank deficient; the engine still returns (pseudo-inverse)
  // numbers so the parallel loop can finish, but the R wrapper stops rather
  // than reporting a min-norm coefficient and a HAC SE that has no
  // interpretation. rank_fail_h is the first offending horizon (-1 if none).
  bool rank_deficient = false;
  int  rank_fail_h    = -1;

  // Only populated when store_full = true:
  std::vector<arma::mat> betas;   // H+1 matrices (kr x n_y) — full coefficients
  std::vector<arma::mat> ses;     // H+1 matrices (kr x n_y) — full SEs
};

// Internal C++ function (callable from other translation units).
// Defaults live in the .cpp definitions (project convention).
LPResult fLP_internal(
    const arma::mat& Y,
    const arma::mat& X,
    int              H,
    int              shock_col,
    double           conf_level,
    int              nw_lags_base,
    bool             store_full,
    bool             cumulative,
    int              n_threads,
    int              nw_offset,
    bool             verbose,
    const arma::mat& Y_pre
);

// R-callable wrapper
Rcpp::List fLP_cpp(
    const arma::mat& Y,
    const arma::mat& X,
    int       H,
    int       shock_col,
    double    conf_level,
    int       nw_lags_base,
    bool      store_full,
    bool      cumulative,
    int       n_threads,
    int       nw_offset,
    bool      verbose,
    Rcpp::Nullable<arma::mat> Y_pre
);

#endif // FLP_H
