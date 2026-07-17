#ifndef FLPDID_H
#define FLPDID_H

#include <RcppArmadillo.h>

// -----------------------------------------------------------------------
// LP-DiD (Local Projections Difference-in-Differences) — C++ engine
//
// Implements Dube, Girardi, Jordà & Taylor (DGJT), "A Local Projections
// Approach to Difference-in-Differences" (NBER WP 31184).
//
// Post horizon h:   y_{i,t+h} - b_{it} = beta_h * dD_{it} + X'_it g + tau_t + e
// Pre  horizon j:   y_{i,t-j} - b_{it} = beta_{-j} * dD_{it} + ...
//   with baseline b_{it} = y_{i,t-1} (default) or the pre-mean (PMD).
//
// Clean-control samples:
//   * absorbing:     treated dD==1; controls untreated at t+h (post) / t (pre)
//   * non-absorbing: recursive clean-control sets over a stabilization
//     window L (Stata missing-value semantics), with CCC flavors:
//       ccc = 0 : no CCS restriction (ANRR-style)
//       ccc = 1 : clean past (CCS_0) for treated and controls
//       ccc = 2 : controls additionally clean up to t+h (CCS_h)
//
// Weighting: OLS (variance-weighted ATT) or reweight = equally-weighted
// ATT via inverse implicit weights w = 1/(1 - p_t) (exogenous case).
//
// Inference: CR1 cluster-robust sandwich with Stata reghdfe conventions,
//   c = G/(G-1) * (N-1)/(N-K), K counting absorbed time FEs (fixef.K full).
//
// The horizon loop is parallelized with OpenMP over horizons; all caches
// (panel index, dD, baselines, CCS tables) are precomputed and read-only.
//
// Author: Dr. Muhsin Ciftci
// -----------------------------------------------------------------------

Rcpp::List fLPDID_cpp(
    const arma::vec&  y,
    const arma::vec&  treat,
    const arma::mat&  X,
    const arma::ivec& i_index,
    const arma::ivec& t_index,
    const arma::ivec& cl_index,
    int  pre_window,
    int  post_window,
    bool nonabsorbing,
    int  Lwin,
    int  ccc,
    bool pmd,
    bool reweight,
    int  n_threads
);

#endif
