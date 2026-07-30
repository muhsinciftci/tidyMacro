#ifndef FLPPANEL_H
#define FLPPANEL_H

#include <RcppArmadillo.h>
#include <vector>

// -----------------------------------------------------------------------
// Panel Local Projections — C++ engine
//
// Port of `panel_LP.m` (Almuzara & Sancibrián, "Micro Responses to Macro
// Shocks", NY Fed Staff Report 1090, 2024).
//
// Model at horizon h (per s-component j = 1..n_s):
//     y_{i,t+h} = alpha_{FE} + sum_j (s_{it,j} * X_it) * beta_h^j
//                 + W'_{it} * G^h
//                 + sum_{k=1..p} lag_k(y, s⊗X)' * D^h_k + u_{i,t+h}
//   (cumulative: replaces y_{i,t+h} by sum_{r=0..h} y_{i,t+r})
//
// Inference:
//   * small_sample = false → asymptotic time-clustered sandwich SEs
//         V = (X'X)^{-1} (sum_t Z_t Z_t') (X'X)^{-1}, Z_t = sum_i X_it u_it
//     with df = Inf (normal critical values).
//   * small_sample = true  → Imbens-Kolesar (2016, REStat) refinement,
//     per s-component, using the within-period projection X_t and the
//     hat-matrix P0 = I_T - X_t * (X_t'X_t)^{-1} * X_t'.
//
// The horizon loop is parallelized with OpenMP; y_h and the lag columns
// used across horizons are precomputed once, so each horizon runs
// independently on a shared read-only cache.
//
// Author: Dr. Muhsin Ciftci
// -----------------------------------------------------------------------

struct LPPanelResult {
  arma::mat  estimate;   // (H+1) x n_s
  arma::mat  SE;         // (H+1) x n_s
  arma::mat  df;         // (H+1) x n_s
  arma::ivec status;     // (H+1)   0 = OK, 1 = no complete rows,
                         //         2 = HDFE non-convergence
  arma::ivec converged;  // (H+1)   1 = HDFE converged, 0 = did not
};

LPPanelResult fLPPanel_internal(
    const arma::mat& y,
    const arma::mat& s,
    const arma::mat& X,
    const arma::mat& W,
    const arma::imat& FE,
    const arma::ivec& i_index,
    const arma::ivec& t_index,
    int  H,
    int  p_max,
    bool small_sample,
    bool cumulative,
    int  n_threads,
    bool verbose
);

Rcpp::List fLPPanel_cpp(
    const arma::mat& y,
    const arma::mat& s,
    const arma::mat& X,
    const arma::mat& W,
    const arma::imat& FE,
    const arma::ivec& i_index,
    const arma::ivec& t_index,
    int  H,
    int  p_max,
    bool small_sample,
    bool cumulative,
    int  n_threads,
    bool verbose
);

#endif // FLPPANEL_H
