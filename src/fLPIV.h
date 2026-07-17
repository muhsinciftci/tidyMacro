#ifndef FLPIV_H
#define FLPIV_H

#include <RcppArmadillo.h>
#include <vector>

// -----------------------------------------------------------------------
// Local Projections with an External Instrument (LP-IV) — C++ engine
//
//   y_{t+h} - y_{t-1} = a_h + b_h * D_t + Gamma_h' * C_t + u_{t+h}    (cumulative)
//   y_{t+h}           = a_h + b_h * D_t + Gamma_h' * C_t + u_{t+h}    (standard)
//
// where D_t is the endogenous treatment (e.g. FFR) and the external
// instrument matrix Z_t (already containing all requested lags) identifies
// b_h via 2SLS. Estimated horizon-by-horizon via Frisch-Waugh-Lovell:
//
//   1. FWL-residualise Y_h, D, Z on [1, C]                     -> Y_r, D_r, Z_r
//   2. First stage:  D_hat = Z_r * (Z_r' Z_r)^{-1} Z_r' D_r
//   3. 2SLS:         b_h  = (D_hat' Y_r) / (D_hat' D_r)
//                    eps  = Y_r - b_h * D_r
//   4. Newey-West HAC delta-method SE:
//          g_t = D_hat_t * eps_t             (scalar score per t)
//          LRV = Gamma_0 + sum_{j=1..bw} w_j (Gamma_j + Gamma_j)
//                on demeaned g_t (Bartlett weight w_j = 1 - j/(bw+1))
//          se  = sqrt(LRV) / |D_hat' D_r|
//
// Bandwidth rule (matches Cesa-Bianchi's LPmodel.m LP-IV path):
//   - nw_lags_iv > 0  : fixed bandwidth = nw_lags_iv for all horizons
//                       (Stata's vce(hac nw <k>))
//   - nw_lags_iv == 0 : horizon-varying bandwidth = h  (just-identified rule)
//
// Author: Dr. Muhsin Ciftci
// -----------------------------------------------------------------------

struct LPIVResult {
  arma::mat irfs;         // (H+1) x n_y — IRF point estimates
  arma::mat irfs_upper;   // (H+1) x n_y
  arma::mat irfs_lower;   // (H+1) x n_y
  arma::mat irfs_se;      // (H+1) x n_y — raw NW HAC SE (unscaled by z_crit)
  arma::vec Fstat_fs;     // (H+1)       — first-stage F stat (per horizon)
  arma::vec rsqr_fs;      // (H+1)       — first-stage R^2 (per horizon)
};

LPIVResult fLPIV_internal(
    const arma::mat& Y,          // T x n_y outcomes
    const arma::vec& D,          // T endogenous treatment
    const arma::mat& Z,          // T x n_z instrument set (with any lags)
    const arma::mat& C,          // T x n_c exogenous controls (no intercept)
    int              H,
    double           conf_level,
    int              nw_lags_iv,
    bool             cumulative,
    int              n_threads
);

Rcpp::List fLPIV_cpp(
    const arma::mat& Y,
    const arma::vec& D,
    const arma::mat& Z,
    const arma::mat& C,
    int       H,
    double    conf_level,
    int       nw_lags_iv,
    bool      cumulative,
    int       n_threads
);

#endif // FLPIV_H
