// -----------------------------------------------------------------------
// fLP.cpp — Classical Local Projections (C++ engine)
//
//   Standard (cumulative = false):
//     y_{t+h} = a_h + B^h * X_t + u_{t+h}
//
//   Cumulative (cumulative = true):
//     y_{t+h} - y_{t-1} = a_h + B^h * X_t + u_{t+h}
//
//   IRF at horizon h = row of B^h on the shock variable.
//
//   Newey-West HAC — full sandwich:
//     V_h = (X'X)^{-1} * G * (X'X)^{-1}
//     G   = w_0*Gamma_0 + sum_{a=1..nwL} w_a*(Gamma_a + Gamma_a')
//     Gamma_a = sum_{t=a..T_h-1} g_t g_{t-a}',    g_t = X_reg_t * u_t
//     w_a = (nwL + 1 - a) / (nwL + 1)             (Bartlett kernel)
//
//   Horizon loop is parallelized with OpenMP (each horizon is independent).
//
// Author: Dr. Muhsin Ciftci
// -----------------------------------------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// CRITICAL: Define this BEFORE including RcppArmadillo to prevent conflicts
#define ARMA_DONT_USE_OPENMP

#include "fLP.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <cstdio>

using namespace arma;


// =======================================================================
// Internal C++ implementation — called from R wrapper and (optionally)
// from other C++ translation units (e.g., bootstrap routines).
// =======================================================================

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
    bool             verbose
) {
  // ---------- dimensions -----------------------------------------------
  const int T  = static_cast<int>(Y.n_rows);
  const int ny = static_cast<int>(Y.n_cols);
  const int k  = static_cast<int>(X.n_cols);
  const int kr = k + 1;                          // +1 for constant (prepended)
  const int sc = shock_col + 1;                  // column index in X_reg (0 = const)

  const int offset = cumulative ? 1 : 0;         // cumulative needs y_{t-1}

  if (T <= H + offset || T - H - offset <= kr) {
    Rcpp::stop("fLP: not enough observations for the largest horizon.");
  }
  if (shock_col < 0 || shock_col >= k) {
    Rcpp::stop("fLP: shock_col must index a column of X.");
  }

  // Critical value for two-sided interval
  const double z_crit = R::qnorm(0.5 * (1.0 + conf_level), 0.0, 1.0, 1, 0);

  // ---------- output storage -------------------------------------------
  LPResult out;
  out.irfs       = arma::mat(H + 1, ny, arma::fill::zeros);
  out.irfs_upper = arma::mat(H + 1, ny, arma::fill::zeros);
  out.irfs_lower = arma::mat(H + 1, ny, arma::fill::zeros);
  out.irfs_se    = arma::mat(H + 1, ny, arma::fill::zeros);
  if (store_full) {
    out.betas.resize(H + 1);
    out.ses.resize(H + 1);
  }

  // ---------- threading setup ------------------------------------------
  // Unified rule across LP engines: n_threads <= 0 → use all available cores.
  int actual_threads = 1;
#ifdef _OPENMP
  actual_threads = (n_threads <= 0)
                       ? std::max(1, omp_get_max_threads())
                       : n_threads;
  if (verbose) {
    std::printf("fLP: using %d thread(s) for parallel horizon loop...\n",
                actual_threads);
  }
#else
  (void) n_threads;
  if (verbose) {
    std::printf("fLP: OpenMP not available. Running single-threaded.\n");
  }
#endif

  // ---------- horizon loop (parallel over h) ---------------------------
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(actual_threads)
#endif
  for (int h = 0; h <= H; h++) {

    // ---- align data: y_{t+h} (or y_{t+h} - y_{t-1}) vs X_t -----------
    arma::mat Yh;
    arma::mat Xh;
    if (cumulative) {
      // t = 1, ..., T - h - 1  →  T - h - 1 rows
      Yh = Y.rows(h + 1, T - 1) - Y.rows(0, T - h - 2);
      Xh = X.rows(1, T - h - 1);
    } else {
      // t = 0, ..., T - h - 1  →  T - h rows
      Yh = Y.rows(h, T - 1);
      Xh = X.rows(0, T - h - 1);
    }
    const int Th = static_cast<int>(Yh.n_rows);

    // ---- cross-products with implicit intercept -----------------------
    arma::mat XpX(kr, kr, arma::fill::zeros);
    XpX(0, 0) = static_cast<double>(Th);
    if (k > 0) {
      const arma::rowvec x_sum = arma::sum(Xh, 0);
      XpX.submat(0, 1, 0, kr - 1) = x_sum;
      XpX.submat(1, 0, kr - 1, 0) = x_sum.t();
      XpX.submat(1, 1, kr - 1, kr - 1) = Xh.t() * Xh;
    }

    arma::mat XpY(kr, ny, arma::fill::zeros);
    XpY.row(0) = arma::sum(Yh, 0);
    if (k > 0) {
      XpY.rows(1, kr - 1) = Xh.t() * Yh;
    }

    arma::mat invXpX;
    const bool inv_ok = arma::inv_sympd(invXpX, XpX);
    if (!inv_ok || !invXpX.is_finite()) {
      invXpX = arma::pinv(XpX);
    }

    // ---- OLS — all equations simultaneously --------------------------
    const arma::mat Beta = invXpX * XpY;             // kr x ny

    // ---- NW bandwidth: nwLags = nw_lags_base + h + nw_offset ---------
    //  nw_offset = +1 (default): Miranda-Agrippino & Ricco rule.
    //  nw_offset =  0: classic Jorda (2005) rule. fLP's h is 0-indexed
    //  (h=0 is the impact horizon), which matches Cesa-Bianchi's VAR
    //  Toolbox LPmodel.m OLSmodel(Y,X,0,hh-1) with 1-indexed hh.
    const int nwL = std::min(std::max(nw_lags_base + h + nw_offset, 0), Th - 1);

    // ---- equation-level storage for this horizon ---------------------
    arma::mat Se_h;
    if (store_full) {
      Se_h.set_size(kr, ny);
    }

    // ---- regressor matrix with explicit intercept column (shared
    //      across equations) --------------------------------------------
    const arma::mat Xreg = arma::join_rows(arma::ones<arma::vec>(Th), Xh);  // Th x kr

    // Row `sc` of invXpX; needed for both the full-sandwich and the
    // shock-only (fast) branches. Computing it once per horizon avoids
    // redoing the extraction inside the equation loop.
    const arma::rowvec invXpX_sc = invXpX.row(sc);

    // ---- loop over equations (each LHS variable) ----------------------
    for (int eq = 0; eq < ny; eq++) {

      arma::vec u = Yh.col(eq) - Beta(0, eq) - Xh * Beta.submat(1, eq, kr - 1, eq);

      // ------------------------------------------------------------------
      // Full Newey-West HAC sandwich (matches Cesa-Bianchi's VAR Toolbox
      // OLSmodel.m). No residual demeaning (numerical no-op under OLS with
      // an intercept column since X'u = 0).
      //
      //   g_t     = X_reg_t * u_t                          (kr x 1 score)
      //   Gamma_a = sum_{t=a}^{Th-1} g_t * g_{t-a}'          (kr x kr)
      //   G       = w_0*Gamma_0 + sum_{a=1}^{nwL} w_a*(Gamma_a + Gamma_a')
      //   w_a     = (nwL + 1 - a) / (nwL + 1)                (Bartlett)
      //   V_h     = (X'X)^{-1} * G * (X'X)^{-1}
      // ------------------------------------------------------------------
      arma::mat XU = Xreg;      // Th x kr
      XU.each_col() %= u;       // row t scaled by u(t)
      const arma::mat S = XU.t();  // kr x Th, S.col(t) = X_reg_t * u_t

      double var_sc;
      if (store_full) {
        // Need every diagonal element of V, so form the full sandwich.
        arma::mat G(kr, kr, arma::fill::zeros);
        for (int a = 0; a <= nwL; a++) {
          const double w = static_cast<double>(nwL + 1 - a) /
                           static_cast<double>(nwL + 1);        // Bartlett
          const arma::mat Gamma_a = S.cols(a, Th - 1) * S.cols(0, Th - 1 - a).t();
          if (a == 0) G += w * Gamma_a;
          else        G += w * (Gamma_a + Gamma_a.t());
        }
        const arma::mat V = invXpX * G * invXpX;
        var_sc = V(sc, sc);
        for (int j = 0; j < kr; j++) {
          Se_h(j, eq) = std::sqrt(std::max(0.0, V(j, j)));
        }
      } else {
        // Fast path: we only need V(sc, sc) = invXpX(sc, :) * G * invXpX(:, sc).
        // Rewrite as sum_a w_a * (invXpX_sc * Gamma_a * invXpX_sc') using
        // a = S' invXpX_sc'  (Th x 1) so invXpX_sc * Gamma_a * invXpX_sc'
        // = sum_{t=a}^{Th-1} a(t) * a(t - a).
        const arma::vec a_vec = S.t() * invXpX_sc.t();   // Th x 1
        double G_sc = arma::dot(a_vec, a_vec);           // Gamma_0
        for (int lag = 1; lag <= nwL; lag++) {
          const double w = static_cast<double>(nwL + 1 - lag) /
                           static_cast<double>(nwL + 1);
          double gamma_lag = 0.0;
          for (int t = lag; t < Th; t++) gamma_lag += a_vec(t) * a_vec(t - lag);
          G_sc += 2.0 * w * gamma_lag;
        }
        var_sc = G_sc;
      }

      const double irf_h = Beta(sc, eq);
      const double se_sh = std::sqrt(std::max(0.0, var_sc));

      out.irfs(h, eq)       = irf_h;
      out.irfs_se(h, eq)    = se_sh;
      out.irfs_upper(h, eq) = irf_h + z_crit * se_sh;
      out.irfs_lower(h, eq) = irf_h - z_crit * se_sh;

    } // end equation loop

    if (store_full) {
      out.betas[h] = Beta;
      out.ses[h]   = Se_h;
    }

  } // end horizon loop

  return out;
}


// =======================================================================
// R-callable wrapper — exported to R via Rcpp
// =======================================================================

//' Classical Local Projections — C++ Engine
//'
//' Low-level engine called by \code{fLP()}. Prefer the high-level wrapper
//' which provides formula-based syntax and macro expansion.
//'
//' @param Y Numeric matrix (T x n_y). LHS variable(s).
//' @param X Numeric matrix (T x k). RHS variables — no constant, no automatic
//'   lags. Pass exactly what you want on the RHS.
//' @param H Integer. Maximum horizon (projections run for h = 0, 1, ..., H).
//' @param shock_col Integer (0-indexed). Which column of \code{X} is the shock.
//' @param conf_level Numeric. Confidence level for bands, e.g. 0.90.
//' @param nw_lags_base Integer. Base Newey-West bandwidth. At horizon h the
//'   effective bandwidth is \code{max(nw_lags_base + h + nw_offset, 0)}
//'   (Bartlett kernel).
//' @param store_full Logical. If \code{TRUE}, return full coefficient and
//'   standard-error matrices for every horizon.
//' @param cumulative Logical. If \code{TRUE}, regress
//'   \eqn{y_{t+h} - y_{t-1}} on the RHS (cumulative response). Default
//'   \code{FALSE}.
//' @param n_threads Integer. OpenMP threads for the horizon loop. 0 = all
//'   cores minus one. Default 0.
//' @param nw_offset Integer. Shift applied to the NW bandwidth rule
//'   (effective bandwidth = \code{nw_lags_base + h + nw_offset}, floored at
//'   0). Default \code{1L}, the Miranda-Agrippino & Ricco convention
//'   (bandwidth grows as h+1). Set to \code{0L} to reproduce the classic
//'   Jorda (2005) rule, bandwidth = h. \code{fLP}'s \code{h} is 0-indexed
//'   (h=0 is the impact horizon), which makes this the exact equivalent of
//'   Cesa-Bianchi's VAR Toolbox \code{LPmodel.m}
//'   (\code{OLSmodel(Y,X,0,hh-1)} with 1-indexed \code{hh}), e.g. for the
//'   Jorda-Taylor (2025) replication.
//'
//' @return A named list with \code{irfs}, \code{irfs_upper}, \code{irfs_lower},
//'   \code{irfs_se} (raw, unscaled SE of the shock coefficient — lets
//'   callers build confidence bands at any level without re-running the
//'   engine), and optionally \code{betas} and \code{ses}.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List fLP_cpp(
    const arma::mat& Y,
    const arma::mat& X,
    int       H,
    int       shock_col,
    double    conf_level,
    int       nw_lags_base,
    bool      store_full   = false,
    bool      cumulative   = false,
    int       n_threads    = 0,
    int       nw_offset    = 1,
    bool      verbose      = false
) {
  if (Y.n_rows != X.n_rows) {
    Rcpp::stop("fLP_cpp: Y and X must have the same number of rows.");
  }
  if (Y.n_cols == 0) {
    Rcpp::stop("fLP_cpp: Y must have at least one column.");
  }
  if (X.n_cols == 0) {
    Rcpp::stop("fLP_cpp: X must have at least one column.");
  }
  if (H < 0) {
    Rcpp::stop("fLP_cpp: H must be non-negative.");
  }
  if (conf_level <= 0.0 || conf_level >= 1.0) {
    Rcpp::stop("fLP_cpp: conf_level must be in (0, 1).");
  }
  if (nw_lags_base < 0) {
    Rcpp::stop("fLP_cpp: nw_lags_base must be non-negative.");
  }
  // Boundary safety for direct .Call users — R wrappers drop NA via
  // complete.cases but +Inf / -Inf slip through.
  if (!Y.is_finite() || !X.is_finite())
    Rcpp::stop("fLP_cpp: Y and X must be finite (no NA/NaN/Inf).");

  LPResult res = fLP_internal(Y, X, H, shock_col, conf_level, nw_lags_base,
                              store_full, cumulative, n_threads, nw_offset,
                              verbose);

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("irfs")       = res.irfs,
    Rcpp::Named("irfs_upper") = res.irfs_upper,
    Rcpp::Named("irfs_lower") = res.irfs_lower,
    Rcpp::Named("irfs_se")    = res.irfs_se
  );

  if (store_full) {
    Rcpp::List betas_list(res.betas.size());
    Rcpp::List ses_list(res.ses.size());
    for (size_t h = 0; h < res.betas.size(); ++h) {
      betas_list[h] = res.betas[h];
      ses_list[h]   = res.ses[h];
    }
    out["betas"] = betas_list;
    out["ses"]   = ses_list;
  }

  return out;
}
