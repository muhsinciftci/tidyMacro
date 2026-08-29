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
    bool             verbose,
    const arma::mat& Y_pre
) {
  // ---------- dimensions -----------------------------------------------
  const int T  = static_cast<int>(Y.n_rows);
  const int ny = static_cast<int>(Y.n_cols);
  const int k  = static_cast<int>(X.n_cols);
  const int kr = k + 1;                          // +1 for constant (prepended)
  const int sc = shock_col + 1;                  // column index in X_reg (0 = const)

  // ---------- cumulative long-difference base ---------------------------
  // Cumulative LP regresses y_{t+h} - y_{t-1}. When Y_pre supplies the
  // outcome at the date immediately before row 0 of Y, every row of Y is an
  // estimable LHS date (t0 = 0) and the first usable observation survives —
  // the alignment used by Cesa-Bianchi's LPmodel.m. Without it, row 0 of Y
  // must be spent as the base and estimation starts at row 1 (t0 = 1).
  if (cumulative && Y_pre.n_rows > 0 &&
      (Y_pre.n_rows != 1 || static_cast<int>(Y_pre.n_cols) != ny)) {
    Rcpp::stop("fLP: Y_pre must be a 1 x ncol(Y) matrix, or empty.");
  }
  const bool has_pre = cumulative && (Y_pre.n_rows == 1);
  const int  t0      = (cumulative && !has_pre) ? 1 : 0;  // first estimable row

  if (T <= H + t0 || T - H - t0 <= kr) {
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

  // ---------- lagged outcome for the cumulative base (horizon-invariant) --
  // Ylag.row(t) holds the outcome one period before Y.row(t). Row 0 comes
  // from Y_pre when supplied; otherwise it is never read (t0 = 1).
  arma::mat Ylag;
  if (cumulative) {
    Ylag.set_size(T, ny);
    if (has_pre) Ylag.row(0) = Y_pre.row(0);
    else         Ylag.row(0).zeros();
    if (T > 1) Ylag.rows(1, T - 1) = Y.rows(0, T - 2);
  }

  // Per-horizon rank flag. Each thread writes only its own element, so no
  // synchronisation is needed; Rcpp::stop() must not be called from inside
  // an OpenMP region, so the diagnostic is raised by the caller afterwards.
  std::vector<char> rank_ok(H + 1, 1);

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
    // Estimable rows are t = t0, ..., T - 1 - h  →  T - h - t0 rows.
    const int t_last = T - 1 - h;
    arma::mat Yh;
    if (cumulative) {
      Yh = Y.rows(t0 + h, T - 1) - Ylag.rows(t0, t_last);
    } else {
      Yh = Y.rows(h, T - 1);
    }
    const int Th = static_cast<int>(Yh.n_rows);

    // ---- design matrix with explicit intercept ------------------------
    // X's horizon block is written straight into Xreg's storage. Slicing it
    // into a separate Xh first would copy the same Th x k block twice, and
    // since the QR path below consumes only Xreg, Xh has no other use.
    arma::mat Xreg(Th, kr);
    Xreg.col(0).ones();
    if (k > 0) Xreg.cols(1, kr - 1) = X.rows(t0, t_last);

    // ---- OLS via QR — all equations simultaneously --------------------
    // Mirrors OLSmodel.m ("Compute (X'X)^{-1} via QR for numerical
    // stability", OLSmodel.m:69-70). X'X is never formed: its condition
    // number is the square of the design's. With [1, X_h] = Q R,
    //   beta      = R^{-1} Q' Y                    (triangular solve)
    //   (X'X)^{-1} = R^{-1} R^{-T}                 (R'R = X'X exactly)
    // so one triangular inverse yields both, and the HAC sandwich below is
    // unchanged. Rank deficiency shows up as a vanishing diagonal entry of
    // R, which is a sharper detector than a failed Cholesky on X'X.
    arma::mat Qx, Rx, invXpX, Beta;
    bool ok = arma::qr_econ(Qx, Rx, Xreg);
    if (ok) {
      const arma::vec rdiag = arma::abs(Rx.diag());
      const double    rmax  = rdiag.max();
      ok = (rmax > 0.0) &&
           (rdiag.min() > rmax * std::max(Th, kr) * arma::datum::eps);
    }
    if (ok) {
      const arma::mat Rinv = arma::inv(arma::trimatu(Rx));
      invXpX = Rinv * Rinv.t();
      Beta   = Rinv * (Qx.t() * Yh);                 // kr x ny
    } else {
      // [1, X_h] is numerically rank deficient at this horizon. Fall back to
      // the pseudo-inverse so the loop completes, but flag it: the resulting
      // coefficient is a min-norm solution and its HAC standard error has no
      // usual interpretation, so the caller stops instead of reporting it.
      rank_ok[h] = 0;
      invXpX = arma::pinv(Xreg.t() * Xreg);
      Beta   = invXpX * (Xreg.t() * Yh);
    }

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

    // Row `sc` of invXpX; needed for both the full-sandwich and the
    // shock-only (fast) branches. Computing it once per horizon avoids
    // redoing the extraction inside the equation loop.
    const arma::rowvec invXpX_sc = invXpX.row(sc);

    // The fast path needs only the projection of X_reg onto row `sc` of
    // (X'X)^{-1}:  a_t = u_t * (X_reg_t . invXpX_sc). The bracketed factor
    // does not depend on the equation, so it is built once per horizon and
    // the equation loop reduces to an elementwise product — no Th x kr score
    // matrices are ever materialised.
    arma::vec proj_sc;   // Th x 1 — fast path only
    if (!store_full) proj_sc = Xreg * invXpX_sc.t();

    // ---- loop over equations (each LHS variable) ----------------------
    for (int eq = 0; eq < ny; eq++) {

      arma::vec u = Yh.col(eq) - Xreg * Beta.col(eq);

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
      double var_sc;
      if (store_full) {
        // Need every diagonal element of V, so form the full sandwich.
        arma::mat XU = Xreg;      // Th x kr
        XU.each_col() %= u;       // row t scaled by u(t)
        const arma::mat S = XU.t();  // kr x Th, S.col(t) = X_reg_t * u_t

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
        // a_t = u_t * (X_reg_t . invXpX_sc) so invXpX_sc * Gamma_a * invXpX_sc'
        // = sum_{t=a}^{Th-1} a(t) * a(t - a).
        const arma::vec a_vec = u % proj_sc;             // Th x 1
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

  // Raise the rank diagnostic outside the parallel region (see rank_ok).
  for (int h = 0; h <= H; h++) {
    if (!rank_ok[h]) {
      out.rank_deficient = true;
      out.rank_fail_h    = h;
      break;
    }
  }

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
//'   available cores. Default 0.
//' @param nw_offset Integer. Shift applied to the NW bandwidth rule
//'   (effective bandwidth = \code{nw_lags_base + h + nw_offset}, floored at
//'   0). Default \code{1L}, the Miranda-Agrippino & Ricco convention
//'   (bandwidth grows as h+1). Set to \code{0L} to reproduce the classic
//'   Jorda (2005) rule, bandwidth = h. \code{fLP}'s \code{h} is 0-indexed
//'   (h=0 is the impact horizon), which makes this the exact equivalent of
//'   Cesa-Bianchi's VAR Toolbox \code{LPmodel.m}
//'   (\code{OLSmodel(Y,X,0,hh-1)} with 1-indexed \code{hh}), e.g. for the
//'   Jorda-Taylor (2025) replication.
//' @param Y_pre Numeric matrix (1 x n_y) or \code{NULL}. Cumulative LP only:
//'   the outcome at the date immediately BEFORE row 1 of \code{Y}, used as
//'   the long-difference base \eqn{y_{t-1}} for the first estimable date.
//'   Supplying it keeps every row of \code{Y} as an estimable LHS date and
//'   matches the \code{endo_lag1} alignment of Cesa-Bianchi's
//'   \code{LPmodel.m}. With \code{NULL} (default) row 1 of \code{Y} is
//'   consumed as the base and one observation is lost — correct only when no
//'   earlier outcome exists. Ignored when \code{cumulative = FALSE}.
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
    bool      verbose      = false,
    Rcpp::Nullable<arma::mat> Y_pre = R_NilValue
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

  arma::mat Y_pre_mat;
  if (Y_pre.isNotNull()) {
    Y_pre_mat = Rcpp::as<arma::mat>(Y_pre.get());
    if (!Y_pre_mat.is_finite())
      Rcpp::stop("fLP_cpp: Y_pre must be finite (no NA/NaN/Inf).");
  }

  LPResult res = fLP_internal(Y, X, H, shock_col, conf_level, nw_lags_base,
                              store_full, cumulative, n_threads, nw_offset,
                              verbose, Y_pre_mat);

  if (res.rank_deficient) {
    Rcpp::stop(
      "fLP_cpp: the regressor matrix [1, X] is rank deficient at horizon %d. "
      "The shock coefficient is not identified and its HAC standard error "
      "has no interpretation. Drop collinear columns from X (check for "
      "duplicated lags, a constant column, or a dummy set that sums to 1).",
      res.rank_fail_h);
  }

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
