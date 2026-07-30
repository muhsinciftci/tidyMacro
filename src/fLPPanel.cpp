// -----------------------------------------------------------------------
// fLPPanel.cpp — Panel Local Projections (C++ engine)
// The horizon loop is parallelized via OpenMP; `y_h` and lag columns are
// precomputed before the parallel region so each horizon works from a
// shared read-only cache.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include "fLPPanel.h"
#include "panel_shift.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <cstdio>

using namespace arma;


// =======================================================================
// Internal helpers (panel index / gap-aware shift live in panel_shift.h,
// shared with fLPDID.cpp).
// =======================================================================


// High-dimensional fixed-effects residualization by iterative alternating
// group-mean demeaning. Matches `regress_HDFE` in the MATLAB / R source.
//
// Convergence uses a RELATIVE Frobenius-norm criterion (scale-invariant);
// max_iter is bumped to match the reference (fixest::demean default 10000).
// The final OLS uses arma::solve with pinv fallback so a rank-deficient
// residualised design does not silently produce garbage coefficients.
static void regress_HDFE(
    arma::mat& y,
    arma::mat& X,
    const arma::imat& FE,
    arma::vec& b_out,
    bool&      converged_out,
    double tol           = 1e-8,
    arma::uword max_iter = 10000
) {
  const arma::uword n_FE = FE.n_cols;
  const arma::uword n_X  = X.n_cols;

  converged_out = true;

  if (n_FE == 0) {
    if (n_X > 0) {
      if (!arma::solve(b_out, X, y, arma::solve_opts::fast)) {
        b_out = arma::pinv(X) * y;
      }
    } else {
      b_out.reset();
    }
    return;
  }

  std::vector<std::vector<std::vector<arma::uword>>> groups(n_FE);
  for (arma::uword f = 0; f < n_FE; ++f) {
    std::unordered_map<int, int> id_map;
    id_map.reserve(64);
    int next_id = 0;
    std::vector<int> compact(FE.n_rows);
    for (arma::uword k = 0; k < FE.n_rows; ++k) {
      int raw = FE(k, f);
      auto it = id_map.find(raw);
      if (it == id_map.end()) {
        id_map.emplace(raw, next_id);
        compact[k] = next_id++;
      } else {
        compact[k] = it->second;
      }
    }
    groups[f].assign(next_id, {});
    for (arma::uword k = 0; k < FE.n_rows; ++k) groups[f][compact[k]].push_back(k);
  }

  // demean_by returns the squared Frobenius norm of the *change* it made
  // to M. Accumulating those in-flight avoids materialising `y_old` and
  // `X_old` (previously two full-matrix copies per iteration).
  auto demean_by = [&](arma::mat& M,
                       const std::vector<std::vector<arma::uword>>& gs)
      -> double {
    double sq_change = 0.0;
    for (const auto& rows : gs) {
      const arma::uword m = rows.size();
      if (m == 0) continue;
      arma::rowvec sums(M.n_cols, arma::fill::zeros);
      for (arma::uword r : rows) sums += M.row(r);
      sums /= static_cast<double>(m);
      // subtracting `sums` from `m` rows changes those rows by `-sums`;
      // the total squared change equals m * dot(sums, sums).
      sq_change += static_cast<double>(m) * arma::dot(sums, sums);
      for (arma::uword r : rows) M.row(r) -= sums;
    }
    return sq_change;
  };

  bool converged = false;
  for (arma::uword iter = 0; iter < max_iter && !converged; ++iter) {
    double sq_diff = 0.0;
    for (arma::uword f = 0; f < n_FE; ++f) {
      sq_diff += demean_by(y, groups[f]);
      if (n_X > 0) sq_diff += demean_by(X, groups[f]);
    }
    const double diff  = std::sqrt(sq_diff);
    const double n1    = arma::norm(y, "fro");
    const double n2    = (n_X > 0) ? arma::norm(X, "fro") : 0.0;
    const double denom = std::sqrt(n1 * n1 + n2 * n2);
    // Relative rule; falls back to absolute if the scale is ~0.
    if (diff <= tol * std::max(denom, 1.0)) converged = true;
  }
  converged_out = converged;

  if (n_X > 0) {
    if (!arma::solve(b_out, X, y, arma::solve_opts::fast)) {
      b_out = arma::pinv(X) * y;
    }
  } else {
    b_out.reset();
  }
}


// =======================================================================
// Main engine
// =======================================================================

LPPanelResult fLPPanel_internal(
    const arma::mat& y_in,
    const arma::mat& s_in,
    const arma::mat& X_in,
    const arma::mat& W_in,
    const arma::imat& FE_in,
    const arma::ivec& i_index,
    const arma::ivec& t_index,
    int  H,
    int  p_max,
    bool small_sample,
    bool cumulative,
    int  n_threads,
    bool verbose
) {
  const arma::uword n_obs = y_in.n_rows;

  if (s_in.n_rows != n_obs || X_in.n_rows != n_obs ||
      i_index.n_elem != n_obs || t_index.n_elem != n_obs) {
    Rcpp::stop("fLPPanel: y, s, X, i_index, t_index must have the same number of rows.");
  }
  if (W_in.n_cols > 0 && W_in.n_rows != n_obs) {
    Rcpp::stop("fLPPanel: W must have the same number of rows as y (or be empty).");
  }
  if (FE_in.n_cols > 0 && FE_in.n_rows != n_obs) {
    Rcpp::stop("fLPPanel: FE must have the same number of rows as y (or be empty).");
  }
  if (H < 0 || p_max < 0) {
    Rcpp::stop("fLPPanel: H and p_max must be non-negative.");
  }

  const arma::vec y = y_in.col(0);
  const arma::mat& s = s_in;
  const arma::mat& X = X_in;

  const arma::uword n_s_cols = s.n_cols;
  const arma::uword n_X_cols = X.n_cols;
  const arma::uword n_s      = n_s_cols * n_X_cols;

  // Interacted regressor sX (columns (jx)*n_s_cols + js)
  arma::mat sX(n_obs, n_s, arma::fill::zeros);
  for (arma::uword jx = 0; jx < n_X_cols; ++jx) {
    for (arma::uword js = 0; js < n_s_cols; ++js) {
      sX.col(jx * n_s_cols + js) = s.col(js) % X.col(jx);
    }
  }

  // Panel index (per-unit rows and t→row maps) — shared helper.
  std::vector<arma::uword>                            unit_of;
  std::vector<std::vector<arma::uword>>               unit_rows;
  std::vector<std::unordered_map<int, arma::uword>>   t_to_row;
  panel_build_index(i_index, t_index, unit_of, unit_rows, t_to_row);

  // --------- Precompute lead / cumulate columns for all horizons -------
  arma::mat y_h_all(n_obs, H + 1);
  {
    arma::vec running(n_obs, arma::fill::zeros);
    for (int h = 0; h <= H; ++h) {
      arma::vec shifted = panel_time_shift(y, unit_rows, t_to_row, t_index, h);
      if (cumulative) {
        running += shifted;
        y_h_all.col(h) = running;
      } else {
        y_h_all.col(h) = shifted;
      }
    }
  }

  // --------- Precompute lag columns of y and sX (once, up to p_max) ----
  arma::mat y_lags(n_obs, p_max);            // column j-1 = lag j of y
  arma::cube sX_lags(n_obs, n_s, p_max);     // slice j-1 = lags of sX cols
  for (int j = 1; j <= p_max; ++j) {
    y_lags.col(j - 1) = panel_time_shift(y, unit_rows, t_to_row, t_index, -j);
    for (arma::uword is = 0; is < n_s; ++is) {
      sX_lags.slice(j - 1).col(is) =
          panel_time_shift(sX.col(is), unit_rows, t_to_row, t_index, -j);
    }
  }

  // Output allocation. CI/pval bands are rebuilt in R from (estimate,
  // SE, df) so they are not computed here — this avoids R::qt / R::pt
  // calls inside the OpenMP region.
  LPPanelResult out;
  out.estimate.set_size(H + 1, n_s);  out.estimate.fill(arma::datum::nan);
  out.SE.set_size(H + 1, n_s);        out.SE.fill(arma::datum::nan);
  out.df.set_size(H + 1, n_s);        out.df.fill(arma::datum::nan);
  out.status.set_size(H + 1);         out.status.zeros();
  out.converged.set_size(H + 1);      out.converged.ones();

  // Threading — unified rule across all LP engines:
  //   n_threads <= 0 -> use all available cores
  //   n_threads  > 0 -> use exactly that many
  int actual_threads = 1;
#ifdef _OPENMP
  actual_threads = (n_threads <= 0)
                       ? std::max(1, omp_get_max_threads())
                       : n_threads;
  if (verbose) {
    std::printf("fLPPanel: using %d thread(s) for parallel horizon loop...\n",
                actual_threads);
  }
#else
  (void) n_threads;
  if (verbose) {
    std::printf("fLPPanel: OpenMP not available. Running single-threaded.\n");
  }
#endif

  // ---- Horizon loop (parallel over h) ----
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(actual_threads)
#endif
  for (int h = 0; h <= H; ++h) {

    const arma::vec y_h = y_h_all.col(h);
    const int p = std::min(h, p_max);

    // Assemble lag panel: p columns of lag(y), then n_s columns per lag j
    // for lag(sX). Match MATLAB reshape ordering: outermost dim is lag j.
    arma::mat W_lag(n_obs, (1 + n_s) * p, arma::fill::zeros);
    for (int j = 1; j <= p; ++j) {
      const arma::uword base = (j - 1) * (1 + n_s);
      W_lag.col(base) = y_lags.col(j - 1);
      for (arma::uword is = 0; is < n_s; ++is) {
        W_lag.col(base + 1 + is) = sX_lags.slice(j - 1).col(is);
      }
    }

    // Keep mask: complete rows over (y_h, sX, W_lag, W)
    arma::uvec keep(n_obs, arma::fill::ones);
    for (arma::uword r = 0; r < n_obs; ++r) {
      if (!std::isfinite(y_h(r))) { keep(r) = 0; continue; }
      for (arma::uword c = 0; c < n_s && keep(r); ++c)
        if (!std::isfinite(sX(r, c))) keep(r) = 0;
      for (arma::uword c = 0; c < W_lag.n_cols && keep(r); ++c)
        if (!std::isfinite(W_lag(r, c))) keep(r) = 0;
      for (arma::uword c = 0; c < W_in.n_cols && keep(r); ++c)
        if (!std::isfinite(W_in(r, c))) keep(r) = 0;
    }
    const arma::uvec keep_idx = arma::find(keep);
    if (keep_idx.n_elem == 0) {
      out.status(h) = 1;   // no complete observations — leave NaN, warn in R
      continue;
    }

    arma::vec y_LP0 = y_h.elem(keep_idx);
    arma::mat X_LP0 = arma::join_rows(sX.rows(keep_idx), W_lag.rows(keep_idx));
    if (W_in.n_cols > 0) X_LP0 = arma::join_rows(X_LP0, W_in.rows(keep_idx));

    arma::imat FE_LP0;
    if (FE_in.n_cols > 0) FE_LP0 = FE_in.rows(keep_idx);
    else                  FE_LP0.set_size(keep_idx.n_elem, 0);

    const arma::uword n_X = X_LP0.n_cols;
    const arma::ivec  t_LP = t_index.elem(keep_idx);

    arma::ivec t_sorted = arma::sort(arma::unique(t_LP));
    const arma::uword T = t_sorted.n_elem;

    // HDFE + OLS
    arma::mat y_LP = y_LP0;
    arma::mat X_LP = X_LP0;
    arma::vec b_LP;
    bool hdfe_ok = true;
    regress_HDFE(y_LP, X_LP, FE_LP0, b_LP, hdfe_ok);
    if (!hdfe_ok) {
      out.converged(h) = 0;
      out.status(h)    = 2;   // flag but continue: coefficients may still be
                              // usable, R wrapper decides how to surface
    }

    for (arma::uword is = 0; is < n_s; ++is) out.estimate(h, is) = b_LP(is);

    // Score aggregation over time
    arma::vec u = y_LP.col(0) - X_LP * b_LP;
    arma::mat Xv_it = X_LP.each_col() % u;

    std::unordered_map<int, arma::uword> t_row;
    t_row.reserve(T * 2);
    for (arma::uword tt = 0; tt < T; ++tt) t_row.emplace(t_sorted(tt), tt);

    arma::mat Xv_t(T, n_X, arma::fill::zeros);
    for (arma::uword r = 0; r < Xv_it.n_rows; ++r) {
      Xv_t.row(t_row[t_LP(r)]) += Xv_it.row(r);
    }

    const arma::mat XX      = X_LP.t() * X_LP;
    const arma::mat XX_pinv = arma::pinv(XX);

    if (small_sample) {

      const arma::mat s_sub = s.rows(keep_idx);
      for (arma::uword is = 0; is < n_s; ++is) {

        const arma::uword js = is % n_s_cols;
        const arma::vec s_LP = s_sub.col(js);

        arma::mat X_t(T, n_X, arma::fill::zeros);
        arma::vec ss_t(T, arma::fill::zeros);
        for (arma::uword r = 0; r < s_LP.n_elem; ++r) {
          const arma::uword tt = t_row[t_LP(r)];
          ss_t(tt) += s_LP(r) * s_LP(r);
          X_t.row(tt) += s_LP(r) * X_LP.row(r);
        }
        for (arma::uword tt = 0; tt < T; ++tt) {
          if (ss_t(tt) != 0.0) X_t.row(tt) /= ss_t(tt);
        }

        if (!X_t.is_finite()) {
          out.SE(h, is) = arma::datum::nan;
          out.df(h, is) = arma::datum::nan;
          continue;
        }

        const arma::mat XtXt      = X_t.t() * X_t;
        const arma::mat XtXt_pinv = arma::pinv(XtXt);
        const arma::mat P0        = arma::eye<arma::mat>(T, T)
                                    - X_t * XtXt_pinv * X_t.t();

        arma::vec diagP0 = P0.diag();
        arma::vec sd     = arma::sqrt(arma::abs(diagP0));
        for (arma::uword tt = 0; tt < T; ++tt)
          if (sd(tt) < 1e-14) sd(tt) = 1e-14;

        arma::mat Xv_scaled = Xv_t.each_col() / sd;
        arma::mat Xv_var    = Xv_scaled.t() * Xv_scaled;
        arma::mat b_var     = XX_pinv * Xv_var * XX_pinv;

        arma::vec XX0_col0 = XtXt_pinv.col(0);
        arma::mat G0(T, T, arma::fill::zeros);
        for (arma::uword tt = 0; tt < T; ++tt) {
          const double scal = arma::as_scalar(X_t.row(tt) * XX0_col0) / sd(tt);
          G0.col(tt) = P0.col(tt) * scal;
        }
        arma::vec lam0 = arma::eig_sym(G0.t() * G0);

        out.SE(h, is) = std::sqrt(std::max(0.0, b_var(is, is)));
        const double s1 = arma::accu(lam0);
        const double s2 = arma::accu(arma::square(lam0));
        out.df(h, is)   = (s2 > 0.0) ? (s1 * s1 / s2) : arma::datum::inf;
      }

    } else {

      const arma::mat Xv_var = Xv_t.t() * Xv_t;
      const arma::mat b_var  = XX_pinv * Xv_var * XX_pinv;
      for (arma::uword is = 0; is < n_s; ++is) {
        out.SE(h, is) = std::sqrt(std::max(0.0, b_var(is, is)));
        out.df(h, is) = arma::datum::inf;
      }
    }

    // CI bands and p-values are computed in R from (estimate, SE, df);
    // keeping them out of the parallel region avoids R::qt / R::pt calls.
  }

  return out;
}


// =======================================================================
// R-callable wrapper — INTERNAL; users call fLPPanel() in R.
// =======================================================================

//' Panel Local Projections — Internal C++ engine
//'
//' Internal engine implementing the estimator of Almuzara & Sancibrián
//' (NY Fed SR 1090, 2024). Users should call \code{fLPPanel()} instead.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
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
    bool small_sample = false,
    bool cumulative   = false,
    int  n_threads    = 0,
    bool verbose      = false
) {
  // Boundary validation — the R wrapper does most of it, but this engine is
  // also exported for direct .Call, so guard against non-finite and shape
  // errors here before entering the parallel region.
  if (y.n_rows == 0) Rcpp::stop("fLPPanel_cpp: y is empty.");
  if (y.n_cols == 0) Rcpp::stop("fLPPanel_cpp: y must have at least one column.");
  if (s.n_rows != y.n_rows || X.n_rows != y.n_rows ||
      (W.n_cols  > 0 && W.n_rows  != y.n_rows) ||
      (FE.n_cols > 0 && FE.n_rows != y.n_rows) ||
      i_index.n_elem != y.n_rows || t_index.n_elem != y.n_rows) {
    Rcpp::stop("fLPPanel_cpp: y, s, X, W, FE, i_index, t_index row counts must match.");
  }
  if (H < 0 || p_max < 0)
    Rcpp::stop("fLPPanel_cpp: H and p_max must be non-negative.");

  LPPanelResult r = fLPPanel_internal(y, s, X, W, FE, i_index, t_index,
                                      H, p_max, small_sample, cumulative,
                                      n_threads, verbose);

  // Surface per-horizon issues collected during the parallel region.
  const arma::uword n_no_obs   = arma::accu(r.status == 1);
  const arma::uword n_hdfe_bad = arma::accu(r.status == 2);
  if (n_no_obs > 0) {
    Rcpp::warning("fLPPanel: %u horizon(s) had no complete observations; "
                  "estimates set to NA.", (unsigned)n_no_obs);
  }
  if (n_hdfe_bad > 0) {
    Rcpp::warning("fLPPanel: HDFE demeaning did not converge for %u horizon(s); "
                  "check for near-collinear fixed effects.", (unsigned)n_hdfe_bad);
  }

  return Rcpp::List::create(
    Rcpp::Named("estimate")  = r.estimate,
    Rcpp::Named("SE")        = r.SE,
    Rcpp::Named("df")        = r.df,
    Rcpp::Named("status")    = r.status,
    Rcpp::Named("converged") = r.converged
  );
}
