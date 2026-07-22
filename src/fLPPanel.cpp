// -----------------------------------------------------------------------
// fLPPanel.cpp — Panel Local Projections (C++ engine)
// The horizon loop is parallelized via OpenMP; `y_h` and lag columns are
// precomputed before the parallel region so each horizon works from a
// shared read-only cache.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include "fLPPanel.h"
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
// Internal helpers
// =======================================================================

// Panel time-shift for a single column: for each unit i, produce
// y1(row) = y0(row') where row' is the unit-i row with time = t + L, or
// NaN if that time is missing. L > 0 → lead ; L < 0 → lag.
static arma::vec panel_time_shift(
    const arma::vec& y0,
    const std::vector<std::vector<arma::uword>>& unit_rows,
    const std::vector<std::unordered_map<int, arma::uword>>& t_to_row,
    const arma::ivec& t_index,
    int L
) {
  const arma::uword n = y0.n_elem;
  arma::vec y1(n, arma::fill::value(arma::datum::nan));

  for (arma::uword u = 0; u < unit_rows.size(); ++u) {
    const auto& rows = unit_rows[u];
    const auto& map  = t_to_row[u];
    for (arma::uword k : rows) {
      auto it = map.find(t_index(k) + L);
      if (it != map.end()) y1(k) = y0(it->second);
    }
  }
  return y1;
}

// Build per-unit row lists and t→row maps used by every panel_time_shift
// invocation. Called once per fLPPanel_internal call.
static void build_panel_index(
    const arma::ivec& i_index,
    const arma::ivec& t_index,
    std::vector<std::vector<arma::uword>>& unit_rows,
    std::vector<std::unordered_map<int, arma::uword>>& t_to_row
) {
  const arma::uword n = i_index.n_elem;
  std::unordered_map<int, arma::uword> unit_ix;
  unit_ix.reserve(256);
  arma::uword next = 0;

  std::vector<arma::uword> unit_of(n);
  for (arma::uword k = 0; k < n; ++k) {
    auto it = unit_ix.find(i_index(k));
    if (it == unit_ix.end()) {
      unit_ix.emplace(i_index(k), next);
      unit_of[k] = next++;
    } else {
      unit_of[k] = it->second;
    }
  }
  unit_rows.assign(next, {});
  for (arma::uword k = 0; k < n; ++k) unit_rows[unit_of[k]].push_back(k);

  t_to_row.assign(next, {});
  for (arma::uword u = 0; u < next; ++u) {
    t_to_row[u].reserve(unit_rows[u].size() * 2);
    for (arma::uword k : unit_rows[u]) t_to_row[u].emplace(t_index(k), k);
  }
}


// High-dimensional fixed-effects residualization by iterative alternating
// group-mean demeaning. Matches `regress_HDFE` in the MATLAB source.
static void regress_HDFE(
    arma::mat& y,
    arma::mat& X,
    const arma::imat& FE,
    arma::vec& b_out,
    double tol      = 1e-8,
    arma::uword max_iter = 500
) {
  const arma::uword n_FE = FE.n_cols;
  const arma::uword n_X  = X.n_cols;

  if (n_FE == 0) {
    if (n_X > 0) b_out = arma::solve(X, y, arma::solve_opts::fast);
    else         b_out.reset();
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

  auto demean_by = [&](arma::mat& M,
                       const std::vector<std::vector<arma::uword>>& gs) {
    for (const auto& rows : gs) {
      const arma::uword m = rows.size();
      if (m == 0) continue;
      arma::rowvec sums(M.n_cols, arma::fill::zeros);
      for (arma::uword r : rows) sums += M.row(r);
      sums /= static_cast<double>(m);
      for (arma::uword r : rows) M.row(r) -= sums;
    }
  };

  bool converged = false;
  for (arma::uword iter = 0; iter < max_iter && !converged; ++iter) {
    arma::mat y_old = y;
    arma::mat X_old = X;
    for (arma::uword f = 0; f < n_FE; ++f) {
      demean_by(y, groups[f]);
      if (n_X > 0) demean_by(X, groups[f]);
    }
    const double d1 = arma::norm(y - y_old, "fro");
    const double d2 = (n_X > 0) ? arma::norm(X - X_old, "fro") : 0.0;
    if (std::sqrt(d1 * d1 + d2 * d2) < tol) converged = true;
  }

  if (n_X > 0) b_out = arma::solve(X, y, arma::solve_opts::fast);
  else         b_out.reset();
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
    int  n_threads
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

  // Panel index (per-unit rows and t→row maps)
  std::vector<std::vector<arma::uword>> unit_rows;
  std::vector<std::unordered_map<int, arma::uword>> t_to_row;
  build_panel_index(i_index, t_index, unit_rows, t_to_row);

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

  // Output allocation
  LPPanelResult out;
  out.estimate.set_size(H + 1, n_s); out.estimate.zeros();
  out.SE.set_size(H + 1, n_s);       out.SE.zeros();
  out.df.set_size(H + 1, n_s);       out.df.zeros();
  out.CI90.set_size(H + 1, n_s, 2);  out.CI90.zeros();
  out.CI95.set_size(H + 1, n_s, 2);  out.CI95.zeros();
  out.CI99.set_size(H + 1, n_s, 2);  out.CI99.zeros();
  out.pval.set_size(H + 1, n_s);     out.pval.zeros();

  // Threading
  int actual_threads = 1;
#ifdef _OPENMP
  actual_threads = (n_threads == 0)
                       ? std::max(1, omp_get_max_threads() - 1)
                       : n_threads;
  omp_set_num_threads(actual_threads);
  std::printf("fLPPanel: using %d thread(s) for parallel horizon loop...\n",
              actual_threads);
#else
  (void) n_threads;
  std::printf("fLPPanel: OpenMP not available. Running single-threaded.\n");
#endif

  // ---- Horizon loop (parallel over h) ----
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
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
      Rcpp::stop("fLPPanel: no complete observations at horizon %d.", h);
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
    regress_HDFE(y_LP, X_LP, FE_LP0, b_LP);

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

    // Fixed 90/95/99 CIs and p-values
    for (arma::uword is = 0; is < n_s; ++is) {
      const double est = out.estimate(h, is);
      const double se  = out.SE(h, is);
      const double df  = out.df(h, is);

      const double q90 = R::qt(0.95, df, 1, 0);
      const double q95 = R::qt(0.975, df, 1, 0);
      const double q99 = R::qt(0.995, df, 1, 0);

      out.CI90(h, is, 0) = est - q90 * se;
      out.CI90(h, is, 1) = est + q90 * se;
      out.CI95(h, is, 0) = est - q95 * se;
      out.CI95(h, is, 1) = est + q95 * se;
      out.CI99(h, is, 0) = est - q99 * se;
      out.CI99(h, is, 1) = est + q99 * se;

      const double t_stat = std::abs(est / se);
      out.pval(h, is) = 2.0 * (1.0 - R::pt(t_stat, df, 1, 0));
    }
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
    int  n_threads    = 0
) {
  LPPanelResult r = fLPPanel_internal(y, s, X, W, FE, i_index, t_index,
                                      H, p_max, small_sample, cumulative,
                                      n_threads);

  return Rcpp::List::create(
    Rcpp::Named("estimate") = r.estimate,
    Rcpp::Named("SE")       = r.SE,
    Rcpp::Named("df")       = r.df,
    Rcpp::Named("CI90")     = r.CI90,
    Rcpp::Named("CI95")     = r.CI95,
    Rcpp::Named("CI99")     = r.CI99,
    Rcpp::Named("pval")     = r.pval
  );
}
