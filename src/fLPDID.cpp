// -----------------------------------------------------------------------
// fLPDID.cpp — LP-DiD (Dube-Girardi-Jorda-Taylor) C++ engine
// See fLPDID.h for the full description. Design mirrors fLPPanel.cpp:
// per-unit time maps for gap-robust shifts, precomputed read-only caches,
// OpenMP parallelization over horizons.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include "fLPDID.h"
#include "panel_shift.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <limits>

using namespace arma;

static inline bool fin(double x) { return std::isfinite(x); }

struct TaskRes {
  double b  = arma::datum::nan;
  double se = arma::datum::nan;
  int    n  = 0;
  int    G  = 0;
  int    dropped = 0;   // collinear control columns omitted (as in Stata)
};

// [[Rcpp::export]]
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
    int  n_threads,
    bool verbose = false
) {
  const uword n = y.n_elem;
  const uword K = X.n_cols;

  // ---- boundary validation (this engine is directly Rcpp-exported) ----
  if (n == 0)
    Rcpp::stop("fLPDID_cpp: y is empty.");
  if (treat.n_elem != n || i_index.n_elem != n ||
      t_index.n_elem != n || cl_index.n_elem != n ||
      X.n_rows != n)
    Rcpp::stop("fLPDID_cpp: y, treat, X, i_index, t_index, cl_index must have the same length.");
  if (pre_window < 0 || post_window < 0)
    Rcpp::stop("fLPDID_cpp: pre_window and post_window must be non-negative.");
  if (ccc < 0 || ccc > 2)
    Rcpp::stop("fLPDID_cpp: ccc must be one of 0, 1, or 2.");
  if (nonabsorbing && ccc > 0 && Lwin < 0)
    Rcpp::stop("fLPDID_cpp: Lwin must be non-negative when nonabsorbing = TRUE and ccc > 0.");

  // =====================================================================
  // 1. Panel index: per-unit row lists and t -> row maps (gap robust) -
  //    shared helper (see src/panel_shift.h).
  // =====================================================================
  std::vector<uword>                          unit_of;
  std::vector<std::vector<uword>>             unit_rows;
  std::vector<std::unordered_map<int, uword>> t2r;
  panel_build_index(i_index, t_index, unit_of, unit_rows, t2r);
  const uword nu = unit_rows.size();

  auto row_at = [&](uword k, int off) -> sword {
    return panel_row_at(k, off, unit_of, t2r, t_index);
  };

  // =====================================================================
  // 2. Caches: treatment change dD_it and baseline b_it
  // =====================================================================
  vec dtr(n,  fill::value(datum::nan));
  vec base(n, fill::value(datum::nan));
  for (uword k = 0; k < n; ++k) {
    sword l = row_at(k, -1);
    if (l >= 0) {
      if (fin(treat(k)) && fin(treat((uword)l)))
        dtr(k) = treat(k) - treat((uword)l);
      if (!pmd) base(k) = y((uword)l);
    }
  }
  if (pmd) {
    // Stata PMD baseline (DGJT replication scripts, e.g.
    //   inst/LPDID/codes/main_text_simulation/first_sim_pmd_LPDiD_estimation.do):
    //
    //   bysort stateid (year) : gen cumulative_y = sum(ln_y)      // sum() treats missing as 0
    //   gen time  = year - start_year + 1                          // start_year = GLOBAL min year
    //   gen aveLY = L.cumulative_y / (time - 1)                    // baseline
    //
    // We reproduce this exactly:
    //   * numerator: running per-unit sum of y, treating NA as 0
    //   * baseline is defined only when the immediately preceding
    //     calendar row exists for the unit (Stata's L.  returns missing
    //     otherwise), i.e. row_at(k, -1) >= 0
    //   * denominator: (t_index(k) - global_t_min), i.e. `time - 1`
    //
    // This matches the published replication numbers on unbalanced /
    // gap-containing panels; the previous observed-past denominator
    // diverged from Stata off the balanced case.
    const int t_min_global = t_index.min();

    vec cum_y(n, fill::value(datum::nan));
    for (uword u = 0; u < nu; ++u) {
      std::vector<uword> rows = unit_rows[u];
      std::sort(rows.begin(), rows.end(),
                [&](uword a, uword b) { return t_index(a) < t_index(b); });
      double s = 0.0;
      for (uword k : rows) {
        if (fin(y(k))) s += y(k);   // Stata sum(): missing -> 0
        cum_y(k) = s;
      }
    }
    for (uword k = 0; k < n; ++k) {
      const sword l = row_at(k, -1);
      const double time_minus_1 = double(t_index(k) - t_min_global);
      if (l >= 0 && fin(cum_y((uword)l)) && time_minus_1 > 0.0) {
        base(k) = cum_y((uword)l) / time_minus_1;
      }
      // else: baseline stays NaN (Stata's aveLY missing case)
    }
  }

  // =====================================================================
  // 3. Clean-control sets (non-absorbing case), Stata missing semantics:
  //    a missing treatment change counts as "no change" (condition OK);
  //    a missing lagged CCS counts as FALSE.
  // =====================================================================
  auto no_change = [&](sword r) -> bool {
    return (r < 0) || !fin(dtr((uword)r)) || std::fabs(dtr((uword)r)) != 1.0;
  };

  std::vector<unsigned char> ccs0;
  std::vector<std::vector<unsigned char>> ccsF;  // [h][row], h = 0..post
  std::vector<std::vector<unsigned char>> ccsm;  // [j][row], j = 1..pre

  if (nonabsorbing && ccc > 0) {
    ccs0.assign(n, 1);
    for (uword k = 0; k < n; ++k)
      for (int j = 1; j <= Lwin; ++j)
        if (!no_change(row_at(k, -j))) { ccs0[k] = 0; break; }

    ccsF.assign((size_t)post_window + 1, std::vector<unsigned char>(n));
    for (uword k = 0; k < n; ++k) ccsF[0][k] = ccs0[k];
    for (int h = 1; h <= post_window; ++h)
      for (uword k = 0; k < n; ++k)
        ccsF[h][k] = ccsF[h - 1][k] && no_change(row_at(k, h));

    ccsm.assign((size_t)pre_window + 1, std::vector<unsigned char>(n));
    if (pre_window >= 1)
      for (uword k = 0; k < n; ++k) ccsm[1][k] = ccs0[k];
    for (int j = 2; j <= pre_window; ++j)
      for (uword k = 0; k < n; ++k) {
        sword l = row_at(k, -1);
        ccsm[j][k] = ccsm[j - 1][k] && (l >= 0) && ccsm[j - 1][(uword)l];
      }
  }

  // =====================================================================
  // 3b. Precomputed global compact time and cluster ids. Every horizon
  //     previously rebuilt an unordered_map<int, uword> from scratch;
  //     with global codes ready the per-horizon remap is a plain array
  //     lookup + a scratch counter of size nT_global / nG_global.
  // =====================================================================
  uword nT_global = 0, nG_global = 0;
  std::vector<uword> t_compact(n), g_compact(n);
  {
    std::unordered_map<int, uword> tm;  tm.reserve(256);
    for (uword k = 0; k < n; ++k) {
      auto it = tm.find(t_index(k));
      if (it == tm.end()) { t_compact[k] = nT_global; tm.emplace(t_index(k), nT_global++); }
      else                { t_compact[k] = it->second; }
    }
  }
  {
    std::unordered_map<int, uword> gm;  gm.reserve(512);
    for (uword k = 0; k < n; ++k) {
      auto it = gm.find(cl_index(k));
      if (it == gm.end()) { g_compact[k] = nG_global; gm.emplace(cl_index(k), nG_global++); }
      else                { g_compact[k] = it->second; }
    }
  }
  const uword NONE = std::numeric_limits<uword>::max();

  // =====================================================================
  // 4. Per-horizon estimation (shared read-only caches; thread safe)
  // =====================================================================
  auto run_task = [&](int hh, bool post) -> TaskRes {
    TaskRes R;

    // ---- 4a. sample selection ----------------------------------------
    std::vector<uword>  rows;  rows.reserve(n / 2);
    std::vector<double> dyv;   dyv.reserve(n / 2);
    std::vector<double> dv;    dv.reserve(n / 2);

    for (uword k = 0; k < n; ++k) {
      const double bs = base(k); if (!fin(bs)) continue;
      const double dk = dtr(k);  if (!fin(dk)) continue;
      const sword  tr_ = row_at(k, post ? hh : -hh);
      if (tr_ < 0) continue;
      const double dy = y((uword)tr_) - bs; if (!fin(dy)) continue;

      const bool treated = (dk == 1.0);
      bool in = false;

      if (!nonabsorbing) {
        if (treated) {
          in = true;
        } else if (dk == 0.0) {
          if (post) {
            sword f = row_at(k, hh);
            in = (f >= 0 && fin(treat((uword)f)) && treat((uword)f) == 0.0);
          } else {
            in = (fin(treat(k)) && treat(k) == 0.0);
          }
        }
      } else {
        if (treated) {
          if (ccc == 0) in = true;
          else          in = post ? (bool)ccs0[k] : (bool)ccsm[hh][k];
        } else if (dk == 0.0 && fin(treat(k)) && treat(k) == 0.0) {
          if (ccc == 0)      in = true;
          else if (post)     in = (ccc == 1) ? (bool)ccs0[k]
                                             : (bool)ccsF[hh][k];
          else               in = (bool)ccsm[hh][k];
        }
      }
      if (!in) continue;

      bool okX = true;
      for (uword c = 0; c < K; ++c) if (!fin(X(k, c))) { okX = false; break; }
      if (!okX) continue;

      rows.push_back(k);
      dyv.push_back(dy);
      dv.push_back(treated ? 1.0 : 0.0);
    }

    const uword m0 = rows.size();
    if (m0 < K + 3) return R;

    // ---- 4b. compact year ids; reweighting; row selection ------------
    // Densify only the years that appear this horizon by scanning
    // t_compact[k] against a scratch remap of size nT_global. Same for
    // clusters. No unordered_map allocations inside the horizon loop.
    std::vector<uword> ymap_scratch(nT_global, NONE);
    std::vector<uword> yid0(m0);
    uword nY0 = 0;
    for (uword s = 0; s < m0; ++s) {
      const uword tg = t_compact[rows[s]];
      if (ymap_scratch[tg] == NONE) ymap_scratch[tg] = nY0++;
      yid0[s] = ymap_scratch[tg];
    }

    std::vector<uword> sel;  sel.reserve(m0);
    std::vector<double> wy;  // per-year weight (reweight case)
    if (reweight) {
      std::vector<double> tot(nY0, 0.0), trt(nY0, 0.0);
      for (uword s = 0; s < m0; ++s) { tot[yid0[s]] += 1.0; trt[yid0[s]] += dv[s]; }
      wy.assign(nY0, datum::nan);
      for (uword v = 0; v < nY0; ++v) {
        const double p = trt[v] / tot[v];
        if (p > 0.0 && p < 1.0) wy[v] = 1.0 / (1.0 - p);
      }
      for (uword s = 0; s < m0; ++s) if (fin(wy[yid0[s]])) sel.push_back(s);
    } else {
      for (uword s = 0; s < m0; ++s) sel.push_back(s);
    }

    // ---- 4b'. drop singleton time cells (reghdfe default) ------------
    // A year cell holding a single selected row demeans to exactly zero in
    // both y and Z, so it moves neither the coefficient nor the cluster
    // score - but it does inflate m in the CR1 finite-sample factor and the
    // reported nobs. reghdfe removes such rows before estimating; match it.
    // One pass suffices: with a single absorbed dimension, removing a
    // size-1 cell cannot shrink any other cell.
    {
      std::vector<uword> ycnt(nY0, 0);
      for (uword s : sel) ycnt[yid0[s]] += 1;
      std::vector<uword> sel_ns; sel_ns.reserve(sel.size());
      for (uword s : sel) if (ycnt[yid0[s]] > 1) sel_ns.push_back(s);
      sel.swap(sel_ns);
    }

    const uword m = sel.size();
    if (m < K + 3) return R;

    // Recompact year & cluster ids over the SELECTED rows (which may be
    // a strict subset when reweight = TRUE drops undefined-weight years).
    // Reset ymap_scratch positions that were touched by yid0 first.
    for (uword s = 0; s < m0; ++s) ymap_scratch[t_compact[rows[s]]] = NONE;

    std::vector<uword> gmap_scratch(nG_global, NONE);
    std::vector<uword> yid(m), gid(m);
    vec w(m, fill::ones);
    uword nY = 0, nG = 0;
    for (uword s = 0; s < m; ++s) {
      const uword s0 = sel[s];
      const uword k  = rows[s0];
      const uword tg = t_compact[k];
      if (ymap_scratch[tg] == NONE) ymap_scratch[tg] = nY++;
      yid[s] = ymap_scratch[tg];
      const uword gg = g_compact[k];
      if (gmap_scratch[gg] == NONE) gmap_scratch[gg] = nG++;
      gid[s] = gmap_scratch[gg];
      if (reweight) w(s) = wy[yid0[s0]];
    }
    if (nG < 2) return R;

    // ---- 4c. build design, weighted within-year demeaning ------------
    const uword P = K + 1;                 // [dD, controls]
    mat Z(m, P);
    vec yv(m);
    for (uword s = 0; s < m; ++s) {
      const uword s0 = sel[s];
      Z(s, 0) = dv[s0];
      const uword k = rows[s0];
      for (uword c = 0; c < K; ++c) Z(s, 1 + c) = X(k, c);
      yv(s) = dyv[s0];
    }

    vec wsum(nY, fill::zeros);
    mat zbar(nY, P, fill::zeros);
    vec ybar(nY, fill::zeros);
    for (uword s = 0; s < m; ++s) {
      wsum(yid[s]) += w(s);
      zbar.row(yid[s]) += w(s) * Z.row(s);
      ybar(yid[s]) += w(s) * yv(s);
    }
    for (uword v = 0; v < nY; ++v) {
      zbar.row(v) /= wsum(v);
      ybar(v)     /= wsum(v);
    }
    for (uword s = 0; s < m; ++s) {
      Z.row(s) -= zbar.row(yid[s]);
      yv(s)    -= ybar(yid[s]);
    }

    // ---- 4d. drop collinear columns (as Stata/reghdfe does), then
    //          weighted OLS + CR1 cluster sandwich ---------------------
    const vec sw = sqrt(w);
    mat Zs = Z.each_col() % sw;
    vec ys = yv % sw;

    const mat A = Zs.t() * Zs;

    // Rank-revealing Cholesky with sequential column dropping. Column 0
    // (the treatment change) is processed first and never dropped unless
    // it is itself degenerate. Later columns whose Schur complement is
    // numerically zero are exactly collinear with the kept set -> omit,
    // mirroring reghdfe's "omitted because of collinearity".
    const double ctol = 1e-10;
    std::vector<uword> keepc; keepc.reserve(P);
    mat Lch(P, P, fill::zeros);
    for (uword j = 0; j < P; ++j) {
      const uword q = keepc.size();
      vec wv(q);
      for (uword r = 0; r < q; ++r) {
        double acc = A(keepc[r], j);
        for (uword c2 = 0; c2 < r; ++c2) acc -= Lch(r, c2) * wv(c2);
        wv(r) = acc / Lch(r, r);
      }
      const double dj = A(j, j) - dot(wv, wv);
      if (!std::isfinite(dj) || A(j, j) <= 0.0 || dj <= ctol * A(j, j)) {
        if (j == 0) return R;          // treatment itself degenerate
        R.dropped += 1;                // collinear control: omit
        continue;
      }
      for (uword c2 = 0; c2 < q; ++c2) Lch(q, c2) = wv(c2);
      Lch(q, q) = std::sqrt(dj);
      keepc.push_back(j);
    }

    const uword Pk = keepc.size();
    const uvec kv  = conv_to<uvec>::from(keepc);
    const mat Zk   = Z.cols(kv);
    const mat Zsk  = (Pk < P) ? mat(Zs.cols(kv)) : Zs;
    // A was already formed for the rank pass; its kept-column submatrix IS
    // Zsk' Zsk, so there is no need to re-accumulate the cross-product.
    const mat Ak   = (Pk < P) ? mat(A.submat(kv, kv)) : A;

    // Only the treatment coefficient and its variance are returned, so the
    // full inverse is unnecessary: q = Ak^{-1} e_0 gives V(0,0) as
    // cadj * q' meat q.
    vec e0(Pk, fill::zeros); e0(0) = 1.0;
    vec coef, q;
    if (!solve(coef, Ak, Zsk.t() * ys, solve_opts::likely_sympd) ||
        !solve(q,    Ak, e0,           solve_opts::likely_sympd)) return R;

    const vec u  = yv - Zk * coef;     // unweighted residuals
    const vec us = u % sw;             // score uses w * z * u = (sw z)(sw u)

    mat S(nG, Pk, fill::zeros);
    for (uword s = 0; s < m; ++s) S.row(gid[s]) += us(s) * Zsk.row(s);
    const vec Sq = S * q;              // meat sandwich needs only S q

    const double Kfull = double(Pk + nY);          // reghdfe fixef.K = "full"
    const double denom = double(m) - Kfull;
    if (denom <= 0.0) return R;
    const double cadj = (double(nG) / double(nG - 1)) *
                        ((double(m) - 1.0) / denom);
    const double v00 = cadj * dot(Sq, Sq);         // q' (S'S) q

    R.b  = coef(0);                    // treatment is keepc[0] == 0
    R.se = std::sqrt(std::max(v00, 0.0));
    R.n  = (int)m;
    R.G  = (int)nG;
    return R;
  };

  // =====================================================================
  // 5. Task list and OpenMP loop over horizons
  // =====================================================================
  const int pre_start = pmd ? 1 : 2;   // -1 is the normalization when !pmd
  std::vector<int>  task_h;
  std::vector<char> task_post;
  for (int h = 0; h <= post_window; ++h) { task_h.push_back(h); task_post.push_back(1); }
  for (int j = pre_start; j <= pre_window; ++j) { task_h.push_back(j); task_post.push_back(0); }
  const int nT = (int)task_h.size();

  std::vector<TaskRes> res(nT);

  int nt = n_threads;
#ifdef _OPENMP
  if (nt <= 0) nt = omp_get_max_threads();
#pragma omp parallel for schedule(dynamic) num_threads(nt)
#endif
  for (int ti = 0; ti < nT; ++ti)
    res[ti] = run_task(task_h[ti], task_post[ti] != 0);

  // =====================================================================
  // 6. Assemble output
  // =====================================================================
  ivec ev(nT);
  vec  b(nT), se(nT);
  ivec nobs(nT), nclust(nT), ndrop(nT);
  for (int ti = 0; ti < nT; ++ti) {
    ev(ti)     = task_post[ti] ? task_h[ti] : -task_h[ti];
    b(ti)      = res[ti].b;
    se(ti)     = res[ti].se;
    nobs(ti)   = res[ti].n;
    nclust(ti) = res[ti].G;
    ndrop(ti)  = res[ti].dropped;
  }

  return Rcpp::List::create(
    Rcpp::Named("event_time") = ev,
    Rcpp::Named("estimate")   = b,
    Rcpp::Named("se")         = se,
    Rcpp::Named("nobs")       = nobs,
    Rcpp::Named("nclust")     = nclust,
    Rcpp::Named("ndrop")      = ndrop
  );
}
