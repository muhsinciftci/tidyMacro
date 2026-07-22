// -----------------------------------------------------------------------
// fLPDID.cpp — LP-DiD (Dube-Girardi-Jorda-Taylor) C++ engine
// See fLPDID.h for the full description. Design mirrors fLPPanel.cpp:
// per-unit time maps for gap-robust shifts, precomputed read-only caches,
// OpenMP parallelization over horizons.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include "fLPDID.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <cmath>

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
    int  n_threads
) {
  const uword n = y.n_elem;
  const uword K = X.n_cols;

  // =====================================================================
  // 1. Panel index: per-unit row lists and t -> row maps (gap robust)
  // =====================================================================
  std::unordered_map<int, uword> uix;  uix.reserve(512);
  std::vector<uword> unit_of(n);
  uword nu = 0;
  for (uword k = 0; k < n; ++k) {
    auto it = uix.find(i_index(k));
    if (it == uix.end()) { uix.emplace(i_index(k), nu); unit_of[k] = nu++; }
    else                 { unit_of[k] = it->second; }
  }
  std::vector<std::vector<uword>> unit_rows(nu);
  for (uword k = 0; k < n; ++k) unit_rows[unit_of[k]].push_back(k);
  std::vector<std::unordered_map<int, uword>> t2r(nu);
  for (uword u = 0; u < nu; ++u) {
    t2r[u].reserve(unit_rows[u].size() * 2);
    for (uword k : unit_rows[u]) t2r[u].emplace(t_index(k), k);
  }

  auto row_at = [&](uword k, int off) -> sword {
    const auto& m = t2r[unit_of[k]];
    auto it = m.find(t_index(k) + off);
    return (it == m.end()) ? sword(-1) : sword(it->second);
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
    // baseline = running mean of ALL past outcomes within the unit
    for (uword u = 0; u < nu; ++u) {
      std::vector<uword> rows = unit_rows[u];
      std::sort(rows.begin(), rows.end(),
                [&](uword a, uword b) { return t_index(a) < t_index(b); });
      double s = 0.0; uword c = 0;
      for (uword k : rows) {
        base(k) = (c > 0) ? s / double(c) : datum::nan;
        if (fin(y(k))) { s += y(k); ++c; }
      }
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
    std::unordered_map<int, uword> ymap0;
    std::vector<uword> yid0(m0);
    for (uword s = 0; s < m0; ++s) {
      const int tv = t_index(rows[s]);
      auto it = ymap0.find(tv);
      if (it == ymap0.end()) { yid0[s] = ymap0.size(); ymap0.emplace(tv, yid0[s]); }
      else                   { yid0[s] = it->second; }
    }
    const uword nY0 = ymap0.size();

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

    const uword m = sel.size();
    if (m < K + 3) return R;

    // recompact year & cluster ids over the selected rows
    std::unordered_map<int, uword> ymap, gmap;
    std::vector<uword> yid(m), gid(m);
    vec w(m, fill::ones);
    for (uword s = 0; s < m; ++s) {
      const uword s0 = sel[s];
      const uword k  = rows[s0];
      const int tv = t_index(k);
      auto it = ymap.find(tv);
      if (it == ymap.end()) { yid[s] = ymap.size(); ymap.emplace(tv, yid[s]); }
      else                  { yid[s] = it->second; }
      const int gv = cl_index(k);
      auto ig = gmap.find(gv);
      if (ig == gmap.end()) { gid[s] = gmap.size(); gmap.emplace(gv, gid[s]); }
      else                  { gid[s] = ig->second; }
      if (reweight) w(s) = wy[yid0[s0]];
    }
    const uword nY = ymap.size();
    const uword nG = gmap.size();
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
    const mat Ak   = Zsk.t() * Zsk;

    vec coef;
    mat Ainv;
    if (!solve(coef, Ak, Zsk.t() * ys, solve_opts::likely_sympd) ||
        !inv_sympd(Ainv, Ak)) return R;

    const vec u  = yv - Zk * coef;     // unweighted residuals
    const vec us = u % sw;             // score uses w * z * u = (sw z)(sw u)

    mat S(nG, Pk, fill::zeros);
    for (uword s = 0; s < m; ++s) S.row(gid[s]) += us(s) * Zsk.row(s);
    const mat meat = S.t() * S;

    const double Kfull = double(Pk + nY);          // reghdfe fixef.K = "full"
    const double denom = double(m) - Kfull;
    if (denom <= 0.0) return R;
    const double cadj = (double(nG) / double(nG - 1)) *
                        ((double(m) - 1.0) / denom);
    const mat V = cadj * Ainv * meat * Ainv;

    R.b  = coef(0);                    // treatment is keepc[0] == 0
    R.se = std::sqrt(std::max(V(0, 0), 0.0));
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
