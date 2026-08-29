// -----------------------------------------------------------------------
// fLPIV.cpp — Local Projections with an External Instrument (LP-IV, C++ engine)
//
// See fLPIV.h for the algorithm; the horizon loop is parallelised with OpenMP.
//
// Author: Dr. Muhsin Ciftci
// -----------------------------------------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include "fLPIV.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <cstdio>

using namespace arma;

// Compute FWL residuals of M on [1, C_h] (constant + C_h).
// Cache-friendly: builds Cext once outside this helper if used repeatedly.
static inline arma::mat fwl_residuals(const arma::mat& M,
                                      const arma::mat& Cext,
                                      const arma::mat& CtC_inv) {
  // beta = (Cext' Cext)^{-1} Cext' M ; resid = M - Cext * beta
  return M - Cext * (CtC_inv * (Cext.t() * M));
}

// =======================================================================
// Internal engine
// =======================================================================

LPIVResult fLPIV_internal(
    const arma::mat& Y,
    const arma::vec& D,
    const arma::mat& Z,
    const arma::mat& C,
    int              H,
    double           conf_level,
    int              nw_lags_iv,
    bool             cumulative,
    int              n_threads,
    bool             verbose,
    const arma::mat& Y_pre
) {
  const int T   = static_cast<int>(Y.n_rows);
  const int ny  = static_cast<int>(Y.n_cols);
  const int nz  = static_cast<int>(Z.n_cols);
  const int nc  = static_cast<int>(C.n_cols);
  const int kc  = nc + 1;                     // controls + intercept

  // Cumulative long-difference base — see fLPIV.h. With a 1 x n_y Y_pre every
  // row of Y is an estimable LHS date (t0 = 0); without it row 0 of Y is spent
  // as the base (t0 = 1) and one observation is lost.
  if (cumulative && Y_pre.n_rows > 0 &&
      (Y_pre.n_rows != 1 || static_cast<int>(Y_pre.n_cols) != ny)) {
    Rcpp::stop("fLPIV: Y_pre must be a 1 x ncol(Y) matrix, or empty.");
  }
  const bool has_pre = cumulative && (Y_pre.n_rows == 1);
  const int  t0      = (cumulative && !has_pre) ? 1 : 0;

  // Alignment with the R wrapper's min_obs guard: need at least
  // kc + nz + a couple of df at the largest horizon.
  if (T <= H + t0 || T - H - t0 <= kc + nz) {
    Rcpp::stop("fLPIV: not enough observations for the largest horizon (need > kc + nz rows).");
  }
  if (D.n_rows != Y.n_rows || Z.n_rows != Y.n_rows ||
      (nc > 0 && C.n_rows != Y.n_rows)) {
    Rcpp::stop("fLPIV: all input matrices must share the same number of rows.");
  }
  if (nz < 1) {
    Rcpp::stop("fLPIV: instrument matrix Z must have at least one column.");
  }

  const double z_crit = R::qnorm(0.5 * (1.0 + conf_level), 0.0, 1.0, 1, 0);

  LPIVResult out;
  out.irfs       = arma::mat(H + 1, ny, arma::fill::zeros);
  out.irfs_upper = arma::mat(H + 1, ny, arma::fill::zeros);
  out.irfs_lower = arma::mat(H + 1, ny, arma::fill::zeros);
  out.irfs_se    = arma::mat(H + 1, ny, arma::fill::zeros);
  out.Fstat_fs   = arma::vec(H + 1, arma::fill::zeros);
  out.rsqr_fs    = arma::vec(H + 1, arma::fill::zeros);

  // Lagged outcome for the cumulative base (horizon-invariant). Row 0 comes
  // from Y_pre when supplied; otherwise it is never read (t0 = 1).
  arma::mat Ylag;
  if (cumulative) {
    Ylag.set_size(T, ny);
    if (has_pre) Ylag.row(0) = Y_pre.row(0);
    else         Ylag.row(0).zeros();
    if (T > 1) Ylag.rows(1, T - 1) = Y.rows(0, T - 2);
  }

  // Per-horizon rank flags (see LPIVResult). Each thread writes only its own
  // element, so no synchronisation is needed; Rcpp::stop() must not be called
  // from inside an OpenMP region, so the caller raises the diagnostic.
  std::vector<char> rank_ok_c(H + 1, 1);   // controls [1, C_h]
  std::vector<char> rank_ok_z(H + 1, 1);   // instruments Z_r

  int actual_threads = 1;
#ifdef _OPENMP
  actual_threads = (n_threads <= 0)
                       ? std::max(1, omp_get_max_threads())
                       : n_threads;
  if (verbose) {
    std::printf("fLPIV: using %d thread(s) for parallel horizon loop...\n",
                actual_threads);
  }
#else
  (void) n_threads;
  if (verbose) {
    std::printf("fLPIV: OpenMP not available. Running single-threaded.\n");
  }
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(actual_threads)
#endif
  for (int h = 0; h <= H; h++) {

    // ---- align data for horizon h --------------------------------------
    arma::mat Yh;
    arma::vec Dh;
    arma::mat Zh;

    // Estimable rows are t = t0, ..., T - 1 - h  →  T - h - t0 rows.
    const int t_last = T - 1 - h;
    if (cumulative) {
      Yh = Y.rows(t0 + h, T - 1) - Ylag.rows(t0, t_last);
    } else {
      Yh = Y.rows(h, T - 1);
    }
    Dh = D.subvec(t0, t_last);
    Zh = Z.rows(t0, t_last);
    const int Th = static_cast<int>(Yh.n_rows);

    // ---- FWL residuals ------------------------------------------------
    // When nc == 0 the projection onto [1, C_h] collapses to a projection
    // onto a constant, i.e. mean-centering. Skip forming Cext / CtC_inv
    // entirely in that case — the common case for LP-IV specifications.
    arma::mat Yr;
    arma::mat Zr;
    arma::vec Dr;
    if (nc == 0) {
      Yr = Yh.each_row() - arma::mean(Yh, 0);
      Zr = Zh.each_row() - arma::mean(Zh, 0);
      Dr = Dh - arma::mean(Dh);
    } else {
      // C's horizon block goes straight into Cext; slicing it into a
      // separate Ch first would copy the same Th x nc block twice.
      arma::mat Cext(Th, kc);
      Cext.col(0).ones();
      Cext.cols(1, kc - 1) = C.rows(t0, t_last);

      // FWL via QR of [1, C_h], mirroring LPmodel.m's backslash
      // (r_y_h = Y - C_h*(C_h\Y)): MATLAB's mldivide on an overdetermined
      // system is a QR least-squares solve, so C'C is never formed. With
      // [1, C_h] = Q R the residual maker is simply I - Q Q', which needs
      // no inverse at all.
      arma::mat Qc, Rc;
      bool ok_c = arma::qr_econ(Qc, Rc, Cext);
      if (ok_c) {
        const arma::vec rd = arma::abs(Rc.diag());
        const double    rm = rd.max();
        ok_c = (rm > 0.0) &&
               (rd.min() > rm * std::max(Th, kc) * arma::datum::eps);
      }
      if (ok_c) {
        Yr = Yh - Qc * (Qc.t() * Yh);
        Zr = Zh - Qc * (Qc.t() * Zh);
        Dr = Dh - Qc * (Qc.t() * Dh);
      } else {
        rank_ok_c[h] = 0;
        const arma::mat CtC_inv = arma::pinv(Cext.t() * Cext);
        Yr = fwl_residuals(Yh, Cext, CtC_inv);
        Zr = fwl_residuals(Zh, Cext, CtC_inv);
        Dr = fwl_residuals(arma::mat(Dh), Cext, CtC_inv);
      }
    }

    // ---- first stage: project D_r onto the instrument space ------------
    // LPmodel.m writes this as D_hat_h = R_Z_h*(R_Z_h\r_d_h), again a QR
    // least-squares solve. With Z_r = Q_z R_z the fitted value is Q_z Q_z' D_r,
    // so Z_r'Z_r — whose condition number is the square of Z_r's — is never
    // formed. This matters most here: weak/collinear instrument sets are
    // exactly where the normal equations lose the most digits.
    arma::mat Qz, Rz;
    bool ok_z = arma::qr_econ(Qz, Rz, Zr);
    if (ok_z) {
      const arma::vec rd = arma::abs(Rz.diag());
      const double    rm = rd.max();
      ok_z = (rm > 0.0) &&
             (rd.min() > rm * std::max(Th, nz) * arma::datum::eps);
    }
    arma::vec Dhat;
    if (ok_z) {
      Dhat = Qz * (Qz.t() * Dr);                             // Th x 1
    } else {
      rank_ok_z[h] = 0;
      Dhat = Zr * (arma::pinv(Zr.t() * Zr) * (Zr.t() * Dr));
    }
    const arma::vec eps_fs = Dr - Dhat;                       // first-stage resid

    // First-stage F-stat and R^2. Dr and Zr are FWL-residualised on
    // [1, C_h] (kc columns), so under H0 the correct residual degrees
    // of freedom for the full regression D ~ [1, C, Z] is Th - kc - nz,
    // NOT Th - nz. Using the latter overstates the denominator df and
    // inflates F, materially when nc is large.
    //
    // NOTE: this is a deliberate divergence from LPmodel.m, which reports
    // ((RSS_res-RSS_unr)/k_z) / (RSS_unr/(n_h - k_z)) — i.e. it ignores the
    // df consumed by the residualised controls. The value below is the
    // textbook exclusion F for D ~ [1, C, Z]; it will not equal the toolbox's
    // reported Fstat_fs. Point estimates, SEs and bands are unaffected.
    //
    // kc and nz are the true ranks here because the caller stops on any
    // rank_ok_* flag, so df2 never overstates the available degrees of
    // freedom in a result that is actually returned.
    const double rss_res = arma::dot(Dr, Dr);          // total SS of Dr
    const double rss_unr = arma::dot(eps_fs, eps_fs);  // resid SS after Zr
    const double kz  = static_cast<double>(nz);
    const int    df2 = Th - kc - nz;
    const double F_h = (df2 > 0 && rss_unr > 0.0)
        ? ((rss_res - rss_unr) / kz) / (rss_unr / static_cast<double>(df2))
        : arma::datum::nan;
    const double R2_fs = (rss_res > 0.0) ? 1.0 - rss_unr / rss_res : 0.0;
    out.Fstat_fs(h) = F_h;
    out.rsqr_fs(h)  = R2_fs;

    // Denominator: D_hat' Dr = D_hat' D_hat (since Zr projects Dr -> D_hat).
    // Guard against numerically-zero denominators with a scale-aware
    // tolerance rather than exact-zero equality — a truly weak instrument
    // will otherwise silently return `0` estimates and SEs.
    const double denom_h  = arma::dot(Dhat, Dr);
    const double denom_tol = std::max(1e-12,
                                      arma::norm(Dhat) * arma::norm(Dr) * 1e-10);
    const bool   denom_ok = std::fabs(denom_h) > denom_tol;

    // ---- Newey-West bandwidth for the score g_t -----------------------
    // nw_lags_iv > 0 : fixed bandwidth = nw_lags_iv
    // nw_lags_iv == 0: horizon-varying bandwidth = h (Miranda hh-1 with 0-indexed h)
    int bw = (nw_lags_iv > 0) ? nw_lags_iv : h;
    bw = std::min(bw, Th - 1);
    if (bw < 0) bw = 0;

    // ---- loop over outcome columns ------------------------------------
    for (int eq = 0; eq < ny; eq++) {

      const arma::vec Yr_eq = Yr.col(eq);
      const double num_h    = arma::dot(Dhat, Yr_eq);
      const double b_h      = denom_ok ? num_h / denom_h : arma::datum::nan;
      const arma::vec eps   = Yr_eq - b_h * Dr;

      // Score and its demeaned version
      arma::vec g = Dhat % eps;
      g -= arma::mean(g);

      // Bartlett HAC long-run variance of g
      double LRV = arma::dot(g, g);   // Gamma_0
      for (int j = 1; j <= bw; j++) {
        const double w = 1.0 - static_cast<double>(j) /
                               static_cast<double>(bw + 1);
        double gamma_j = 0.0;
        for (int t = j; t < Th; t++) {
          gamma_j += g(t) * g(t - j);
        }
        LRV += 2.0 * w * gamma_j;
      }

      const double se_h = denom_ok
          ? std::sqrt(std::max(0.0, LRV)) / std::fabs(denom_h)
          : arma::datum::nan;

      out.irfs(h, eq)       = b_h;
      out.irfs_se(h, eq)    = se_h;
      out.irfs_upper(h, eq) = b_h + z_crit * se_h;
      out.irfs_lower(h, eq) = b_h - z_crit * se_h;
    }

  } // end horizon loop

  // Raise the rank diagnostic outside the parallel region (see rank_ok_*).
  for (int h = 0; h <= H; h++) {
    if (!rank_ok_c[h] || !rank_ok_z[h]) {
      out.rank_deficient  = true;
      out.rank_fail_h     = h;
      out.rank_fail_is_iv = (rank_ok_z[h] == 0);
      break;
    }
  }

  return out;
}


// =======================================================================
// R-callable wrapper
// =======================================================================

//' Local Projections with External Instrument — C++ Engine
//'
//' Low-level engine called by \code{fLPIV()}. Prefer the high-level wrapper.
//'
//' @param Y Numeric matrix (T x n_y). Outcome variable(s).
//' @param D Numeric vector (T). Endogenous treatment.
//' @param Z Numeric matrix (T x n_z). External instrument set — build any
//'   lags externally (e.g. via \code{l(z, 0:k)} in the wrapper).
//' @param C Numeric matrix (T x n_c). Exogenous controls; NO intercept (the
//'   engine adds one internally). May be empty (0 columns).
//' @param H Integer. Maximum horizon (h = 0, 1, ..., H).
//' @param conf_level Numeric. Confidence level for bands, e.g. 0.90.
//' @param nw_lags_iv Integer. Newey-West Bartlett bandwidth for the IV score.
//'   \code{> 0}: fixed bandwidth (Stata's \code{vce(hac nw <k>)});
//'   \code{== 0}: horizon-varying bandwidth = h (just-identified LP rule).
//' @param cumulative Logical. \code{TRUE} regresses \eqn{y_{t+h} - y_{t-1}};
//'   \code{FALSE} regresses the level \eqn{y_{t+h}}.
//' @param n_threads Integer. OpenMP threads. \code{0} = all available cores.
//' @param verbose Logical. Print the thread count. Default \code{FALSE}.
//' @param Y_pre Numeric matrix (1 x n_y) or \code{NULL}. Cumulative LP only:
//'   the outcome at the date immediately BEFORE row 1 of \code{Y}, used as the
//'   long-difference base for the first estimable date. Supplying it keeps
//'   every row of \code{Y} estimable and matches the \code{endo_lag1}
//'   alignment of Cesa-Bianchi's \code{LPmodel.m}. With \code{NULL} (default)
//'   row 1 of \code{Y} is consumed as the base and one observation is lost.
//'
//' @return A named list with \code{irfs}, \code{irfs_upper}, \code{irfs_lower},
//'   \code{irfs_se}, \code{Fstat_fs} (first-stage F per horizon) and
//'   \code{rsqr_fs} (first-stage R^2 per horizon). \code{Fstat_fs} uses
//'   \eqn{T_h - k_c - n_z} residual degrees of freedom, which differs by
//'   construction from the value reported by \code{LPmodel.m}.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List fLPIV_cpp(
    const arma::mat& Y,
    const arma::vec& D,
    const arma::mat& Z,
    const arma::mat& C,
    int       H,
    double    conf_level,
    int       nw_lags_iv,
    bool      cumulative = false,
    int       n_threads  = 0,
    bool      verbose    = false,
    Rcpp::Nullable<arma::mat> Y_pre = R_NilValue
) {
  if (Y.n_rows == 0) Rcpp::stop("fLPIV_cpp: Y is empty.");
  if (H < 0) Rcpp::stop("fLPIV_cpp: H must be non-negative.");
  if (conf_level <= 0.0 || conf_level >= 1.0)
    Rcpp::stop("fLPIV_cpp: conf_level must be in (0, 1).");
  if (nw_lags_iv < 0)
    Rcpp::stop("fLPIV_cpp: nw_lags_iv must be non-negative.");
  // Boundary safety for direct .Call users — reject non-finite inputs
  // early rather than letting NaNs propagate through the horizon loop.
  if (!Y.is_finite() || !D.is_finite() || !Z.is_finite() ||
      (C.n_elem > 0 && !C.is_finite()))
    Rcpp::stop("fLPIV_cpp: Y, D, Z, C must be finite (no NA/NaN/Inf).");

  arma::mat Y_pre_mat;
  if (Y_pre.isNotNull()) {
    Y_pre_mat = Rcpp::as<arma::mat>(Y_pre.get());
    if (!Y_pre_mat.is_finite())
      Rcpp::stop("fLPIV_cpp: Y_pre must be finite (no NA/NaN/Inf).");
  }

  LPIVResult res = fLPIV_internal(Y, D, Z, C, H, conf_level, nw_lags_iv,
                                  cumulative, n_threads, verbose, Y_pre_mat);

  if (res.rank_deficient) {
    Rcpp::stop(
      "fLPIV_cpp: the %s are rank deficient at horizon %d. The 2SLS "
      "coefficient is not identified and the first-stage F degrees of freedom "
      "are undefined. Drop the collinear columns from %s.",
      res.rank_fail_is_iv ? "instruments Z" : "controls [1, C]",
      res.rank_fail_h,
      res.rank_fail_is_iv ? "Z (check for duplicated instrument lags)"
                          : "C (check for duplicated lags or a constant column)");
  }

  return Rcpp::List::create(
    Rcpp::Named("irfs")       = res.irfs,
    Rcpp::Named("irfs_upper") = res.irfs_upper,
    Rcpp::Named("irfs_lower") = res.irfs_lower,
    Rcpp::Named("irfs_se")    = res.irfs_se,
    Rcpp::Named("Fstat_fs")   = res.Fstat_fs,
    Rcpp::Named("rsqr_fs")    = res.rsqr_fs
  );
}
