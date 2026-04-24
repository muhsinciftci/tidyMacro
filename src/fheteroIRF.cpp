// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include "fheteroIRF.h"
#include "fVAR.h"
#include "fgenerateVARdata.h"
#include "fmbb_var.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// =========================================================================
// Internal helper: proxy-IV structural impact vector
// Matches MATLAB: b1 = [uhat] \ yiIV  (treatment group only)
//   XrIV = ZrIV = proxy(indsR1) - mean(proxy(indsR1))   [centred proxy]
//   yiIV = U(indsR1) - mean(U(indsR1))                  [centred residuals]
//   b1   = (XrIV' XrIV)^{-1} XrIV' yiIV                [OLS, no intercept]
//   then normalise so b1[nvar-1] = scale
// =========================================================================
static arma::vec hetero_iv_b1(const arma::mat& eta,    // T1 x N
                               const arma::mat& Z,      // T1 x 1
                               const arma::ivec& indsR1,
                               int nvar, double scale) {
    int T1 = (int)eta.n_rows;

    // Collect OPEC-month indices
    std::vector<arma::uword> idx_vec;
    for (int t = 0; t < T1; ++t)
        if (indsR1(t) == 1) idx_vec.push_back((arma::uword)t);

    if (idx_vec.empty())
        Rcpp::stop("fheteroIRF: indsR1 contains no treatment observations.");

    arma::uvec idx(idx_vec);

    arma::mat U_opec = eta.rows(idx);          // T_opec x N
    arma::vec z_full = Z.col(0);
    arma::vec z_opec = z_full.elem(idx);       // T_opec x 1

    // Centre (MATLAB: proxy - mean(proxy), U - mean(U))
    z_opec   -= arma::mean(z_opec);
    U_opec.each_row() -= arma::mean(U_opec, 0);

    // b1 = (z'z)^{-1} z'U  →  N-vector
    double zz = arma::dot(z_opec, z_opec);
    if (zz < 1e-14)
        Rcpp::stop("fheteroIRF: proxy has zero variance in treatment months.");

    arma::vec b1 = (U_opec.t() * z_opec) / zz;

    // Normalise: b1[nvar-1] = scale
    double nf = b1(nvar - 1);
    if (std::abs(nf) < 1e-14)
        Rcpp::stop("fheteroIRF: normalization element is numerically zero.");
    b1 *= (scale / nf);

    return b1;
}

// =========================================================================
// fheteroIRF_cpp: point estimate
// =========================================================================
HeteroIRFResult fheteroIRF_cpp(const arma::mat& AL,
                                const arma::mat& eta,
                                const arma::mat& Z,
                                const arma::ivec& indsR1,
                                int p, int hor, int nvar, double scale) {
    int N = (int)AL.n_rows;
    int H = hor + 1;

    arma::vec b1 = hetero_iv_b1(eta, Z, indsR1, nvar, scale);

    // MA representation: C_0 = I, C_h = sum_{j=1}^{min(h,p)} A_j C_{h-j}
    std::vector<arma::mat> C(H);
    C[0] = arma::eye(N, N);
    for (int h = 1; h < H; ++h) {
        C[h].zeros(N, N);
        for (int j = 1; j <= std::min(h, p); ++j)
            C[h] += AL.cols((j - 1) * N, j * N - 1) * C[h - j];
    }

    arma::mat IRF(N, H);
    for (int h = 0; h < H; ++h)
        IRF.col(h) = C[h] * b1;

    HeteroIRFResult res;
    res.IRF = IRF;
    res.b1  = b1;
    return res;
}

// =========================================================================
// fbootstrapHetero_cpp: MBB bootstrap
// Resamples (residuals, Z) jointly in the proxy window — same as
// fbootstrapIV_mbb — but identification uses OPEC months only.
// =========================================================================
BootHeteroResult fbootstrapHetero_cpp(const arma::mat& y,
                                      const VARResult& var_result,
                                      const arma::mat& Z,
                                      const arma::ivec& indsR1,
                                      const arma::ivec& adjustu,
                                      int nboot, int blocksize,
                                      int hor, int nvar, double scale,
                                      double prc, double prc2,
                                      int n_threads) {
    const arma::mat& beta      = var_result.beta;
    const arma::mat& residuals = var_result.residuals;
    const int p = var_result.p;
    const int c = var_result.c;
    const int N = (int)y.n_cols;
    const int H = hor + 1;

    const int u_start = adjustu(0) - 1;
    const int u_end   = adjustu(1) - 1;

    // AL from original VAR (used for point-estimate recentering at end)
    arma::mat AL = beta.rows(c, c + N * p - 1).t();

    // Proxy-window residuals and proxy (for MBB)
    arma::mat Epsstar = residuals.rows(u_start, u_end);  // T1 x N
    // Z is already proxy-window aligned (length = u_end - u_start + 1)

    int actual_threads = 1;
#ifdef _OPENMP
    actual_threads = (n_threads == 0) ?
        std::max(1, omp_get_max_threads() - 1) : n_threads;
    omp_set_num_threads(actual_threads);
    Rprintf("Using %d thread(s) for heteroskedasticity bootstrap...\n",
            actual_threads);
#else
    Rprintf("OpenMP not available. Running in single-threaded mode.\n");
#endif

    arma::cube irf_boot(N, H, nboot, arma::fill::zeros);
    // Track valid bootstrap iterations (quality check: >= 15 treatment obs)
    std::vector<int> valid(nboot, 1);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int b = 0; b < nboot; ++b) {

        // Resample (eps, Z) jointly in proxy window
        MBBVARResult mbb = fmbb_var_cpp(Epsstar, p, blocksize, Z);

        // Pad bootstrap residuals to full sample length
        arma::mat res_boot(residuals.n_rows, N, arma::fill::none);
        if (u_start > 0)
            res_boot.rows(0, u_start - 1) = residuals.rows(0, u_start - 1);
        res_boot.rows(u_start, u_end) = mbb.eps_boot;
        if (u_end < (int)residuals.n_rows - 1)
            res_boot.rows(u_end + 1, residuals.n_rows - 1) =
                residuals.rows(u_end + 1, residuals.n_rows - 1);

        // Bootstrap proxy (proxy window only)
        const arma::mat& Z_boot = mbb.M_boot;

        // Derive treatment indicator from bootstrapped proxy (non-zero = OPEC month)
        int T1b = (int)Z_boot.n_rows;
        arma::ivec indsR1_b(T1b);
        int ntreat = 0;
        for (int t = 0; t < T1b; ++t) {
            indsR1_b(t) = (std::abs(Z_boot(t, 0)) > 1e-14) ? 1 : 0;
            ntreat += indsR1_b(t);
        }

        // Quality check: skip if too few treatment observations
        if (ntreat < 15) {
            valid[b] = 0;
            continue;
        }

        // Regenerate VAR data and re-estimate
        arma::mat yboot = fgenerateVARdata(y, p, c, beta, res_boot);
        VARResult vb    = fVAR_cpp(yboot, p, c, R_NilValue);

        arma::mat eta_b = vb.residuals.rows(u_start, u_end);
        arma::mat AL_b  = vb.beta.rows(c, c + N * p - 1).t();

        // Identify using proxy IV in bootstrap OPEC months
        HeteroIRFResult hr =
            fheteroIRF_cpp(AL_b, eta_b, Z_boot, indsR1_b, p, hor, nvar, scale);
        irf_boot.slice(b) = hr.IRF;
    }

    // Point estimate (for recentering)
    arma::mat eta_pt = residuals.rows(u_start, u_end);
    HeteroIRFResult pt =
        fheteroIRF_cpp(AL, eta_pt, Z, indsR1, p, hor, nvar, scale);
    const arma::mat& irf_pt = pt.IRF;

    // Percentile bounds
    const double up_pct   = 50.0 + prc  * 0.5;
    const double low_pct  = 50.0 - prc  * 0.5;
    const double up_pct2  = 50.0 + prc2 * 0.5;
    const double low_pct2 = 50.0 - prc2 * 0.5;

    arma::mat upper(N, H),   lower(N, H);
    arma::mat upper2(N, H),  lower2(N, H);
    arma::mat meanirf(N, H), medianirf(N, H);

    // Collect valid bootstrap indices
    std::vector<int> valid_idx;
    valid_idx.reserve(nboot);
    for (int b = 0; b < nboot; ++b)
        if (valid[b]) valid_idx.push_back(b);
    int nv = (int)valid_idx.size();

    if (nv < 10)
        Rcpp::stop("fbootstrapHetero: fewer than 10 valid bootstrap iterations. "
                   "Check proxy and adjustu settings.");

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int i = 0; i < N; ++i) {
        for (int h = 0; h < H; ++h) {
            arma::vec bv(nv);
            for (int k = 0; k < nv; ++k)
                bv(k) = irf_boot(i, h, valid_idx[k]);

            arma::vec sv = arma::sort(bv);
            double bmed  = sv((int)(nv * 0.5));

            auto interp = [&](double pct) -> double {
                double id = (pct / 100.0) * (nv - 1);
                int    fl = (int)id;
                double fr = id - fl;
                return (fl < nv - 1) ?
                    sv(fl) * (1 - fr) + sv(fl + 1) * fr : sv(nv - 1);
            };

            // Hall (1992) recentered percentile bands
            upper(i, h)     = interp(up_pct)   - bmed + irf_pt(i, h);
            lower(i, h)     = interp(low_pct)  - bmed + irf_pt(i, h);
            upper2(i, h)    = interp(up_pct2)  - bmed + irf_pt(i, h);
            lower2(i, h)    = interp(low_pct2) - bmed + irf_pt(i, h);
            medianirf(i, h) = bmed;
            meanirf(i, h)   = arma::mean(bv);
        }
    }

    BootHeteroResult res;
    res.upper     = upper;
    res.lower     = lower;
    res.upper2    = upper2;
    res.lower2    = lower2;
    res.meanirf   = meanirf;
    res.medianirf = medianirf;
    res.point     = irf_pt;
    return res;
}

// =========================================================================
// R wrappers
// =========================================================================

//' Heteroskedasticity-Based VAR Identification (Proxy IV, Treatment Months)
//'
//' Identifies a structural shock using the external instrument restricted to
//' high-variance (treatment) months, following the proxy-SVAR approach of
//' Kaenzig (2021). The structural impact vector is estimated via OLS of the
//' centred residuals on the centred proxy, restricted to OPEC announcement
//' months (\code{indsR1 == 1}).
//'
//' @param var_result List returned by \code{fVAR()}.
//' @param Z          Proxy matrix (T1 x 1), aligned to the identification window
//'   defined by \code{adjustu}.
//' @param adjustu    Integer vector \code{c(start, end)} (1-based) giving the
//'   rows of \code{var_result$residuals} corresponding to the proxy window.
//' @param indsR1     Integer 0/1 vector of length \code{adjustu[2]-adjustu[1]+1};
//'   1 = treatment month (e.g. OPEC announcement), 0 otherwise.
//' @param hor        IRF horizon.
//' @param nvar       1-based normalization variable index.
//' @param scale      Shock size for normalization (default 10).
//'
//' @return A list with \code{IRF} (N x hor+1) and \code{b1} (N-vector).
//' @export
// [[Rcpp::export]]
Rcpp::List fheteroIRF(const Rcpp::List& var_result,
                      const arma::mat& Z,
                      const arma::ivec& adjustu,
                      const arma::ivec& indsR1,
                      int hor,
                      int nvar,
                      double scale = 10.0) {

    arma::mat beta      = Rcpp::as<arma::mat>(var_result["beta"]);
    arma::mat residuals = Rcpp::as<arma::mat>(var_result["residuals"]);
    int p = Rcpp::as<int>(var_result["p"]);
    int c = Rcpp::as<int>(var_result["c"]);
    int N = (int)beta.n_cols;

    arma::mat AL  = beta.rows(c, c + N * p - 1).t();
    arma::mat eta = residuals.rows(adjustu(0) - 1, adjustu(1) - 1);

    HeteroIRFResult res = fheteroIRF_cpp(AL, eta, Z, indsR1, p, hor, nvar, scale);

    return Rcpp::List::create(
        Rcpp::Named("IRF") = res.IRF,
        Rcpp::Named("b1")  = res.b1
    );
}

//' MBB Bootstrap for Heteroskedasticity-Based VAR Identification
//'
//' Moving-block bootstrap confidence bands for structural IRFs identified via
//' proxy IV restricted to OPEC announcement months (Kaenzig 2021, treatment
//' group only). Residuals and proxy are resampled jointly.
//'
//' @param y          T x N data matrix.
//' @param var_result List returned by \code{fVAR()}.
//' @param Z          Proxy matrix (T1 x 1), aligned to \code{adjustu}.
//' @param indsR1     Integer 0/1 vector of length T1 (1 = treatment month).
//' @param adjustu    Integer vector \code{c(start, end)} (1-based).
//' @param nboot      Bootstrap replications.
//' @param blocksize  MBB block size (0 = automatic \code{5.03*T^0.25}).
//' @param hor        IRF horizon.
//' @param nvar       1-based normalization variable.
//' @param scale      Shock size.
//' @param prc        Primary confidence level (e.g. 90).
//' @param prc2       Secondary confidence level (e.g. 68).
//' @param n_threads  OpenMP threads (0 = all available - 1).
//'
//' @return A list with \code{upper}, \code{lower}, \code{upper2},
//'   \code{lower2}, \code{meanirf}, \code{medianirf} (each N x hor+1).
//' @export
// [[Rcpp::export]]
Rcpp::List fbootstrapHetero(const arma::mat& y,
                             const Rcpp::List& var_result,
                             const arma::mat& Z,
                             const arma::ivec& indsR1,
                             const arma::ivec& adjustu,
                             int nboot      = 1000,
                             int blocksize  = 0,
                             int hor        = 48,
                             int nvar       = 1,
                             double scale   = 10.0,
                             double prc     = 90.0,
                             double prc2    = 68.0,
                             int n_threads  = 0) {

    VARResult vr;
    vr.beta      = Rcpp::as<arma::mat>(var_result["beta"]);
    vr.residuals = Rcpp::as<arma::mat>(var_result["residuals"]);
    vr.sigma     = Rcpp::as<arma::mat>(var_result["sigma"]);
    vr.p         = Rcpp::as<int>(var_result["p"]);
    vr.c         = Rcpp::as<int>(var_result["c"]);
    vr.n_exog    = 0;

    BootHeteroResult res = fbootstrapHetero_cpp(y, vr, Z, indsR1, adjustu,
                                                 nboot, blocksize, hor,
                                                 nvar, scale, prc, prc2,
                                                 n_threads);

    return Rcpp::List::create(
        Rcpp::Named("upper")     = res.upper,
        Rcpp::Named("lower")     = res.lower,
        Rcpp::Named("upper2")    = res.upper2,
        Rcpp::Named("lower2")    = res.lower2,
        Rcpp::Named("meanirf")   = res.meanirf,
        Rcpp::Named("medianirf") = res.medianirf,
        Rcpp::Named("point")     = res.point
    );
}
