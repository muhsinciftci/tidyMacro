// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include "fbootstrapMax.h"
#include "fVAR.h"
#include "fbootstrapVAR.h"
#include "fwoldIRF.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstdio>
#include <algorithm>

// Closed-form LR_max solution (replaces MATLAB fminsearch on LR_max.m).
// Maximises sum_h wold.slice(h).row(var_idx) * S * h2
// subject to h2(0) = 0, ||h2|| = 1.
// Optimal: h2 = [0; M_sub / ||M_sub||] where
// M_sub = (sum_h wold.slice(h).row(var_idx) * S).cols(1, N-1).
static arma::vec lr_max_solve(const arma::cube& wold, const arma::mat& S,
                               int var_idx) {
    const int H = static_cast<int>(wold.n_slices);
    const int N = static_cast<int>(wold.n_rows);

    arma::rowvec M_full(N, arma::fill::zeros);
    for (int h = 0; h < H; ++h) {
        M_full += wold.slice(h).row(var_idx) * S;
    }

    arma::rowvec M_sub = M_full.cols(1, N - 1);
    double M_norm = arma::norm(M_sub);

    arma::vec h2(N, arma::fill::zeros);
    if (M_norm > 1e-14) {
        h2.rows(1, N - 1) = M_sub.t() / M_norm;
    }
    return h2;
}

static double pctile_interp(const arma::vec& sv, double pct) {
    const int n = static_cast<int>(sv.n_elem);
    double idx_d = (pct / 100.0) * (n - 1);
    int lo = static_cast<int>(idx_d);
    double frac = idx_d - lo;
    return (lo < n - 1) ? sv(lo) * (1.0 - frac) + sv(lo + 1) * frac : sv(n - 1);
}

BootstrapMaxResult
fbootstrapMax_cpp(const arma::mat& y, const VARResult& var_result,
                  int nboot, int horizon, int var_idx, double prc,
                  const arma::uvec& cumulate,
                  Rcpp::Nullable<arma::vec> scaling,
                  int n_threads) {

    const int p        = var_result.p;
    const int c        = var_result.c;
    const int T        = static_cast<int>(y.n_rows);
    const int N        = static_cast<int>(y.n_cols);
    const int H        = horizon + 1;
    const int n_coef   = static_cast<int>(var_result.beta.n_rows);
    const int slice_sz = N * H;

    const double df = static_cast<double>(T - 1 - p - N * p);
    if (df <= 0.0)
        Rcpp::stop("Degrees of freedom T-1-p-N*p = %.0f <= 0.", df);

    arma::mat bootmax_flat(slice_sz, nboot, arma::fill::zeros);
    arma::cube boot_beta(N, n_coef, nboot, arma::fill::zeros);

    int actual_threads = 1;
#ifdef _OPENMP
    actual_threads = (n_threads == 0)
                         ? std::max(1, omp_get_max_threads() - 1)
                         : n_threads;
    omp_set_num_threads(actual_threads);
    std::printf("Using %d thread(s) for LR-Max bootstrap...\n", actual_threads);
#else
    std::printf("OpenMP not available. Running single-threaded LR-Max bootstrap.\n");
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int b = 0; b < nboot; ++b) {

        BootstrapVARResult boot_data = fbootstrapVAR_cpp(y, var_result, "residual");
        VARResult var_loop = fVAR_cpp(boot_data.ynext, p, c, R_NilValue);
        boot_beta.slice(b) = var_loop.beta.t();

        WoldIRFResult wold_res = fwoldIRF_cpp(var_loop, horizon);
        const arma::cube& wold_loop = wold_res.irfwold;

        arma::mat S_loop;
        if (!arma::chol(S_loop, var_loop.sigma_full, "lower")) {
            arma::mat reg = var_loop.sigma_full;
            reg.diag() += 1e-8 * arma::trace(reg) / N;
            arma::chol(S_loop, reg, "lower");
        }

        arma::vec h2 = lr_max_solve(wold_loop, S_loop, var_idx);

        // Sign: ensure last-horizon non-structural response of var 0 is non-negative
        arma::mat last_ns = wold_loop.slice(H - 1) * S_loop;
        arma::rowvec last_ns_row = last_ns.row(0);
        if (arma::dot(last_ns_row, h2) < 0.0) h2 = -h2;

        arma::mat struct_irf(N, H, arma::fill::zeros);
        for (int hh = 0; hh < H; ++hh) {
            struct_irf.col(hh) = wold_loop.slice(hh) * S_loop * h2;
        }

        if (scaling.isNotNull()) {
            arma::vec s = Rcpp::as<arma::vec>(scaling);
            int s_idx = static_cast<int>(s(0)) - 1;
            double scale_val = struct_irf(s_idx, 0) / s(1);
            if (std::abs(scale_val) > 1e-14) struct_irf /= scale_val;
        }

        for (arma::uword ci = 0; ci < cumulate.n_elem; ++ci) {
            arma::uword ri = cumulate(ci);
            struct_irf.row(ri) = arma::cumsum(struct_irf.row(ri));
        }

        bootmax_flat.col(b) = arma::vectorise(struct_irf);
    }

    const double up_pct  = 50.0 + prc * 0.5;
    const double low_pct = 50.0 - prc * 0.5;

    arma::mat upper(N, H, arma::fill::zeros);
    arma::mat lower(N, H, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < slice_sz; ++i) {
        arma::vec sv = arma::sort(bootmax_flat.row(i).t());
        int row_i = i % N;
        int col_i = i / N;
        upper(row_i, col_i) = pctile_interp(sv, up_pct);
        lower(row_i, col_i) = pctile_interp(sv, low_pct);
    }

    BootstrapMaxResult result;
    result.bootmax_flat = bootmax_flat;
    result.upper        = upper;
    result.lower        = lower;
    result.boot_beta    = boot_beta;
    result.N            = N;
    result.H            = H;
    return result;
}

//' @export
// [[Rcpp::export]]
Rcpp::List fbootstrapMax(const arma::mat& y, const Rcpp::List& var_result,
                         int nboot, int horizon, int var_idx, double prc,
                         const arma::uvec& cumulate,
                         Rcpp::Nullable<arma::vec> scaling = R_NilValue,
                         int n_threads = 0) {
    VARResult vr;
    vr.beta       = Rcpp::as<arma::mat>(var_result["beta"]);
    vr.residuals  = Rcpp::as<arma::mat>(var_result["residuals"]);
    vr.sigma_full = Rcpp::as<arma::mat>(var_result["sigma_full"]);
    vr.p          = Rcpp::as<int>(var_result["p"]);
    vr.c          = Rcpp::as<int>(var_result["c"]);
    vr.n_exog     = var_result.containsElementNamed("n_exog")
                        ? Rcpp::as<int>(var_result["n_exog"]) : 0;

    arma::uvec cumulate_cpp = cumulate - 1;
    BootstrapMaxResult res = fbootstrapMax_cpp(y, vr, nboot, horizon,
                                               var_idx - 1, prc,
                                               cumulate_cpp, scaling, n_threads);

    Rcpp::NumericVector bootmax_out(res.bootmax_flat.begin(), res.bootmax_flat.end());
    bootmax_out.attr("dim") = Rcpp::IntegerVector::create(res.N, res.H, nboot);

    return Rcpp::List::create(
        Rcpp::Named("bootmax")   = bootmax_out,
        Rcpp::Named("upper")     = res.upper,
        Rcpp::Named("lower")     = res.lower,
        Rcpp::Named("boot_beta") = res.boot_beta
    );
}
