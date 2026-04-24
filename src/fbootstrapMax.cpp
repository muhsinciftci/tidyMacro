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
// Maximises wold.slice(H-1).row(var_idx) * S * h2  (last-horizon response only,
// matching MATLAB irfs(var,:,end)) subject to h2(0) = 0, ||h2|| = 1.
// Optimal: h2 = [0; M_sub / ||M_sub||] where
// M_sub = (wold.slice(H-1).row(var_idx) * S).cols(1, N-1).
static arma::vec lr_max_solve(const arma::cube& wold, const arma::mat& S,
                               int var_idx) {
    const int H = static_cast<int>(wold.n_slices);
    const int N = static_cast<int>(wold.n_rows);

    // Maximise the last-horizon response of var_idx (MATLAB LR_max.m uses irfs(var,:,end))
    arma::rowvec M_full = wold.slice(H - 1).row(var_idx) * S;

    arma::rowvec M_sub = M_full.cols(1, N - 1);
    double M_norm = arma::norm(M_sub);

    arma::vec h2(N, arma::fill::zeros);
    if (M_norm > 1e-14) {
        h2.rows(1, N - 1) = M_sub.t() / M_norm;
    }
    return h2;
}


BootstrapMaxResult
fbootstrapMax_cpp(const arma::mat& y, const VARResult& var_result,
                  int nboot, int horizon, int var_idx, double prc, double prc2,
                  const arma::uvec& cumulate,
                  Rcpp::Nullable<arma::vec> scaling,
                  Rcpp::Nullable<arma::mat> exog,
                  int n_threads) {

    const int p        = var_result.p;
    const int c        = var_result.c;
    const int n_exog   = var_result.n_exog;
    const int T        = static_cast<int>(y.n_rows);
    const int N        = static_cast<int>(y.n_cols);
    const int H        = horizon + 1;
    const int n_coef   = static_cast<int>(var_result.beta.n_rows);
    const int slice_sz = N * H;

    if (n_exog > 0 && exog.isNull())
        Rcpp::stop("Original VAR used exogenous variables. You must provide the 'exog' parameter.");
    if (n_exog == 0 && exog.isNotNull())
        Rcpp::stop("Original VAR did not use exogenous variables. Do not provide the 'exog' parameter.");

    const double df = static_cast<double>(T - 1 - p - N * p);
    if (df <= 0.0)
        Rcpp::stop("Degrees of freedom T-1-p-N*p = %.0f <= 0.", df);

    arma::mat bootmax_flat(slice_sz, nboot, arma::fill::zeros);
    arma::cube boot_beta(N, n_coef, nboot, arma::fill::zeros);

    const bool has_scaling = scaling.isNotNull();
    arma::vec scaling_vec;
    if (has_scaling) scaling_vec = Rcpp::as<arma::vec>(scaling);

    const bool has_exog = exog.isNotNull();
    arma::mat exog_mat;
    if (has_exog) exog_mat = Rcpp::as<arma::mat>(exog);

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

    for (int b = 0; b < nboot; ++b) {

        BootstrapVARResult boot_data = fbootstrapVAR_cpp(y, var_result, "residual",
                                                          has_exog ? &exog_mat : nullptr);

        VARResult var_loop = has_exog
                                 ? fVAR_cpp(boot_data.ynext, p, c, exog_mat)
                                 : fVAR_cpp(boot_data.ynext, p, c, R_NilValue);
        boot_beta.slice(b) = var_loop.beta.t();

        WoldIRFResult wold_res = fwoldIRF_cpp(var_loop, horizon);
        const arma::cube& wold_loop = wold_res.irfwold;

        arma::mat S_loop;
        if (!arma::chol(S_loop, var_loop.sigma, "lower")) {
            arma::mat reg = var_loop.sigma;
            reg.diag() += 1e-8 * arma::trace(reg) / N;
            arma::chol(S_loop, reg, "lower");
        }

        arma::vec h2 = lr_max_solve(wold_loop, S_loop, var_idx);

        // Precompute impact = S_loop * h2 once per draw — reused for sign
        // check and every horizon below.
        arma::vec impact = S_loop * h2;

        // Sign: ensure last-horizon non-structural response of var 0 is non-negative
        if (arma::dot(wold_loop.slice(H - 1).row(0), impact) < 0.0)
            impact = -impact;

        arma::mat struct_irf(N, H, arma::fill::none);
        for (int hh = 0; hh < H; ++hh) {
            struct_irf.col(hh) = wold_loop.slice(hh) * impact;
        }

        if (has_scaling) {
            int s_idx = static_cast<int>(scaling_vec(0)) - 1;
            double scale_val = struct_irf(s_idx, 0) * scaling_vec(1);
            if (std::abs(scale_val) > 1e-14) struct_irf /= scale_val;
        }

        for (arma::uword ci = 0; ci < cumulate.n_elem; ++ci) {
            arma::uword ri = cumulate(ci);
            struct_irf.row(ri) = arma::cumsum(struct_irf.row(ri));
        }

        bootmax_flat.col(b) = arma::vectorise(struct_irf);
    }

    const double up_pct   = 50.0 + prc  * 0.5;
    const double low_pct  = 50.0 - prc  * 0.5;
    const double up_pct2  = 50.0 + prc2 * 0.5;
    const double low_pct2 = 50.0 - prc2 * 0.5;

    arma::mat upper (N, H, arma::fill::zeros);
    arma::mat lower (N, H, arma::fill::zeros);
    arma::mat upper2(N, H, arma::fill::zeros);
    arma::mat lower2(N, H, arma::fill::zeros);

    auto nth_pct = [](std::vector<double>& v, double pct) -> double {
        const int n = static_cast<int>(v.size());
        const double raw = (pct / 100.0) * (n - 1);
        const int lo = static_cast<int>(raw);
        const double frac = raw - lo;
        std::nth_element(v.begin(), v.begin() + lo, v.end());
        const double lo_val = v[lo];
        if (frac < 1e-12 || lo + 1 >= n) return lo_val;
        return lo_val * (1.0 - frac) +
               *std::min_element(v.begin() + lo + 1, v.end()) * frac;
    };

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < slice_sz; ++i) {
        arma::rowvec row_rv = bootmax_flat.row(i);
        std::vector<double> sv(row_rv.begin(), row_rv.end());
        int row_i = i % N;
        int col_i = i / N;
        upper (row_i, col_i) = nth_pct(sv, up_pct);
        lower (row_i, col_i) = nth_pct(sv, low_pct);
        upper2(row_i, col_i) = nth_pct(sv, up_pct2);
        lower2(row_i, col_i) = nth_pct(sv, low_pct2);
    }

    BootstrapMaxResult result;
    result.bootmax_flat = bootmax_flat;
    result.upper        = upper;
    result.lower        = lower;
    result.upper2       = upper2;
    result.lower2       = lower2;
    result.boot_beta    = boot_beta;
    result.N            = N;
    result.H            = H;
    return result;
}

//' @export
// [[Rcpp::export]]
Rcpp::List fbootstrapMax(const arma::mat& y, const Rcpp::List& var_result,
                         int nboot, int horizon, int var_idx,
                         double prc = 90.0, double prc2 = 68.0,
                         Rcpp::IntegerVector cumulate = Rcpp::IntegerVector(),
                         Rcpp::Nullable<arma::vec> scaling = R_NilValue,
                         Rcpp::Nullable<arma::mat> exog    = R_NilValue,
                         int n_threads = 0) {
    VARResult vr;
    vr.beta       = Rcpp::as<arma::mat>(var_result["beta"]);
    vr.residuals  = Rcpp::as<arma::mat>(var_result["residuals"]);
    vr.sigma = Rcpp::as<arma::mat>(var_result["sigma"]);
    vr.p          = Rcpp::as<int>(var_result["p"]);
    vr.c          = Rcpp::as<int>(var_result["c"]);
    vr.n_exog     = var_result.containsElementNamed("n_exog")
                        ? Rcpp::as<int>(var_result["n_exog"]) : 0;

    arma::uvec cumulate_cpp = (cumulate.size() == 0)
        ? arma::uvec()
        : Rcpp::as<arma::uvec>(cumulate) - 1;
    BootstrapMaxResult res = fbootstrapMax_cpp(y, vr, nboot, horizon,
                                               var_idx - 1, prc, prc2,
                                               cumulate_cpp, scaling,
                                               exog, n_threads);

    Rcpp::NumericVector bootmax_out(res.bootmax_flat.begin(), res.bootmax_flat.end());
    bootmax_out.attr("dim") = Rcpp::IntegerVector::create(res.N, res.H, nboot);

    return Rcpp::List::create(
        Rcpp::Named("bootmax")   = bootmax_out,
        Rcpp::Named("upper")     = res.upper,
        Rcpp::Named("lower")     = res.lower,
        Rcpp::Named("upper2")    = res.upper2,
        Rcpp::Named("lower2")    = res.lower2,
        Rcpp::Named("boot_beta") = res.boot_beta
    );
}
