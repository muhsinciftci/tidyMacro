// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include "fbootstrapUhligCorrected.h"
#include "fbootstrapUhlig.h"
#include "fVAR.h"
#include "fbootstrapVAR.h"
#include "fcompanionMatrix.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstdio>
#include <algorithm>

static void bias_correct_uhlig(const arma::mat& beta, int c, int p,
                                const arma::cube& boot_betas,
                                arma::mat& Beta_t, int& corrections) {
    const int n_coef = static_cast<int>(beta.n_rows);
    const int N      = static_cast<int>(beta.n_cols);
    const int nboot1 = static_cast<int>(boot_betas.n_slices);

    arma::mat boot_mean(N, n_coef, arma::fill::zeros);
    for (int b = 0; b < nboot1; ++b) {
        boot_mean += boot_betas.slice(b);
    }
    boot_mean /= static_cast<double>(nboot1);

    const arma::mat beta_t = beta.t();
    arma::mat bias = boot_mean - beta_t;

    CompanionMatrixResult comp0 = fcompanionMatrix_cpp(beta, c, p);
    arma::cx_vec ev0 = arma::eig_gen(comp0.comp);
    double max_ev = arma::max(arma::abs(ev0));

    corrections = 1;

    if (max_ev >= 1.0) {
        Beta_t = beta_t;
        return;
    }

    Beta_t = beta_t - bias;

    {
        CompanionMatrixResult comp1 = fcompanionMatrix_cpp(Beta_t.t(), c, p);
        arma::cx_vec ev1 = arma::eig_gen(comp1.comp);
        max_ev = arma::max(arma::abs(ev1));
    }

    double delta = 1.0;
    while (max_ev >= 1.0) {
        delta -= 0.01;
        corrections += 1;
        if (delta < 0.0 || corrections > 200) {
            Beta_t = beta_t;
            break;
        }
        Beta_t = beta_t - bias * delta;
        CompanionMatrixResult comp_loop = fcompanionMatrix_cpp(Beta_t.t(), c, p);
        arma::cx_vec ev_loop = arma::eig_gen(comp_loop.comp);
        max_ev = arma::max(arma::abs(ev_loop));
    }
}

BootstrapUhligResult
fbootstrapUhligCorrected_cpp(const arma::mat& y, const VARResult& var_result,
                              int nboot1, int nboot2, int horizon, int idx,
                              double prc, const arma::uvec& cumulate,
                              int n_threads) {

    const int p      = var_result.p;
    const int c      = var_result.c;
    const int N      = static_cast<int>(y.n_cols);
    const int n_coef = static_cast<int>(var_result.beta.n_rows);

    int actual_threads = 1;
#ifdef _OPENMP
    actual_threads = (n_threads == 0)
                         ? std::max(1, omp_get_max_threads() - 1)
                         : n_threads;
    omp_set_num_threads(actual_threads);
    std::printf("[Pass 1] Bias estimation: %d reps, %d thread(s)\n",
                nboot1, actual_threads);
#else
    std::printf("[Pass 1] Bias estimation: %d reps, single-threaded\n", nboot1);
#endif

    arma::cube boot_beta_1(N, n_coef, nboot1, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int b = 0; b < nboot1; ++b) {
        BootstrapVARResult boot_data = fbootstrapVAR_cpp(y, var_result, "residual");
        VARResult var_loop = fVAR_cpp(boot_data.ynext, p, c, R_NilValue);
        boot_beta_1.slice(b) = var_loop.beta.t();
    }

    arma::mat Beta_t;
    int corrections = 1;
    bias_correct_uhlig(var_result.beta, c, p, boot_beta_1, Beta_t, corrections);

    if (corrections > 1)
        std::printf("[Bias correction] %d shrinkage iteration(s)\n", corrections);
    else
        std::printf("[Bias correction] Full correction applied\n");

#ifdef _OPENMP
    std::printf("[Pass 2] Bias-corrected bootstrap: %d reps, %d thread(s)\n",
                nboot2, actual_threads);
#else
    std::printf("[Pass 2] Bias-corrected bootstrap: %d reps, single-threaded\n", nboot2);
#endif

    VARResult corrected_var = var_result;
    corrected_var.beta      = Beta_t.t();

    return fbootstrapUhlig_cpp(y, corrected_var, nboot2, horizon, idx,
                               prc, cumulate, actual_threads);
}

//' @export
// [[Rcpp::export]]
Rcpp::List fbootstrapUhligCorrected(const arma::mat& y, const Rcpp::List& var_result,
                                     int nboot1, int nboot2, int horizon, int idx,
                                     double prc, const arma::uvec& cumulate,
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
    BootstrapUhligResult res = fbootstrapUhligCorrected_cpp(y, vr, nboot1, nboot2,
                                                            horizon, idx - 1,
                                                            prc, cumulate_cpp,
                                                            n_threads);

    Rcpp::NumericVector bootuhlig_out(res.bootuhlig_flat.begin(), res.bootuhlig_flat.end());
    bootuhlig_out.attr("dim") = Rcpp::IntegerVector::create(res.N, res.H, nboot2);

    return Rcpp::List::create(
        Rcpp::Named("bootuhlig") = bootuhlig_out,
        Rcpp::Named("upper")     = res.upper,
        Rcpp::Named("lower")     = res.lower,
        Rcpp::Named("boot_beta") = res.boot_beta
    );
}
