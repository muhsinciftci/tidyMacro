// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include "fbootstrapMaxCorrected.h"
#include "fbootstrapMax.h"
#include "fVAR.h"
#include "fbootstrapVAR.h"
#include "fcompanionMatrix.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstdio>
#include <algorithm>

// Takes the already-averaged Pass-1 coefficient mean (boot_mean, N x n_coef).
static void bias_correct_max(const arma::mat& beta, int c, int p,
                              const arma::mat& boot_mean,
                              arma::mat& Beta_t, int& corrections) {
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

BootstrapMaxResult
fbootstrapMaxCorrected_cpp(const arma::mat& y, const VARResult& var_result,
                           int nboot1, int nboot2, int horizon, int var_idx,
                           double prc, double prc2, const arma::uvec& cumulate,
                           Rcpp::Nullable<arma::vec> scaling,
                           Rcpp::Nullable<arma::mat> exog,
                           int n_threads) {

    const int p      = var_result.p;
    const int c      = var_result.c;
    const int n_exog = var_result.n_exog;
    const int N      = static_cast<int>(y.n_cols);
    const int n_coef = static_cast<int>(var_result.beta.n_rows);

    if (n_exog > 0 && exog.isNull())
        Rcpp::stop("Original VAR used exogenous variables. You must provide the 'exog' parameter.");
    if (n_exog == 0 && exog.isNotNull())
        Rcpp::stop("Original VAR did not use exogenous variables. Do not provide the 'exog' parameter.");

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

    const bool has_exog = exog.isNotNull();
    arma::mat exog_mat;
    if (has_exog) exog_mat = Rcpp::as<arma::mat>(exog);

    arma::mat boot_mean(N, n_coef, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp parallel
    {
        arma::mat thread_sum(N, n_coef, arma::fill::zeros);
#pragma omp for schedule(dynamic) nowait
        for (int b = 0; b < nboot1; ++b) {
            BootstrapVARResult boot_data = fbootstrapVAR_cpp(y, var_result, "residual");
            VARResult var_loop = has_exog
                                     ? fVAR_cpp(boot_data.ynext, p, c, exog_mat)
                                     : fVAR_cpp(boot_data.ynext, p, c, R_NilValue);
            thread_sum += var_loop.beta.t();
        }
#pragma omp critical
        boot_mean += thread_sum;
    }
#else
    for (int b = 0; b < nboot1; ++b) {
        BootstrapVARResult boot_data = fbootstrapVAR_cpp(y, var_result, "residual");
        VARResult var_loop = has_exog
                                 ? fVAR_cpp(boot_data.ynext, p, c, exog_mat)
                                 : fVAR_cpp(boot_data.ynext, p, c, R_NilValue);
        boot_mean += var_loop.beta.t();
    }
#endif
    boot_mean /= static_cast<double>(nboot1);

    arma::mat Beta_t;
    int corrections = 1;
    bias_correct_max(var_result.beta, c, p, boot_mean, Beta_t, corrections);

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

    return fbootstrapMax_cpp(y, corrected_var, nboot2, horizon, var_idx,
                             prc, prc2, cumulate, scaling, exog, actual_threads);
}

//' @export
// [[Rcpp::export]]
Rcpp::List fbootstrapMaxCorrected(const arma::mat& y, const Rcpp::List& var_result,
                                   int nboot1, int nboot2, int horizon, int var_idx,
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
    BootstrapMaxResult res = fbootstrapMaxCorrected_cpp(y, vr, nboot1, nboot2,
                                                        horizon, var_idx - 1,
                                                        prc, prc2, cumulate_cpp,
                                                        scaling, exog, n_threads);

    Rcpp::NumericVector bootmax_out(res.bootmax_flat.begin(), res.bootmax_flat.end());
    bootmax_out.attr("dim") = Rcpp::IntegerVector::create(res.N, res.H, nboot2);

    return Rcpp::List::create(
        Rcpp::Named("bootmax")   = bootmax_out,
        Rcpp::Named("upper")     = res.upper,
        Rcpp::Named("lower")     = res.lower,
        Rcpp::Named("upper2")    = res.upper2,
        Rcpp::Named("lower2")    = res.lower2,
        Rcpp::Named("boot_beta") = res.boot_beta
    );
}
