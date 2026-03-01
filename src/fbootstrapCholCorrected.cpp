// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#define ARMA_DONT_USE_OPENMP

#include "fbootstrapCholCorrected.h"
#include "fbootstrapChol.h"
#include "fbootstrapVAR.h"
#include "fcompanionMatrix.h"
#include "fcholeskyIRF.h"
#include "fwoldIRF.h"
#include "fVAR.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstdio>
#include <algorithm>

// Internal bias-correction helper (mirrors remove_bias.m)
static void bias_correct(const arma::mat&  beta,
                         int               c,
                         int               p,
                         const arma::cube& boot_betas,
                         arma::mat&        Beta_t,
                         int&              corrections)
{
    const int n_coef  = beta.n_rows;
    const int N       = beta.n_cols;
    const int nboot1  = boot_betas.n_slices;

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


BootstrapCholCorrectedResult
fbootstrapCholCorrected_cpp(const arma::mat&          y,
                            const VARResult&           var_result,
                            int                        nboot1,
                            int                        nboot2,
                            int                        horizon,
                            double                     prc,
                            const std::string&         bootscheme,
                            Rcpp::Nullable<arma::mat>  exog,
                            int                        n_threads = 0)
{
    const arma::mat& beta   = var_result.beta;
    const int        p      = var_result.p;
    const int        c      = var_result.c;
    const int        n_exog = var_result.n_exog;
    const int        N      = y.n_cols;
    const int        n_coef = beta.n_rows;

    if (n_exog > 0 && exog.isNull()) {
        Rcpp::stop("Original VAR used exogenous variables — provide 'exog'.");
    }
    if (n_exog == 0 && exog.isNotNull()) {
        Rcpp::stop("Original VAR did not use exogenous variables — omit 'exog'.");
    }

    int actual_threads = 1;
#ifdef _OPENMP
    actual_threads = (n_threads == 0)
                         ? std::max(1, omp_get_max_threads() - 1)
                         : n_threads;
    omp_set_num_threads(actual_threads);
    std::printf("[Pass 1] Bias estimation: %d boot reps, %d thread(s)\n",
                nboot1, actual_threads);
#else
    std::printf("[Pass 1] Bias estimation: %d boot reps, single-threaded\n", nboot1);
#endif

    // ================================================================== //
    //  PASS 1 — Bootstrap to estimate bias
    // ================================================================== //
    arma::cube boot_beta_1(N, n_coef, nboot1, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int b = 0; b < nboot1; ++b) {
        BootstrapVARResult boot_data = fbootstrapVAR_cpp(y, var_result, bootscheme);

        VARResult var_loop;
        if (exog.isNotNull()) {
            var_loop = fVAR_cpp(boot_data.ynext, p, c, exog);
        } else {
            var_loop = fVAR_cpp(boot_data.ynext, p, c, R_NilValue);
        }

        boot_beta_1.slice(b) = var_loop.beta.t();
    }

    // ================================================================== //
    //  BIAS CORRECTION
    // ================================================================== //
    arma::mat Beta_t;
    int       corrections = 1;

    bias_correct(beta, c, p, boot_beta_1, Beta_t, corrections);

    if (corrections > 1) {
        std::printf("[Bias correction] %d shrinkage iteration(s) applied\n", corrections);
    } else {
        std::printf("[Bias correction] Full correction applied (stable on first attempt)\n");
    }

    // ================================================================== //
    //  PASS 2 — Bias-corrected bootstrap for confidence bands
    // ================================================================== //
#ifdef _OPENMP
    std::printf("[Pass 2] Bias-corrected bootstrap: %d boot reps, %d thread(s)\n",
                nboot2, actual_threads);
#else
    std::printf("[Pass 2] Bias-corrected bootstrap: %d boot reps, single-threaded\n", nboot2);
#endif

    VARResult corrected_var  = var_result;
    corrected_var.beta       = Beta_t.t();

    BootstrapCholResult pass2 = fbootstrapChol_cpp(
        y, corrected_var, nboot2, horizon, prc, bootscheme, exog, actual_threads);

    BootstrapCholCorrectedResult result;
    result.bootchol_flat = pass2.bootchol_flat;
    result.upper         = pass2.upper;
    result.lower         = pass2.lower;
    result.boot_beta     = pass2.boot_beta;
    result.Beta          = Beta_t;
    result.corrections   = corrections;
    result.N             = pass2.N;
    result.H             = pass2.H;
    return result;
}


//' Bias-Corrected Bootstrap Cholesky IRF Confidence Bands (Kilian 1998)
//'
//' @param y T x N matrix of endogenous variables.
//' @param var_result List returned by \code{fVAR()}.
//' @param nboot1 First-pass bootstrap replications for bias estimation.
//' @param nboot2 Second-pass bootstrap replications for confidence bands.
//' @param horizon Maximum IRF horizon.
//' @param prc Confidence level in percent (e.g. 68 or 90).
//' @param bootscheme \code{"residual"} or \code{"wild"}.
//' @param exog Optional T x M matrix of exogenous variables. Default NULL.
//' @param n_threads OpenMP threads. 0 = all cores minus one.
//'
//' @return List: \code{bootchol}, \code{upper}, \code{lower},
//'   \code{boot_beta}, \code{Beta}, \code{corrections}.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fbootstrapCholCorrected(const arma::mat&          y,
                                   const Rcpp::List&          var_result,
                                   int                        nboot1,
                                   int                        nboot2,
                                   int                        horizon,
                                   double                     prc,
                                   const std::string&         bootscheme,
                                   Rcpp::Nullable<arma::mat>  exog      = R_NilValue,
                                   int                        n_threads = 0)
{
    VARResult vr;
    vr.beta       = Rcpp::as<arma::mat>(var_result["beta"]);
    vr.residuals  = Rcpp::as<arma::mat>(var_result["residuals"]);
    vr.sigma_full = Rcpp::as<arma::mat>(var_result["sigma_full"]);
    vr.p          = Rcpp::as<int>(var_result["p"]);
    vr.c          = Rcpp::as<int>(var_result["c"]);
    vr.n_exog     = var_result.containsElementNamed("n_exog")
                        ? Rcpp::as<int>(var_result["n_exog"]) : 0;

    BootstrapCholCorrectedResult result = fbootstrapCholCorrected_cpp(
        y, vr, nboot1, nboot2, horizon, prc, bootscheme, exog, n_threads);

    Rcpp::NumericVector bootchol_out(result.bootchol_flat.begin(),
                                     result.bootchol_flat.end());
    bootchol_out.attr("dim") = Rcpp::IntegerVector::create(
        result.N, result.N, result.H, nboot2);

    return Rcpp::List::create(
        Rcpp::Named("bootchol")    = bootchol_out,
        Rcpp::Named("upper")       = result.upper,
        Rcpp::Named("lower")       = result.lower,
        Rcpp::Named("boot_beta")   = result.boot_beta,
        Rcpp::Named("Beta")        = result.Beta,
        Rcpp::Named("corrections") = result.corrections
    );
}
