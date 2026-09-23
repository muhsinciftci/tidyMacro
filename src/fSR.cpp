// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// CRITICAL: Define this BEFORE including RcppArmadillo to prevent conflicts
#define ARMA_DONT_USE_OPENMP

// Bayesian sign- and narrative-restriction inference for SVARs.
// Port of SR.m from Cesa-Bianchi's VAR Toolbox 4.0, restructured so that the
// accepted draws are produced in parallel.
//
// The toolbox runs one sequential while-loop that keeps drawing until `ndraws`
// rotations have been accepted.  Here each accepted draw is a slot with its own
// deterministic RNG stream, so slots are independent, the loop parallelises
// cleanly, and the output does not depend on the number of threads.

#include "fVAR.h"
#include "fWoldIRF.h"
#include "fVARPosterior.h"
#include "fSignRestrictions.h"
#include "fCheckNarrative.h"
#include "rng_tidymacro.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <cstdio>
#include <vector>
#include <limits>

namespace {

// Linear-interpolated percentile, matching the convention already used by the
// bootstrap routines in this package.  `v` is reordered in place.
double nth_pct(std::vector<double>& v, double pct) {
    const int n = static_cast<int>(v.size());
    if (n == 1) return v[0];
    const double raw = (pct / 100.0) * (n - 1);
    const int    lo  = static_cast<int>(raw);
    const double frac = raw - lo;
    std::nth_element(v.begin(), v.begin() + lo, v.end());
    const double lo_val = v[lo];
    if (frac < 1e-12 || lo + 1 >= n) return lo_val;
    return lo_val * (1.0 - frac) +
           *std::min_element(v.begin() + lo + 1, v.end()) * frac;
}

// Weighted percentile used when ADRR importance weights are switched on.
// `idx` is a scratch permutation buffer reused across calls.
double nth_pct_w(const std::vector<double>& v, const std::vector<double>& w,
                 std::vector<int>& idx, double pct) {
    const int n = static_cast<int>(v.size());
    double total = 0.0;
    for (int i = 0; i < n; ++i) total += w[idx[i]];
    if (!(total > 0.0)) return v[idx[n / 2]];

    const double target = (pct / 100.0) * total;
    double cum = 0.0;
    for (int i = 0; i < n; ++i) {
        const double next = cum + w[idx[i]];
        if (next >= target) {
            if (i == 0 || w[idx[i]] <= 0.0) return v[idx[i]];
            // Interpolate inside the cell that straddles the target mass.
            const double frac = (target - cum) / w[idx[i]];
            const double prev = v[idx[i - 1]];
            return prev + frac * (v[idx[i]] - prev);
        }
        cum = next;
    }
    return v[idx[n - 1]];
}

// Structural IRFs and the implied FEVD for one impact matrix.
// irf(i, j, h) = (Psi_h B)(i, j); vd(i, j, h) is the share of the h-step
// forecast error variance of variable i explained by shock j.
void ir_and_vd(const arma::cube& wold, const arma::mat& B, int nsteps,
               arma::cube& irf, arma::cube& vd, arma::mat& cum_sq,
               bool want_vd) {

    const arma::uword N = B.n_rows;
    if (irf.n_slices != static_cast<arma::uword>(nsteps)) irf.set_size(N, N, nsteps);

    for (int h = 0; h < nsteps; ++h) irf.slice(h) = wold.slice(h) * B;

    if (!want_vd) return;
    if (vd.n_slices != static_cast<arma::uword>(nsteps)) vd.set_size(N, N, nsteps);
    cum_sq.zeros(N, N);

    for (int h = 0; h < nsteps; ++h) {
        cum_sq += arma::square(irf.slice(h));
        for (arma::uword i = 0; i < N; ++i) {
            double tot = 0.0;
            for (arma::uword j = 0; j < N; ++j) tot += cum_sq(i, j);
            if (tot > 1e-15) {
                for (arma::uword j = 0; j < N; ++j) vd(i, j, h) = cum_sq(i, j) / tot;
            } else {
                for (arma::uword j = 0; j < N; ++j) vd(i, j, h) = 1.0 / static_cast<double>(N);
            }
        }
    }
}

} // namespace

//' Bayesian Sign and Narrative Restrictions for SVARs
//'
//' Identifies a structural VAR by sign restrictions, optionally combined with
//' Antolin-Diaz and Rubio-Ramirez narrative restrictions and with an external
//' instrument that pins down the first impact column. Inference is Bayesian:
//' reduced-form parameters are drawn from their flat-prior Normal-inverse-Wishart
//' posterior and each draw is paired with a Haar-uniform rotation.
//'
//' @param y A T x N numeric matrix of endogenous variables.
//' @param p Integer lag order.
//' @param c Integer intercept indicator (1 = include, 0 = exclude).
//' @param SIGN An N x ds sign-restriction matrix: \code{+1} the response must be
//'   non-negative, \code{-1} non-positive, \code{0} unrestricted. \code{ds} must
//'   equal N minus the number of columns of \code{Bfix}.
//' @param nsteps Integer number of horizons to report, counting impact as the
//'   first (default 40).
//' @param ndraws Integer number of accepted draws to collect (default 500).
//' @param sr_hor Integer. Sign restrictions are imposed at horizons
//'   \code{0, ..., sr_hor - 1} (default 1, impact only).
//' @param sr_rot Integer maximum rotations attempted per parameter draw
//'   (default 500).
//' @param max_post_draws Integer maximum parameter draws attempted per accepted
//'   draw before that slot is abandoned (default 1000).
//' @param conf Numeric coverage of the reported credible bands, in percent
//'   (default 90).
//' @param inference Integer. \code{1} (default) draws reduced-form parameters
//'   from the posterior, so bands reflect both parameter and identification
//'   uncertainty; \code{0} holds them at OLS, leaving only the set of admissible
//'   rotations.
//' @param Bfix Optional N x q matrix of impact columns identified elsewhere and
//'   held fixed, which selects the \code{sign+iv} scheme (default NULL).
//' @param narr_sign_shock,narr_sign_period,narr_sign_sign Type-1 narrative
//'   restrictions: the structural shock \code{shock} at residual-sample row
//'   \code{period} must have the given \code{sign}. All 1-indexed; default NULL.
//' @param narr_dom_shock,narr_dom_period,narr_dom_var Type-2 narrative
//'   restrictions: at row \code{period}, shock \code{shock} must contribute more
//'   to the unexpected movement in variable \code{var} than all other shocks
//'   combined. All 1-indexed; default NULL.
//' @param narr_weight_mc Integer. When positive, each accepted draw is weighted
//'   by the ADRR importance weight, estimated with this many Monte Carlo
//'   replications, and bands become weighted percentiles. \code{0} (default)
//'   reproduces the plain rejection sampling of the VAR Toolbox.
//' @param resid_from_draw Logical. \code{FALSE} (default) evaluates narrative
//'   restrictions on the OLS residuals, as the VAR Toolbox does; \code{TRUE}
//'   recomputes residuals from each parameter draw.
//' @param store_draws Logical. \code{TRUE} (default) returns the full IRF and
//'   FEVD distributions across accepted draws; \code{FALSE} returns only medians
//'   and bands, which avoids copying two
//'   \code{(N * N * nsteps) x ndraws} matrices back into R.
//' @param exog Optional T x M matrix of exogenous regressors (default NULL).
//' @param n_threads Integer. \code{0} (default) uses all available OpenMP threads.
//' @param seed Integer base seed. Accepted draw \code{d} uses a stream derived
//'   from \code{seed} and \code{d}, so results are independent of the thread
//'   count.
//' @param verbose Logical; print thread count and acceptance diagnostics
//'   (default FALSE).
//'
//' @param bands_conf Optional vector of confidence levels, with the first equal
//'   to conf. All bands are computed in the same C++ summary pass.
//' @param fitted_var Optional precomputed fVAR result on exactly y, p, c and
//'   exog, supplied by the R wrapper to reuse its IV first-stage VAR fit.
//'
//' @return A list with medians and credible bands for the IRFs (\code{IRmed},
//'   \code{IRinf}, \code{IRsup}) and the FEVD (\code{VDmed}, \code{VDinf},
//'   \code{VDsup}), each an N x N x nsteps cube indexed
//'   \code{[variable, shock, horizon]}; the accepted impact matrices
//'   \code{Ball} and coefficient draws \code{beta_all}; the median impact matrix
//'   \code{Bmed} and the Fry-Pagan draw \code{Bfp} with its \code{IRfp},
//'   \code{VDfp} and index \code{fp_index}; the importance \code{weights};
//'   \code{accept_rate}, \code{ndraws_tried}, \code{n_tried} and
//'   \code{n_failed}; and, when \code{store_draws = TRUE}, the flattened draw
//'   distributions \code{IRall} and \code{VDall}, each
//'   \code{(N * N * nsteps) x ndraws} and reshapeable to
//'   \code{c(N, N, nsteps, ndraws)}.
//'
//' @details
//' Each accepted draw is produced independently: parameters are drawn from the
//' posterior, a Haar rotation is sought that satisfies \code{SIGN}, and the
//' candidate is then screened against the narrative restrictions. A slot that
//' exhausts \code{max_post_draws} parameter draws is reported in
//' \code{n_failed} rather than silently dropped.
//'
//' FEVD shares are returned on the \code{[0, 1]} scale, not in percent.
//'
//' The instrument-identified column is fixed at
//' its OLS point estimate while the rest of the system is redrawn, so the
//' reported bands omit the instrument's own sampling uncertainty. This
//' reproduces the VAR Toolbox exactly. The fixed column is rescaled by
//' \eqn{1/\|L^{-1}b_1\|} against each draw's Cholesky factor \eqn{L}, so its
//' direction is held constant but its length is not; this too matches the
//' toolbox.
//'
//' @references
//' Uhlig, H. (2005). What are the effects of monetary policy on output?
//' \emph{Journal of Monetary Economics}, 52(2), 381--419.
//'
//' Rubio-Ramirez, J. F., Waggoner, D. F., & Zha, T. (2010). Structural vector
//' autoregressions. \emph{Review of Economic Studies}, 77(2), 665--696.
//'
//' Antolin-Diaz, J., & Rubio-Ramirez, J. F. (2018). Narrative sign restrictions
//' for SVARs. \emph{American Economic Review}, 108(10), 2802--2829.
//'
//' @seealso \code{\link{fSignRestrictions_cpp}}, \code{\link{fCheckNarrative_cpp}},
//'   \code{\link{fRecoverBIV_cpp}}, \code{\link{fVARPosterior_cpp}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fSR_cpp(const arma::mat& y, int p, int c, const arma::mat& SIGN,
                   int nsteps = 40,
                   int ndraws = 500,
                   int sr_hor = 1,
                   int sr_rot = 500,
                   int max_post_draws = 1000,
                   double conf = 90.0,
                   int inference = 1,
                   Rcpp::Nullable<arma::mat> Bfix = R_NilValue,
                   Rcpp::Nullable<arma::uvec> narr_sign_shock = R_NilValue,
                   Rcpp::Nullable<arma::uvec> narr_sign_period = R_NilValue,
                   Rcpp::Nullable<arma::vec>  narr_sign_sign = R_NilValue,
                   Rcpp::Nullable<arma::uvec> narr_dom_shock = R_NilValue,
                   Rcpp::Nullable<arma::uvec> narr_dom_period = R_NilValue,
                   Rcpp::Nullable<arma::uvec> narr_dom_var = R_NilValue,
                   int narr_weight_mc = 0,
                   bool resid_from_draw = false,
                   bool store_draws = true,
                   Rcpp::Nullable<arma::mat> exog = R_NilValue,
                   int n_threads = 0,
                   int seed = 42,
                   bool verbose = false,
                   Rcpp::Nullable<arma::vec> bands_conf = R_NilValue,
                   Rcpp::Nullable<Rcpp::List> fitted_var = R_NilValue) {

    // ---- 0. validation ------------------------------------------------
    if (nsteps < 1)          Rcpp::stop("nsteps must be at least 1.");
    if (ndraws < 1)          Rcpp::stop("ndraws must be at least 1.");
    if (sr_hor < 1)          Rcpp::stop("sr_hor must be at least 1.");
    if (sr_rot < 1)          Rcpp::stop("sr_rot must be at least 1.");
    if (max_post_draws < 1)  Rcpp::stop("max_post_draws must be at least 1.");
    if (conf <= 0.0 || conf >= 100.0) Rcpp::stop("conf must lie strictly between 0 and 100.");
    if (inference != 0 && inference != 1) Rcpp::stop("inference must be 0 or 1.");
    if (p < 1 || (c != 0 && c != 1) || y.n_cols < 1 ||
        y.n_rows <= static_cast<arma::uword>(p))
        Rcpp::stop("Invalid VAR dimensions, lag order or intercept.");
    if (!std::isfinite(conf)) Rcpp::stop("conf must be finite.");
    if (!SIGN.is_finite() || arma::any(arma::vectorise((SIGN != -1.0) % (SIGN != 0.0) % (SIGN != 1.0))))
        Rcpp::stop("SIGN must contain only -1, 0 or 1.");
    arma::vec levels = bands_conf.isNotNull() ? Rcpp::as<arma::vec>(bands_conf)
                                             : arma::vec{conf};
    if (levels.is_empty() || !levels.is_finite() || arma::any(levels <= 0.0) ||
        arma::any(levels >= 100.0)) Rcpp::stop("Invalid band coverage levels.");
    if (levels(0) != conf) Rcpp::stop("The first bands_conf level must equal conf.");
    if (!y.is_finite())      Rcpp::stop("y contains non-finite values.");

    const arma::uword N = y.n_cols;
    if (SIGN.n_rows != N) Rcpp::stop("SIGN must have one row per variable.");

    arma::mat exog_mat;
    if (exog.isNotNull()) exog_mat = Rcpp::as<arma::mat>(exog);

    arma::mat Bfix_mat;
    if (Bfix.isNotNull()) Bfix_mat = Rcpp::as<arma::mat>(Bfix);
    const int n_fixed = static_cast<int>(Bfix_mat.n_cols);
    if (n_fixed > 0 && Bfix_mat.n_rows != N) {
        Rcpp::stop("Bfix must have one row per variable.");
    }
    if (static_cast<int>(SIGN.n_cols) + n_fixed != static_cast<int>(N)) {
        Rcpp::stop("SIGN must have N - ncol(Bfix) = %d columns, but has %d.",
                   static_cast<int>(N) - n_fixed, static_cast<int>(SIGN.n_cols));
    }

    // ---- 1. OLS VAR and posterior factorisation -----------------------
    if (exog.isNotNull() && (exog_mat.n_rows != y.n_rows || !exog_mat.is_finite()))
        Rcpp::stop("exog must be finite and have nrow(y) rows.");
    if (y.n_rows - p <= N * p + c + exog_mat.n_cols)
        Rcpp::stop("Too few observations for the VAR coefficients.");
    VARResult var_ols;
    if (fitted_var.isNotNull()) {
        const Rcpp::List fit(fitted_var.get());
        var_ols.beta = Rcpp::as<arma::mat>(fit["beta"]);
        var_ols.residuals = Rcpp::as<arma::mat>(fit["residuals"]);
        var_ols.sigma = Rcpp::as<arma::mat>(fit["sigma"]);
        var_ols.p = p; var_ols.c = c; var_ols.n_exog = exog_mat.n_cols;
        if (var_ols.beta.n_rows != N*p+c+exog_mat.n_cols || var_ols.beta.n_cols != N ||
            var_ols.residuals.n_rows != y.n_rows-p || var_ols.residuals.n_cols != N ||
            var_ols.sigma.n_rows != N || var_ols.sigma.n_cols != N ||
            !var_ols.beta.is_finite() || !var_ols.residuals.is_finite() || !var_ols.sigma.is_finite())
            Rcpp::stop("fitted_var has incompatible dimensions or non-finite values.");
    } else {
        var_ols = exog.isNotNull() ? fVAR_cpp(y, p, c, exog_mat)
                                  : fVAR_cpp(y, p, c, R_NilValue);
    }
    arma::mat Y, X;
    fVARDesign_cpp(y, p, c, exog_mat, Y, X);

    const int nobs      = static_cast<int>(X.n_rows);
    const int k         = static_cast<int>(X.n_cols);
    const int H_wold    = std::max(nsteps, sr_hor);

    NIWPosterior post;
    if (inference == 1) {
        post = fNIWPosteriorPrep_cpp(var_ols.beta, var_ols.sigma, X.t() * X, nobs);
    }

    // ---- 2. narrative restrictions ------------------------------------
    auto uvec_or_empty = [](Rcpp::Nullable<arma::uvec> v) {
        return v.isNotNull() ? Rcpp::as<arma::uvec>(v) : arma::uvec();
    };
    auto vec_or_empty = [](Rcpp::Nullable<arma::vec> v) {
        return v.isNotNull() ? Rcpp::as<arma::vec>(v) : arma::vec();
    };
    NarrativeRestrictions narr = fNarrativePrep_cpp(
        uvec_or_empty(narr_sign_shock), uvec_or_empty(narr_sign_period),
        vec_or_empty(narr_sign_sign),
        uvec_or_empty(narr_dom_shock), uvec_or_empty(narr_dom_period),
        uvec_or_empty(narr_dom_var), N, static_cast<arma::uword>(nobs));
    const bool use_weights = narr.active && narr_weight_mc > 0;

    // Wold multipliers are draw-invariant when the coefficients are held at OLS.
    arma::cube wold_fixed;
    if (inference == 0) wold_fixed = fWoldIRF_cpp(var_ols, H_wold - 1).irfwold;

    // ---- 4. storage ---------------------------------------------------
    const double size_needed = static_cast<double>(N) * N * nsteps;
    if (size_needed > std::numeric_limits<int>::max())
        Rcpp::stop("Requested horizon exceeds the matrix dimension limit.");
    const int slice_sz = static_cast<int>(size_needed);
    // R owns these buffers: workers write through Armadillo views, and returning
    // stored draws does not create a second copy of the two largest arrays.
    Rcpp::NumericMatrix IR_storage(Rcpp::no_init(slice_sz, ndraws));
    Rcpp::NumericMatrix VD_storage(Rcpp::no_init(slice_sz, ndraws));
    arma::mat IRall(IR_storage.begin(), slice_sz, ndraws, false, true);
    arma::mat VDall(VD_storage.begin(), slice_sz, ndraws, false, true);
    arma::cube Ball(N, N, ndraws, arma::fill::zeros);
    arma::cube beta_all(k, N, ndraws, arma::fill::zeros);
    arma::vec  weights(ndraws, arma::fill::ones);
    arma::ivec n_tried_v(ndraws, arma::fill::zeros);
    std::vector<char> slot_ok(ndraws, 0);
    std::vector<std::string> slot_error(ndraws);

    int actual_threads = 1;
#ifdef _OPENMP
    actual_threads = (n_threads <= 0) ? omp_get_max_threads()
                                      : n_threads;
    actual_threads = std::max(1, std::min(actual_threads, ndraws));
    if (verbose) std::printf("Using %d thread(s) for sign-restriction draws...\n",
                             actual_threads);
#else
    if (verbose) std::printf("OpenMP not available. Running single-threaded.\n");
#endif

    double total_rot = 0.0;

    // ---- 5. draw loop -------------------------------------------------
#ifdef _OPENMP
#pragma omp parallel num_threads(actual_threads) reduction(+ : total_rot)
#endif
    {
        NIWScratch     niw_scr;
        SignRotScratch rot_scr;
        NarrativeScratch narr_scr;
        arma::mat G, sigma_d, beta_d, B, resid_draw, cum_sq, E_mc;
        arma::cube wold_local, irf, vd;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
        for (int d = 0; d < ndraws; ++d) {
        // Each slot owns its stream; scratch is reused by the worker.
        tidymacro::RNG rng(static_cast<std::uint64_t>(seed) * 1000003ULL +
                           static_cast<std::uint64_t>(d) * 7919ULL + 1ULL);
        beta_d = var_ols.beta;
        sigma_d = var_ols.sigma;
        bool accepted = false;
        try {
        for (int t = 0; t < max_post_draws && !accepted; ++t) {

            const arma::cube* woldp;
            if (inference == 1) {
                fNIWPosteriorDraw_cpp(post, rng, niw_scr, G, sigma_d, beta_d);
                VARResult vd_draw;
                vd_draw.beta = beta_d;
                vd_draw.p = p; vd_draw.c = c; vd_draw.n_exog = var_ols.n_exog;
                wold_local = fWoldIRF_cpp(vd_draw, sr_hor - 1).irfwold;
                woldp = &wold_local;
            } else {
                woldp = &wold_fixed;
            }

            fSignRotationPrep_cpp(sigma_d, Bfix_mat, rng, rot_scr);

            int n_rot = 0;
            const bool found = fSignRotation_cpp(SIGN, *woldp, sr_hor, sr_rot,
                                                 n_fixed, rng, rot_scr, B, n_rot);
            total_rot += static_cast<double>(n_rot);
            n_tried_v(d) += n_rot;
            if (!found) continue;

            // The narrative screen follows `resid_from_draw`: the toolbox
            // screens on OLS residuals even when the coefficients have been
            // redrawn.
            if (narr.active && resid_from_draw && inference == 1) {
                resid_draw = Y - X * beta_d;
            }
            const arma::mat& resid_narr =
                (resid_from_draw && inference == 1) ? resid_draw : var_ols.residuals;

            if (narr.active && !fNarrativeCheck_cpp(B, resid_narr, narr, narr_scr)) continue;

            // Only accepted draws need the full reporting horizon.
            if (inference == 1 && nsteps > sr_hor) {
                VARResult vd_draw;
                vd_draw.beta = beta_d; vd_draw.p = p; vd_draw.c = c;
                vd_draw.n_exog = var_ols.n_exog;
                wold_local = fWoldIRF_cpp(vd_draw, nsteps - 1).irfwold;
            }
            // ---- accepted ---------------------------------------------
            ir_and_vd(*woldp, B, nsteps, irf, vd, cum_sq, true);
            IRall.col(d)      = arma::vectorise(irf);
            VDall.col(d)      = arma::vectorise(vd);
            Ball.slice(d)     = B;
            beta_all.slice(d) = beta_d;
            if (use_weights) {
                weights(d) = fNarrativeWeight_cpp(B, narr, narr_weight_mc, rng, E_mc);
            }
            slot_ok[d] = 1;
            accepted   = true;
        }
        } catch (const std::exception& e) {
            slot_error[d] = e.what();
            slot_ok[d] = 0;
        } catch (...) {
            slot_error[d] = "Unknown numerical error";
            // Numerically degenerate slot: leave it unaccepted and report it
            // through n_failed.  Throwing out of an OpenMP region is not safe.
            slot_ok[d] = 0;
        }
    }

    } // worker scratch

    // ---- 6. drop failed slots -----------------------------------------
    std::vector<arma::uword> keep_v;
    keep_v.reserve(ndraws);
    for (int d = 0; d < ndraws; ++d) if (slot_ok[d]) keep_v.push_back(static_cast<arma::uword>(d));

    const int n_ok = static_cast<int>(keep_v.size());
    int n_errors = 0;
    std::string first_error;
    for (const auto& error : slot_error) if (!error.empty()) {
        if (first_error.empty()) first_error = error;
        ++n_errors;
    }
    if (n_ok == 0 && n_errors > 0)
        Rcpp::stop("No accepted draw; %d slots failed numerically. First error: %s",
                   n_errors, first_error.c_str());
    if (n_ok == 0) {
        Rcpp::stop("No draw satisfied the restrictions. Loosen them, or raise "
                   "sr_rot / max_post_draws.");
    }
    const int n_failed = ndraws - n_ok;

    if (n_failed > 0) {
        const arma::uvec keep(keep_v);
        // Pack successful columns in place; their order is unchanged.
        for (int i = 0; i < n_ok; ++i) {
            if (keep_v[i] != static_cast<arma::uword>(i)) {
                std::copy_n(IRall.colptr(keep_v[i]), slice_sz, IRall.colptr(i));
                std::copy_n(VDall.colptr(keep_v[i]), slice_sz, VDall.colptr(i));
            }
        }

        arma::cube Ball_k(N, N, n_ok, arma::fill::none);
        arma::cube beta_k(k, N, n_ok, arma::fill::none);
        arma::vec  w_k(n_ok);
        arma::ivec t_k(n_ok);
        for (int i = 0; i < n_ok; ++i) {
            Ball_k.slice(i) = Ball.slice(keep_v[i]);
            beta_k.slice(i) = beta_all.slice(keep_v[i]);
            w_k(i)          = weights(keep_v[i]);
            t_k(i)          = n_tried_v(keep_v[i]);
        }
        Ball      = std::move(Ball_k);
        beta_all  = std::move(beta_k);
        weights   = std::move(w_k);
        n_tried_v = std::move(t_k);
    }

    // ---- 7. medians, bands, Fry-Pagan ---------------------------------
    arma::cube IRmed(N, N, nsteps), VDmed(N, N, nsteps);
    std::vector<arma::cube> IRlo(levels.n_elem), IRhi(levels.n_elem),
                            VDlo(levels.n_elem), VDhi(levels.n_elem);
    for (arma::uword l = 0; l < levels.n_elem; ++l) {
        IRlo[l].set_size(N, N, nsteps); IRhi[l].set_size(N, N, nsteps);
        VDlo[l].set_size(N, N, nsteps); VDhi[l].set_size(N, N, nsteps);
    }
    std::vector<double> wv(weights.begin(), weights.end());
#ifdef _OPENMP
#pragma omp parallel num_threads(actual_threads)
#endif
    {
        std::vector<double> values(n_ok);
        std::vector<int> idx(n_ok);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (int i = 0; i < slice_sz; ++i) {
            for (int which = 0; which < 2; ++which) {
                const arma::mat& draws = which == 0 ? IRall : VDall;
                for (int d = 0; d < n_ok; ++d) values[d] = draws(i, d);
                if (use_weights) {
                    for (int d = 0; d < n_ok; ++d) idx[d] = d;
                    std::sort(idx.begin(), idx.end(),
                              [&](int a, int b) { return values[a] < values[b]; });
                }
                auto percentile = [&](double pct) {
                    return use_weights ? nth_pct_w(values, wv, idx, pct)
                                       : nth_pct(values, pct);
                };
                (which == 0 ? IRmed : VDmed).at(i) = percentile(50.0);
                for (arma::uword l = 0; l < levels.n_elem; ++l) {
                    const double lo = (100.0 - levels(l)) / 2.0;
                    (which == 0 ? IRlo[l] : VDlo[l]).at(i) = percentile(lo);
                    (which == 0 ? IRhi[l] : VDhi[l]).at(i) = percentile(100.0 - lo);
                }
            }
        }
    }

    // Element-wise median of the accepted impact matrices (arma has no
    // cube-wise median reduction).
    arma::mat Bmed(N, N, arma::fill::none);
    {
        std::vector<double> col(Ball.n_slices);
        for (arma::uword j = 0; j < N; ++j) {
            for (arma::uword i = 0; i < N; ++i) {
                for (arma::uword d = 0; d < Ball.n_slices; ++d) col[d] = Ball(i, j, d);
                Bmed(i, j) = nth_pct(col, 50.0);
            }
        }
    }
    double best = arma::datum::inf;
    int    sel  = 0;
    for (int d = 0; d < static_cast<int>(Ball.n_slices); ++d) {
        const double dist = arma::accu(arma::square(Ball.slice(d) - Bmed));
        if (dist < best) { best = dist; sel = d; }
    }

    arma::cube IRfp(N, N, nsteps), VDfp(N, N, nsteps);
    std::copy(IRall.colptr(sel), IRall.colptr(sel) + slice_sz, IRfp.memptr());
    std::copy(VDall.colptr(sel), VDall.colptr(sel) + slice_sz, VDfp.memptr());

    const double accept_rate = (total_rot > 0.0)
                                 ? static_cast<double>(n_ok) / total_rot : 0.0;
    if (verbose) {
        std::printf("Accepted %d of %d slots; %.0f rotations tried (acceptance %.3f%%).\n",
                    n_ok, ndraws, total_rot, 100.0 * accept_rate);
    }

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("IRmed") = IRmed, Rcpp::Named("IRinf") = IRlo[0],
        Rcpp::Named("IRsup") = IRhi[0], Rcpp::Named("VDmed") = VDmed,
        Rcpp::Named("VDinf") = VDlo[0], Rcpp::Named("VDsup") = VDhi[0],
        Rcpp::Named("Ball") = Ball, Rcpp::Named("beta_all") = beta_all,
        Rcpp::Named("Bmed") = Bmed, Rcpp::Named("Bfp") = Ball.slice(sel),
        Rcpp::Named("IRfp") = IRfp, Rcpp::Named("VDfp") = VDfp,
        Rcpp::Named("fp_index") = sel + 1,
        Rcpp::Named("weights") = weights,
        Rcpp::Named("n_tried") = n_tried_v,
        Rcpp::Named("accept_rate") = accept_rate,
        Rcpp::Named("ndraws_tried") = total_rot,
        Rcpp::Named("n_failed") = n_failed,
        Rcpp::Named("n_errors") = n_errors,
        Rcpp::Named("first_error") = first_error);

    Rcpp::List bands(levels.n_elem);
    for (arma::uword l = 0; l < levels.n_elem; ++l) {
        bands[l] = Rcpp::List::create(Rcpp::Named("IRinf") = IRlo[l],
            Rcpp::Named("IRsup") = IRhi[l], Rcpp::Named("VDinf") = VDlo[l],
            Rcpp::Named("VDsup") = VDhi[l]);
    }
    out["bands"] = bands;
    out["var"] = Rcpp::List::create(Rcpp::Named("beta") = var_ols.beta,
        Rcpp::Named("residuals") = var_ols.residuals, Rcpp::Named("sigma") = var_ols.sigma,
        Rcpp::Named("p") = p, Rcpp::Named("c") = c, Rcpp::Named("n_exog") = var_ols.n_exog);
    if (store_draws) {
        if (n_failed == 0) {
            out["IRall"] = IR_storage;
            out["VDall"] = VD_storage;
        } else {
            out["IRall"] = arma::mat(IRall.cols(0, n_ok - 1));
            out["VDall"] = arma::mat(VDall.cols(0, n_ok - 1));
        }
    }
    return out;
}
