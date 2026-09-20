// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_DONT_USE_OPENMP

#include "fCheckNarrative.h"
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>

namespace {

// Position of `p` in the ascending unique-period list.
arma::uword slot_of(const arma::uvec& periods, arma::uword p) {
    const arma::uword* first = periods.memptr();
    const arma::uword* last  = first + periods.n_elem;
    return static_cast<arma::uword>(std::lower_bound(first, last, p) - first);
}

// One restriction set evaluated against a matrix of dated structural shocks
// (column s = shocks at the s-th unique period).
inline bool restrictions_hold(const arma::mat& B,
                              const NarrativeRestrictions& R,
                              const arma::mat& E) {

    for (arma::uword k = 0; k < R.ns_shock.n_elem; ++k) {
        if (R.ns_sign(k) * E(R.ns_shock(k), R.ns_slot(k)) <= 0.0) return false;
    }

    const arma::uword N = B.n_rows;
    for (arma::uword k = 0; k < R.nd_shock.n_elem; ++k) {
        const arma::uword i = R.nd_var(k);
        const arma::uword j = R.nd_shock(k);
        const double*     e = E.colptr(R.nd_slot(k));

        double own = 0.0, others = 0.0;
        for (arma::uword s = 0; s < N; ++s) {
            const double contrib = std::fabs(B(i, s) * e[s]);
            if (s == j) own = contrib; else others += contrib;
        }
        if (own <= others) return false;
    }
    return true;
}

} // namespace

// Nullable -> empty-vector helpers used by the R-facing wrapper.
arma::uvec as_uvec_or_empty(Rcpp::Nullable<arma::uvec> x) {
    return x.isNotNull() ? Rcpp::as<arma::uvec>(x) : arma::uvec();
}
arma::vec as_vec_or_empty(Rcpp::Nullable<arma::vec> x) {
    return x.isNotNull() ? Rcpp::as<arma::vec>(x) : arma::vec();
}

NarrativeRestrictions fNarrativePrep_cpp(const arma::uvec& ns_shock,
                                         const arma::uvec& ns_period,
                                         const arma::vec&  ns_sign,
                                         const arma::uvec& nd_shock,
                                         const arma::uvec& nd_period,
                                         const arma::uvec& nd_var,
                                         arma::uword N,
                                         arma::uword nobs) {

    if (ns_shock.n_elem != ns_period.n_elem || ns_shock.n_elem != ns_sign.n_elem) {
        Rcpp::stop("narr_sign: shock, period and sign must have the same length.");
    }
    if (nd_shock.n_elem != nd_period.n_elem || nd_shock.n_elem != nd_var.n_elem) {
        Rcpp::stop("narr_dom: shock, period and var must have the same length.");
    }

    NarrativeRestrictions R;
    R.active = (ns_shock.n_elem + nd_shock.n_elem) > 0u;
    if (!R.active) return R;

    auto to0 = [&](const arma::uvec& v, const char* what, arma::uword hi) {
        arma::uvec out(v.n_elem);
        for (arma::uword k = 0; k < v.n_elem; ++k) {
            if (v(k) < 1u || v(k) > hi) {
                Rcpp::stop("%s index %u is out of range (allowed: 1..%u).",
                           what, static_cast<unsigned>(v(k)), static_cast<unsigned>(hi));
            }
            out(k) = v(k) - 1u;
        }
        return out;
    };

    R.ns_shock  = to0(ns_shock,  "narr_sign shock",  N);
    R.ns_period = to0(ns_period, "narr_sign period", nobs);
    R.ns_sign   = ns_sign;
    R.nd_shock  = to0(nd_shock,  "narr_dom shock",   N);
    R.nd_period = to0(nd_period, "narr_dom period",  nobs);
    R.nd_var    = to0(nd_var,    "narr_dom var",     N);

    for (arma::uword k = 0; k < R.ns_sign.n_elem; ++k) {
        if (R.ns_sign(k) == 0.0) Rcpp::stop("narr_sign sign entries must be +1 or -1.");
    }

    arma::uvec all = arma::join_cols(R.ns_period, R.nd_period);
    R.periods = arma::unique(all);          // arma::unique returns ascending

    R.ns_slot.set_size(R.ns_period.n_elem);
    for (arma::uword k = 0; k < R.ns_period.n_elem; ++k)
        R.ns_slot(k) = slot_of(R.periods, R.ns_period(k));
    R.nd_slot.set_size(R.nd_period.n_elem);
    for (arma::uword k = 0; k < R.nd_period.n_elem; ++k)
        R.nd_slot(k) = slot_of(R.periods, R.nd_period(k));

    return R;
}

bool fNarrativeCheck_cpp(const arma::mat&             B,
                         const arma::mat&             resid,
                         const NarrativeRestrictions& R,
                         NarrativeScratch&            scr) {

    if (!R.active) return true;

    const arma::uword N  = B.n_rows;
    const arma::uword ns = R.periods.n_elem;

    if (!arma::inv(scr.Binv, B)) return false;   // singular draw: reject

    if (scr.E.n_rows != N || scr.E.n_cols != ns) scr.E.set_size(N, ns);
    for (arma::uword s = 0; s < ns; ++s) {
        scr.E.col(s) = scr.Binv * resid.row(R.periods(s)).t();
    }

    return restrictions_hold(B, R, scr.E);
}

double fNarrativeWeight_cpp(const arma::mat&             B,
                            const NarrativeRestrictions& R,
                            int                          n_mc,
                            tidymacro::RNG&              rng,
                            arma::mat&                   E_scratch) {

    if (!R.active || n_mc <= 0) return 1.0;

    const arma::uword N  = B.n_rows;
    const arma::uword ns = R.periods.n_elem;
    if (E_scratch.n_rows != N || E_scratch.n_cols != ns) E_scratch.set_size(N, ns);

    int hits = 0;
    for (int m = 0; m < n_mc; ++m) {
        double* e = E_scratch.memptr();
        const arma::uword n_elem = N * ns;
        for (arma::uword i = 0; i < n_elem; ++i) e[i] = rng.norm();
        if (restrictions_hold(B, R, E_scratch)) ++hits;
    }

    // Floor the probability at one Monte Carlo cell so the weight stays finite.
    const double prob = std::max(static_cast<double>(hits), 1.0) /
                        static_cast<double>(n_mc);
    return 1.0 / prob;
}

//' Check Narrative Sign Restrictions for One Draw
//'
//' Evaluates Antolin-Diaz and Rubio-Ramirez (2018) narrative restrictions for a
//' candidate structural impact matrix, and optionally returns the importance
//' weight that makes plain rejection sampling agree with their algorithm.
//'
//' @param B N x N structural impact matrix.
//' @param resid T x N matrix of reduced-form VAR residuals. Row 1 is the first
//'   post-lag observation.
//' @param narr_sign_shock Integer vector (1-indexed) of restricted shocks,
//'   type 1. Default \code{integer(0)}.
//' @param narr_sign_period Integer vector (1-indexed rows of \code{resid}) of
//'   restricted dates, type 1. Default \code{integer(0)}.
//' @param narr_sign_sign Numeric vector of \code{+1} / \code{-1} required signs,
//'   type 1. Default \code{numeric(0)}.
//' @param narr_dom_shock Integer vector (1-indexed) of dominant shocks, type 2.
//'   Default \code{integer(0)}.
//' @param narr_dom_period Integer vector (1-indexed rows of \code{resid}) of
//'   restricted dates, type 2. Default \code{integer(0)}.
//' @param narr_dom_var Integer vector (1-indexed) of variables whose unexpected
//'   movement the shock must dominate, type 2. Default \code{integer(0)}.
//' @param n_mc Integer number of Monte Carlo replications used to estimate the
//'   importance weight. \code{0} (default) skips the estimate and returns a
//'   weight of 1.
//' @param seed Integer seed for the Monte Carlo replications (default 42).
//'
//' @return A list with \code{pass} (logical: all restrictions hold) and
//'   \code{weight} (the ADRR importance weight, or 1 when \code{n_mc = 0}).
//'
//' @details
//' Type 1 requires \code{sign * e[t, j] > 0} for the structural shocks
//' \eqn{e_t = B^{-1} u_t}. Type 2 requires shock \code{j} to contribute more to
//' the unexpected movement in variable \code{i} at date \code{t} than all other
//' shocks combined. The importance weight is the reciprocal of the probability
//' that the restrictions hold when shocks are drawn from their unconditional
//' \eqn{N(0, I)} distribution; the VAR Toolbox omits it, which makes plain
//' rejection sampling an approximation to the ADRR posterior.
//'
//' @references
//' Antolin-Diaz, J., & Rubio-Ramirez, J. F. (2018). Narrative sign restrictions
//' for SVARs. \emph{American Economic Review}, 108(10), 2802--2829.
//'
//' @seealso \code{\link{fSR_cpp}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fCheckNarrative_cpp(const arma::mat& B,
                               const arma::mat& resid,
                               Rcpp::Nullable<arma::uvec> narr_sign_shock = R_NilValue,
                               Rcpp::Nullable<arma::uvec> narr_sign_period = R_NilValue,
                               Rcpp::Nullable<arma::vec>  narr_sign_sign = R_NilValue,
                               Rcpp::Nullable<arma::uvec> narr_dom_shock = R_NilValue,
                               Rcpp::Nullable<arma::uvec> narr_dom_period = R_NilValue,
                               Rcpp::Nullable<arma::uvec> narr_dom_var = R_NilValue,
                               int n_mc = 0,
                               int seed = 42) {

    NarrativeRestrictions R = fNarrativePrep_cpp(
        as_uvec_or_empty(narr_sign_shock), as_uvec_or_empty(narr_sign_period),
        as_vec_or_empty(narr_sign_sign),
        as_uvec_or_empty(narr_dom_shock), as_uvec_or_empty(narr_dom_period),
        as_uvec_or_empty(narr_dom_var),
        B.n_rows, resid.n_rows);

    NarrativeScratch scr;
    const bool pass = fNarrativeCheck_cpp(B, resid, R, scr);

    double weight = 1.0;
    if (n_mc > 0) {
        tidymacro::RNG rng(static_cast<std::uint64_t>(seed));
        arma::mat E;
        weight = fNarrativeWeight_cpp(B, R, n_mc, rng, E);
    }

    return Rcpp::List::create(Rcpp::Named("pass")   = pass,
                              Rcpp::Named("weight") = weight);
}
