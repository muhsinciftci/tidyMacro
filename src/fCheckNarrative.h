#ifndef FCHECKNARRATIVE_H
#define FCHECKNARRATIVE_H

#include <RcppArmadillo.h>
#include "rng_tidymacro.h"

// Narrative sign restrictions of Antolin-Diaz & Rubio-Ramirez (2018).
//
//   Type 1 (narr_sign): the structural shock j at date t* has a given sign,
//                       sign_k * e(t*, j) > 0.
//   Type 2 (narr_dom):  shock j is the overwhelming driver of the unexpected
//                       movement in variable i at date t*,
//                       |B(i,j) e(t*,j)| > sum_{k != j} |B(i,k) e(t*,k)|.
//
// Structural shocks are e_t = B^{-1} u_t.  Only the dated rows are ever
// formed: solve B * E = dated residuals for the handful of required dates,
// rather than constructing B inverse or solving for the whole sample.

// Restrictions in 0-based internal form, with the dated residual rows
// de-duplicated so each is transformed only once per draw.
struct NarrativeRestrictions {
    arma::uvec ns_shock, ns_period, ns_slot;   // Type 1
    arma::vec  ns_sign;
    arma::uvec nd_shock, nd_period, nd_var, nd_slot;   // Type 2
    arma::uvec periods;        // unique residual rows, ascending
    bool active;
};

// Per-worker scratch for the narrative check.
struct NarrativeScratch {
    arma::mat dated;  // N x n_periods, selected residual vectors
    arma::mat E;      // N x n_periods, column s = shocks at periods(s)
};

// Build the internal representation.  Period indices arrive 1-based (rows of
// the residual sample); shock, variable and sign vectors must be the same
// length within each restriction type.
NarrativeRestrictions fNarrativePrep_cpp(const arma::uvec& ns_shock,
                                         const arma::uvec& ns_period,
                                         const arma::vec&  ns_sign,
                                         const arma::uvec& nd_shock,
                                         const arma::uvec& nd_period,
                                         const arma::uvec& nd_var,
                                         arma::uword N,
                                         arma::uword nobs);

// True when every narrative restriction holds for this draw.
bool fNarrativeCheck_cpp(const arma::mat&             B,
                         const arma::mat&             resid,
                         const NarrativeRestrictions& R,
                         NarrativeScratch&            scr);

// Antolin-Diaz & Rubio-Ramirez importance weight 1 / Pr(narrative restrictions
// hold), estimated by Monte Carlo under the unconditional shock distribution
// e ~ N(0, I).  Returns 1 when no restrictions are active.
double fNarrativeWeight_cpp(const arma::mat&             B,
                            const NarrativeRestrictions& R,
                            int                          n_mc,
                            tidymacro::RNG&              rng,
                            arma::mat&                   E_scratch);

#endif // FCHECKNARRATIVE_H
