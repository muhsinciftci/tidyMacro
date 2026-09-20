#ifndef FSIGNRESTRICTIONS_H
#define FSIGNRESTRICTIONS_H

#include <RcppArmadillo.h>
#include "rng_tidymacro.h"

// Search for an orthonormal rotation of a reduced-form covariance matrix whose
// implied impact (or short-horizon) responses satisfy a sign-restriction
// pattern.  Port of SignRestrictions.m from Cesa-Bianchi's VAR Toolbox 4.0.
//
// Column convention.  `Bfix` holds columns of B that are already identified by
// another scheme (an external instrument under ident = "sign+iv"); those
// columns are held fixed and only the remaining ones are rotated.  Pass an
// empty matrix for plain sign restrictions.

// Per-worker scratch.  Everything that would otherwise be allocated inside the
// rotation loop lives here and is resized at most once per posterior draw.
struct SignRotScratch {
    arma::mat  C;            // N x N   lower Cholesky of sigma
    arma::mat  startingMat;  // N x N   fixed columns + orthonormal completion
    arma::mat  rotated;      // N x m   rotated free block
    arma::mat  termaa;       // N x N   candidate B before column reordering
    arma::mat  Qs, Rs, Gs;   //         Haar-rotation scratch
    arma::cube irfchk;       // N x m x sr_hor  responses used by the sign check
    std::vector<char>       used;   // free column already matched to a shock
    std::vector<arma::uword> order; // column permutation applied at the end
};

// Prepare the parts of the search that depend only on (sigma, Bfix): the
// Cholesky factor and the orthonormal completion.  Call once per posterior
// draw, then call fSignRotation_cpp repeatedly.
void fSignRotationPrep_cpp(const arma::mat& sigma,
                           const arma::mat& Bfix,
                           tidymacro::RNG&  rng,
                           SignRotScratch&  scr);

// One rotation search.  Returns true and fills `B_out` when a rotation
// satisfying SIGN is found within `sr_rot` attempts; `n_tried` always reports
// how many rotations were drawn.  `wold` supplies the Wold multipliers of the
// current draw and is read only for slices 0 .. sr_hor-1 (slice 0 must be the
// identity); it is ignored when sr_hor == 1.
bool fSignRotation_cpp(const arma::mat&  SIGN,
                       const arma::cube& wold,
                       int               sr_hor,
                       int               sr_rot,
                       int               n_fixed,
                       tidymacro::RNG&   rng,
                       SignRotScratch&   scr,
                       arma::mat&        B_out,
                       int&              n_tried);

#endif // FSIGNRESTRICTIONS_H
