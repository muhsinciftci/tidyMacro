#ifndef FCHECKRESTRICTIONS_H
#define FCHECKRESTRICTIONS_H

#include <RcppArmadillo.h>

// Check whether the IRF column for one shock satisfies sign restrictions.
//
// shock_idx : 0-indexed shock position (column of Q / column of irf)
// restr     : length-N row vector of signs (-1, 0, +1) for each variable
// hor_vec   : 0-indexed horizon indices to check (length rr)
//
// Returns
//   +1  all non-generic restrictions pass at every checked horizon
//    0  all non-generic restrictions fail at every checked horizon
//         (flipping the sign of this column fixes them all)
//   -1  mixed outcome across horizons or variables; draw must be discarded
int fCheckRestrictions_cpp(const arma::cube& irf,
                           int               shock_idx,
                           const arma::rowvec& restr,
                           const arma::uvec&   hor_vec);

// R wrapper (shock: 1-indexed; hor_vec: 1-indexed horizon numbers)
int fCheckRestrictions(const arma::cube&   irf,
                       int                 shock,
                       const arma::rowvec& restr,
                       const arma::uvec&   hor_vec);

#endif // FCHECKRESTRICTIONS_H
