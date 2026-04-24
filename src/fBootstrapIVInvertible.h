#ifndef FBOOTSTRAPIVINVERTIBLE_H
#define FBOOTSTRAPIVINVERTIBLE_H

#include <RcppArmadillo.h>
#include "fVAR.h"

arma::cube fBootstrapIVInvertible_cpp(
    const arma::mat& y,
    const arma::vec& instr,
    const VARResult& var_result,
    int nboot, int p, int c, int hor,
    const arma::ivec& cumu);

#endif
