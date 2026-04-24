#ifndef FBOOTSTRAPIVRECOVER_H
#define FBOOTSTRAPIVRECOVER_H

#include <RcppArmadillo.h>
#include "fVAR.h"

arma::cube fBootstrapIVRecover_cpp(
    const arma::mat& y,
    const arma::vec& instr,
    const VARResult& var_result,
    const arma::vec& noise,
    const arma::vec& delta,
    int nboot, int p, int c, int r, int hor,
    const arma::ivec& cumu);

#endif
