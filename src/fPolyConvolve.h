#ifndef FPOLYCONVOLVE_H
#define FPOLYCONVOLVE_H

#include <RcppArmadillo.h>

arma::cube fPolyConvolve_cpp(const arma::cube& A, const arma::cube& B, int nlags);

#endif
