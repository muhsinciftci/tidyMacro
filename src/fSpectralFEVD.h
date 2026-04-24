#ifndef FSPECTRALFEVD_H
#define FSPECTRALFEVD_H

#include <RcppArmadillo.h>

arma::vec fSpectralFEVD_cpp(const arma::cube& D, const arma::mat& irf_s,
                              const arma::vec& band, int J, bool fourier);

#endif
