#ifndef FGETBANDS_H
#define FGETBANDS_H

#include <RcppArmadillo.h>

struct FGetBandsResult {
    arma::mat upper;
    arma::mat lower;
    arma::mat median;
};

struct FGetBands2Result {
    arma::mat upper;
    arma::mat lower;
    arma::mat upper2;
    arma::mat lower2;
    arma::mat median;
};

FGetBandsResult  fGetBands_cpp(const arma::cube& bootirf, double prc);
FGetBands2Result fGetBands2_cpp(const arma::cube& bootirf, double prc, double prc2);

#endif
