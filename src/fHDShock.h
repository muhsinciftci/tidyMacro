#ifndef FHDSHOCK_H
#define FHDSHOCK_H

#include <RcppArmadillo.h>

// Historical decomposition of a structurally identified VAR: the contribution
// of each structural shock to each series, plus the initial condition, the
// intercept and the exogenous regressors.  Port of compute_HD.m from the
// VAR Toolbox, without the linear-trend block (fVAR carries an intercept only).

struct HDShockResult {
    arma::cube shock;   // T x N x N   [time, variable, shock]
    arma::mat  init;    // T x N       initial condition
    arma::mat  cons;    // T x N       intercept
    arma::cube exo;     // T x N x M   exogenous regressors
    arma::mat  endo;    // T x N       sum of all components
};

HDShockResult fHDCompute_cpp(const arma::mat& y,
                             const arma::mat& beta,
                             const arma::mat& B,
                             int p, int c,
                             const arma::mat& exog);

#endif // FHDSHOCK_H
