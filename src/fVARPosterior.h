#ifndef FVARPOSTERIOR_H
#define FVARPOSTERIOR_H

#include <RcppArmadillo.h>
#include "rng_tidymacro.h"

// Normal-inverse-Wishart posterior for a VAR under a flat (diffuse) prior.
//
//   sigma | y      ~ IW( nobs * sigma_hat , nobs )
//   vec(beta) | .  ~ N( vec(beta_hat) , kron(sigma, (X'X)^{-1}) )
//
// which is the sampler in VARdrawpost.m of Cesa-Bianchi's VAR Toolbox.
//
// The Kronecker covariance is never formed.  With Lx Lx' = (X'X)^{-1} and
// G G' = sigma_draw, drawing Z (k x N) i.i.d. standard normal and setting
//   beta_draw = beta_hat + Lx * Z * G'
// reproduces exactly the required covariance at O(k^2 N + k N^2) cost and
// O(k^2) memory, instead of the O(k^2 N^2) of an explicit kron.

// Reusable factorisation: built once, then sampled many times.
struct NIWPosterior {
    arma::mat beta_hat;   // k x N  OLS coefficients (posterior mean)
    arma::mat Lsig;       // N x N  lower Cholesky of nobs * sigma_hat
    arma::mat Lx;         // k x k  upper-triangular factor, Lx * Lx' = (X'X)^{-1}
    int nobs;             // residual-sample length (posterior degrees of freedom)
    int N;                // number of endogenous variables
    int k;                // number of regressors per equation
};

// Scratch buffers reused across draws by one worker (no allocation in the loop).
struct NIWScratch {
    arma::mat A;      // N x N  Bartlett factor
    arma::mat Ainv;   // N x N
    arma::mat Z;      // k x N
    arma::mat tmp;    // k x N
};

// Build the factorisation from the OLS design.  Throws if nobs <= N.
NIWPosterior fNIWPosteriorPrep_cpp(const arma::mat& beta_hat,
                                   const arma::mat& sigma_hat,
                                   const arma::mat& XtX,
                                   int nobs);

// One posterior draw.  Writes sigma_draw (N x N), its lower factor G
// (sigma_draw = G G') and beta_draw (k x N).  Thread safe: all state is in
// `rng` and `scratch`.
void fNIWPosteriorDraw_cpp(const NIWPosterior& post,
                           tidymacro::RNG&     rng,
                           NIWScratch&         scratch,
                           arma::mat&          G,
                           arma::mat&          sigma_draw,
                           arma::mat&          beta_draw);

// Build the VAR design matrices Y = X * beta + u, matching fVAR's layout:
// X = [const | y_{t-1} ... y_{t-p} | exog_t].
void fVARDesign_cpp(const arma::mat& y, int p, int c, const arma::mat& exog,
                    arma::mat& Y, arma::mat& X);

#endif // FVARPOSTERIOR_H
