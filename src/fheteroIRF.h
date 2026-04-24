#ifndef FHETERORF_H
#define FHETERORF_H

#include "fVAR.h"
#include <RcppArmadillo.h>

// Point-estimate result
struct HeteroIRFResult {
    arma::mat IRF;  // N x (hor+1) impulse responses
    arma::vec b1;   // N structural impact vector (normalized)
};

// Bootstrap result
struct BootHeteroResult {
    arma::mat upper;     // N x (hor+1)
    arma::mat lower;     // N x (hor+1)
    arma::mat upper2;    // N x (hor+1)
    arma::mat lower2;    // N x (hor+1)
    arma::mat meanirf;   // N x (hor+1)
    arma::mat medianirf; // N x (hor+1)
    arma::mat point;     // N x (hor+1) plug-in point estimate (center of Hall bands)
};

// Internal C++ point estimate
// AL     : N x N*p VAR coefficients (no constant)
// eta    : T1 x N residuals in proxy window
// Z      : T1 x 1 proxy in proxy window
// indsR1 : T1 integer 0/1 (1 = OPEC announcement month)
// nvar   : 1-based normalization variable
// scale  : shock size
HeteroIRFResult fheteroIRF_cpp(const arma::mat& AL,
                                const arma::mat& eta,
                                const arma::mat& Z,
                                const arma::ivec& indsR1,
                                int p, int hor, int nvar, double scale);

// Internal C++ bootstrap
BootHeteroResult fbootstrapHetero_cpp(const arma::mat& y,
                                      const VARResult& var_result,
                                      const arma::mat& Z,
                                      const arma::ivec& indsR1,
                                      const arma::ivec& adjustu,
                                      int nboot, int blocksize,
                                      int hor, int nvar, double scale,
                                      double prc, double prc2,
                                      int n_threads);

// R wrapper: point estimate
Rcpp::List fheteroIRF(const Rcpp::List& var_result,
                      const arma::mat& Z,
                      const arma::ivec& adjustu,
                      const arma::ivec& indsR1,
                      int hor, int nvar, double scale);

// R wrapper: bootstrap
Rcpp::List fbootstrapHetero(const arma::mat& y,
                             const Rcpp::List& var_result,
                             const arma::mat& Z,
                             const arma::ivec& indsR1,
                             const arma::ivec& adjustu,
                             int nboot, int blocksize,
                             int hor, int nvar, double scale,
                             double prc, double prc2, int n_threads);

#endif // FHETERORF_H
