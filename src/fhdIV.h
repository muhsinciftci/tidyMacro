#ifndef FHDIV_H
#define FHDIV_H

#include <RcppArmadillo.h>
#include "fVAR.h"

// Struct for point-estimate companion-form HD
struct HDIVResult {
    arma::mat HDshock;  // T_eff x N matrix (contribution of IV shock to each variable)
};

// Struct for bootstrap HD confidence bands
struct BootstrapHDIVResult {
    arma::mat HDshock;  // T_eff x N point-estimate HD (recentering anchor)
    arma::mat upper;    // T_eff x N upper confidence bands
    arma::mat lower;    // T_eff x N lower confidence bands
};

// Internal C++ functions
HDIVResult fhdIV_cpp(const arma::mat& residuals,
                     const arma::mat& sigma,
                     const arma::vec& s,
                     const arma::mat& beta,
                     int c, int p);

BootstrapHDIVResult fbootstrapHDIV_cpp(const arma::mat& y,
                                        const VARResult& var_result,
                                        const arma::mat& Z,
                                        const arma::vec& s,
                                        int nboot, int blocksize,
                                        const arma::ivec& adjustZ,
                                        const arma::ivec& adjustu,
                                        int policyvar,
                                        double prc,
                                        int n_threads);

// R wrapper functions
Rcpp::List fhdIV(const arma::mat& residuals,
                 const arma::mat& sigma,
                 const arma::vec& s,
                 const arma::mat& beta,
                 int c, int p);

Rcpp::List fbootstrapHDIV(const arma::mat& y,
                           const Rcpp::List& var_result,
                           const arma::mat& Z,
                           const arma::vec& s,
                           int nboot, int blocksize,
                           const arma::ivec& adjustZ,
                           const arma::ivec& adjustu,
                           int policyvar,
                           double prc,
                           int n_threads);

#endif // FHDIV_H
