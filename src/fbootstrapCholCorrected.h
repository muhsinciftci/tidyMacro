#ifndef FBOOTSTRAPCHOLCORRECTED_H
#define FBOOTSTRAPCHOLCORRECTED_H

#include "fVAR.h"
#include "fbootstrapChol.h"
#include <RcppArmadillo.h>

// Struct to hold bias-corrected bootstrap results
// Extends BootstrapCholResult with bias diagnostics
struct BootstrapCholCorrectedResult {
    arma::mat bootchol_flat; // slice_sz x nboot2 flattened bootstrap IRFs
    arma::cube upper;        // N x N x (horizon+1) upper confidence bands
    arma::cube lower;        // N x N x (horizon+1) lower confidence bands
    arma::cube boot_beta;    // N x n_coef x nboot2 bootstrapped coefficients
    arma::mat  Beta;         // N x n_coef bias-corrected coefficient matrix
    int        corrections;  // Number of iterative shrinkage steps applied
    int        N;
    int        H;
};

// Internal C++ function
BootstrapCholCorrectedResult
fbootstrapCholCorrected_cpp(const arma::mat&    y,
                            const VARResult&    var_result,
                            int                 nboot1,
                            int                 nboot2,
                            int                 horizon,
                            double              prc,
                            const std::string&  bootscheme,
                            Rcpp::Nullable<arma::mat> exog,
                            int                 n_threads);

// R wrapper
Rcpp::List fbootstrapCholCorrected(const arma::mat&    y,
                                   const Rcpp::List&   var_result,
                                   int                 nboot1,
                                   int                 nboot2,
                                   int                 horizon,
                                   double              prc,
                                   const std::string&  bootscheme,
                                   Rcpp::Nullable<arma::mat> exog,
                                   int                 n_threads);

#endif // FBOOTSTRAPCHOLCORRECTED_H
