#ifndef FRECOVERBIV_H
#define FRECOVERBIV_H

#include <RcppArmadillo.h>

// External-instrument (proxy SVAR) identification of the first structural
// impact column, following Mertens & Ravn (2013) and Gertler & Karadi (2015)
// as implemented in recover_B.m of the VAR Toolbox, plus the Cholesky-QR
// completion of the remaining columns.

struct IVColumnResult {
    arma::vec b1;        // N x 1  identified impact column
    arma::mat B;         // N x N  completed impact matrix, B B' = sigma
    arma::mat sigma_b;   // N x N  residual covariance on the instrument subsample
    arma::vec fs_beta;   // first-stage coefficients (intercept first)
    double    fs_F;      // first-stage F statistic on the excluded instruments
    double    fs_r2;     // first-stage R-squared
    double    shock_sd;  // sp, the implied standard deviation of the shock
    int       n_iv;      // instrument-subsample length
};

// Complete one identified impact column to a full invertible B with
// B B' = sigma exactly and B(:,1) = b1 exactly.  Columns 2:N carry no
// economic content; they exist so that B can be inverted.
arma::mat fCompleteB_cpp(const arma::vec& b1, const arma::mat& sigma);

// resid_sub and Z_sub must already be row-aligned on the instrument subsample.
// The instrumented variable is the first column of resid_sub.
IVColumnResult fIVColumn_cpp(const arma::mat& resid_sub,
                             const arma::mat& Z_sub,
                             const arma::mat& sigma,
                             int ntotcoeff);

#endif // FRECOVERBIV_H
