#ifndef FRECOVERBIVMULTI_H
#define FRECOVERBIVMULTI_H

#include <RcppArmadillo.h>

// Joint external-instrument (proxy SVAR) identification of k >= 1 structural
// impact columns from k instruments, following the block generalisation of
// Mertens & Ravn (2013) used in VARiriv_B.m (Cesa-Bianchi & Sokol 2022
// replication code). Reduces to the single-instrument case of
// fIVColumn_cpp/fRecoverBIV_cpp when k == 1, but is not required to match it
// bit-for-bit since VARiriv_B.m does not apply the extra Gertler-Karadi
// sign/scale normalisation that fIVColumn_cpp adds for k == 1.
//
// The first k columns of resid_sub are the instrumented ("up") residuals, in
// the same order as the k columns of Z_sub; columns k+1..N are the remaining
// ("uq") residuals identified elsewhere (sign restrictions).

struct IVColumnsResult {
    arma::mat B;         // N x k  identified impact columns (b11 stacked on b21)
    arma::mat sigma_b;   // N x N  residual covariance on the instrument subsample (df-corrected)
    arma::mat fs_beta;   // (1+k) x k  first-stage coefficients (intercept row first)
    arma::vec fs_F;      // k x 1  first-stage F statistic per instrumented equation
    arma::vec fs_r2;     // k x 1  first-stage R-squared per instrumented equation
    arma::vec relEig;    // k x 1  eigenvalues of the reliability matrix, descending
    int       n_iv;      // instrument-subsample length
};

// resid_sub and Z_sub must already be row-aligned on the (contiguous)
// instrument subsample, as fSignRestr() aligns them before calling this.
IVColumnsResult fIVColumns_cpp(const arma::mat& resid_sub,
                               const arma::mat& Z_sub,
                               int ntotcoeff);

#endif // FRECOVERBIVMULTI_H
