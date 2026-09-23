// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_DONT_USE_OPENMP

#include "fRecoverBIVMulti.h"
#include <RcppArmadillo.h>

namespace {

// Port of inv2.m: m \ eye(size(m)), with entries below 1e-16 zeroed.
arma::mat inv2_cpp(const arma::mat& m) {
    arma::mat out = arma::solve(m, arma::eye<arma::mat>(m.n_rows, m.n_cols));
    out.elem(arma::find(arma::abs(out) < 1e-16)).zeros();
    return out;
}

// chol(A, 'lower'), falling back to an eigenvalue-clipped PSD projection when
// A is not (numerically) positive definite - same fallback fCompleteB_cpp uses.
arma::mat chol_lower_safe(const arma::mat& A) {
    arma::mat L;
    if (arma::chol(L, A, "lower")) return L;
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);
    eigval.elem(arma::find(eigval < 0.0)).fill(1e-10);
    arma::mat A_psd = eigvec * arma::diagmat(eigval) * eigvec.t();
    if (!arma::chol(L, A_psd, "lower")) {
        Rcpp::stop("Could not Cholesky-factorise S1S1p even after PSD projection.");
    }
    return L;
}

} // namespace

//' Proxy-SVAR Impact Columns from k Jointly-Identified External Instruments
//'
//' Recovers the structural impact columns of k >= 1 instrumented shocks by
//' the block generalisation of the Mertens-Ravn procedure used in
//' \code{VARiriv_B.m} (Cesa-Bianchi & Sokol 2022 replication code): first-
//' and second-stage regressions on the joint instrument set, followed by the
//' block partition of the instrument-subsample covariance into
//' \code{S11}/\code{S21}/\code{S22} blocks that pins down \code{b11} and
//' \code{b21}. With \code{k = 1} this collapses to the same regressions as
//' \code{\link{fRecoverBIV_cpp}}, but without that function's extra
//' Gertler-Karadi sign/scale normalisation, so the two need not agree
//' numerically at \code{k = 1}.
//'
//' @param resid_sub T x N matrix of reduced-form VAR residuals restricted to
//'   the rows overlapping the instruments. The instrumented variables must be
//'   the first k columns, in the same order as \code{Z_sub}.
//' @param Z_sub T x k matrix of instruments, row-aligned with \code{resid_sub}.
//' @param ntotcoeff Integer number of coefficients per VAR equation
//'   (\code{N * p + c + n_exog}), used in the degrees-of-freedom correction of
//'   the instrument-subsample covariance.
//'
//' @return A list with \code{B} (N x k identified impact columns, ready to
//'   pass as \code{Bfix} to \code{\link{fSR_cpp}}), \code{sigma_b}
//'   (instrument-subsample covariance), \code{fs_beta}, \code{fs_F} and
//'   \code{fs_r2} (first-stage coefficients, F statistics and R-squared, one
//'   per instrumented equation), \code{relEig} (eigenvalues of the
//'   reliability matrix, descending) and \code{n_iv}.
//'
//' @references
//' Mertens, K., & Ravn, M. O. (2013). The dynamic effects of personal and
//' corporate income tax changes in the United States. \emph{American Economic
//' Review}, 103(4), 1212--1247.
//'
//' Cesa-Bianchi, A., & Sokol, A. (2022). Financial shocks, credit spreads,
//' and the international credit channel. \emph{Journal of International
//' Economics}, 135, 103543.
//'
//' @seealso \code{\link{fSR_cpp}}, \code{\link{fRecoverBIV_cpp}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fRecoverBIVMulti_cpp(const arma::mat& resid_sub,
                                const arma::mat& Z_sub,
                                int ntotcoeff) {

    IVColumnsResult r = fIVColumns_cpp(resid_sub, Z_sub, ntotcoeff);

    return Rcpp::List::create(Rcpp::Named("B")       = r.B,
                              Rcpp::Named("sigma_b") = r.sigma_b,
                              Rcpp::Named("fs_beta") = r.fs_beta,
                              Rcpp::Named("fs_F")    = r.fs_F,
                              Rcpp::Named("fs_r2")   = r.fs_r2,
                              Rcpp::Named("relEig")  = r.relEig,
                              Rcpp::Named("n_iv")    = r.n_iv);
}

IVColumnsResult fIVColumns_cpp(const arma::mat& resid_sub,
                               const arma::mat& Z_sub,
                               int ntotcoeff) {

    const arma::uword T = resid_sub.n_rows;
    const arma::uword N = resid_sub.n_cols;
    const arma::uword k = Z_sub.n_cols;

    if (Z_sub.n_rows != T) Rcpp::stop("resid_sub and Z_sub must have the same number of rows.");
    if (k == 0)             Rcpp::stop("Z_sub must contain at least one instrument.");
    if (N <= k)             Rcpp::stop("The VAR must contain more variables than instrumented shocks.");
    if (ntotcoeff < 0)      Rcpp::stop("ntotcoeff must be non-negative.");
    if (!resid_sub.is_finite() || !Z_sub.is_finite()) {
        Rcpp::stop("resid_sub and Z_sub must contain only finite values.");
    }
    if (static_cast<int>(T) <= ntotcoeff + static_cast<int>(k)) {
        Rcpp::stop("The instrument subsample (%d rows) is too short for %d VAR coefficients "
                   "and %d instrumented shocks.", static_cast<int>(T), ntotcoeff, static_cast<int>(k));
    }

    arma::mat u1 = resid_sub.cols(0, k - 1);             // T x k   "up"
    arma::mat u2 = resid_sub.cols(k, N - 1);             // T x (N-k) "uq"
    const arma::mat& m = Z_sub;                          // T x k

    // VARiriv_B.m explicitly removes the residual means before both stages.
    // This matters for its reported first-stage and reliability diagnostics.
    u1.each_row() -= arma::mean(u1, 0);
    u2.each_row() -= arma::mean(u2, 0);

    // ---- first stage: u1 on [1, m] --------------------------------------
    arma::mat X1 = arma::join_rows(arma::ones<arma::mat>(T, 1), m);
    arma::mat fs_beta = arma::solve(X1, u1);              // (1+k) x k
    arma::mat u1_hat  = X1 * fs_beta;                     // T x k

    // ---- second stage: u2 on u1_hat, without an intercept -----------------
    // This is MATLAB's `u1Hat \\ u2` exactly.
    arma::mat invSmu1Smu2 = arma::solve(u1_hat, u2);      // k x (N-k)

    // ---- instrument-subsample covariance, df-corrected -------------------
    arma::mat pq_dm  = resid_sub.each_row() - arma::mean(resid_sub, 0);
    arma::mat sigma_b = (pq_dm.t() * pq_dm) / (static_cast<double>(T) - ntotcoeff);

    const arma::mat S11 = sigma_b.submat(0, 0, k - 1, k - 1);
    const arma::mat S21 = sigma_b.submat(k, 0, N - 1, k - 1);
    const arma::mat S22 = sigma_b.submat(k, k, N - 1, N - 1);

    // ---- block algebra (VARiriv_B.m, generalised to k >= 1) --------------
    const arma::mat b21invb11  = invSmu1Smu2.t();         // (N-k) x k
    const arma::mat b21invb11p = invSmu1Smu2;             // k x (N-k)

    const arma::mat Zmat = (b21invb11 * S11) * b21invb11p
                          - (S21 * b21invb11p + b21invb11 * S21.t())
                          + S22;                                    // (N-k) x (N-k)

    const arma::mat d = S21 - b21invb11 * S11;                      // (N-k) x k
    const arma::mat b12b12p = d.t() * inv2_cpp(Zmat) * d;           // k x k

    const arma::mat b22b22p = S22 + (b21invb11 * (b12b12p - S11)) * b21invb11p; // (N-k)x(N-k)

    const arma::mat b12invb22 = (b12b12p * b21invb11p + d.t()) * inv2_cpp(b22b22p); // k x (N-k)

    const arma::mat b11b11p = S11 - b12b12p;                        // k x k

    const arma::mat I_k = arma::eye<arma::mat>(k, k);
    const arma::mat M   = I_k - b12invb22 * b21invb11;              // k x k

    const arma::mat b11invS1 = inv2_cpp(M);                         // k x k
    const arma::mat b21invS1 = b21invb11 * b11invS1;                // (N-k) x k

    const arma::mat S1S1p = (M * b11b11p) * M.t();                  // k x k
    const arma::mat S1    = chol_lower_safe(S1S1p);                 // k x k

    const arma::mat b11 = b11invS1 * S1;                            // k x k
    const arma::mat b21 = b21invS1 * S1;                            // (N-k) x k

    arma::mat B(N, k, arma::fill::none);
    B.rows(0, k - 1)   = b11;
    B.rows(k, N - 1)   = b21;

    // ---- reliability statistics (Mertens & Ravn 2013, Table p.1229) ------
    const arma::mat u1_res = u1 - u1_hat;
    arma::uvec nonCensored = arma::find(arma::sum(arma::abs(m), 1) > 0.0);
    if (nonCensored.n_elem <= k + 1) {
        Rcpp::stop("Too few non-zero instrument observations for first-stage diagnostics.");
    }
    const double dfrac = static_cast<double>(nonCensored.n_elem) / static_cast<double>(T);

    const arma::mat Smm  = (m.t() * m) / static_cast<double>(T);           // k x k
    const arma::mat Smu1 = (m.t() * u1) / static_cast<double>(T);          // k x k
    const arma::mat Lambda = inv2_cpp(Smm) * Smu1 * (inv2_cpp(b11b11p) * Smu1.t()) / dfrac;

    arma::vec relEig = arma::sort(arma::real(arma::eig_gen(Lambda)), "descend");

    const arma::mat u1_nc  = u1.rows(nonCensored);
    const arma::mat res_nc = u1_res.rows(nonCensored);
    const arma::mat ss_tot = u1_nc.t() * u1_nc;
    const arma::mat ss_res = res_nc.t() * res_nc;
    const double n_nc = static_cast<double>(nonCensored.n_elem);

    arma::vec fs_r2(k), fs_F(k);
    for (arma::uword i = 0; i < k; ++i) {
        fs_r2(i) = 1.0 - ss_res(i, i) / ss_tot(i, i);
        fs_F(i)  = ((ss_tot(i, i) - ss_res(i, i)) / static_cast<double>(k)) /
                   (ss_res(i, i) / (n_nc - static_cast<double>(k) - 1.0));
    }

    IVColumnsResult out;
    out.B       = B;
    out.sigma_b = sigma_b;
    out.fs_beta = fs_beta;
    out.fs_F    = fs_F;
    out.fs_r2   = fs_r2;
    out.relEig  = relEig;
    out.n_iv    = static_cast<int>(T);
    return out;
}
