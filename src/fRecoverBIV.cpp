// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_DONT_USE_OPENMP

#include "fRecoverBIV.h"
#include <RcppArmadillo.h>

arma::mat fCompleteB_cpp(const arma::vec& b1, const arma::mat& sigma) {

    const arma::uword N = sigma.n_rows;
    if (b1.n_elem != N) Rcpp::stop("b1 and sigma have incompatible dimensions.");

    arma::mat P;
    if (!arma::chol(P, sigma, "lower")) {
        Rcpp::stop("The reduced-form covariance matrix is not positive definite.");
    }

    // b1 in the Cholesky basis.  Its norm is only approximately one because b1
    // is calibrated on the instrument subsample while sigma uses the full
    // sample, so normalise before completing the basis.
    arma::vec q1 = arma::solve(arma::trimatl(P), b1);
    const double nrm = arma::norm(q1, 2);
    if (nrm < 1e-12) Rcpp::stop("The IV-identified impact column is numerically zero.");
    q1 /= nrm;

    arma::mat Q, R;
    arma::qr_econ(Q, R, arma::join_rows(q1, arma::eye<arma::mat>(N, N)));
    if (arma::dot(Q.col(0), q1) < 0.0) Q.col(0) *= -1.0;

    arma::mat B = P * Q;    // B B' = P Q Q' P' = sigma, exactly
    B.col(0) = b1;          // restore the exact identified column
    return B;
}

IVColumnResult fIVColumn_cpp(const arma::mat& resid_sub,
                             const arma::mat& Z_sub,
                             const arma::mat& sigma,
                             int ntotcoeff) {

    const arma::uword T = resid_sub.n_rows;
    const arma::uword N = resid_sub.n_cols;
    const arma::uword k = Z_sub.n_cols;

    if (Z_sub.n_rows != T) Rcpp::stop("resid_sub and Z_sub must have the same number of rows.");
    if (N < 2u)            Rcpp::stop("The VAR must contain at least two variables.");
    if (static_cast<int>(T) <= ntotcoeff) {
        Rcpp::stop("The instrument subsample (%d rows) is too short for %d VAR coefficients.",
                   static_cast<int>(T), ntotcoeff);
    }

    const arma::vec p = resid_sub.col(0);
    const arma::mat q = resid_sub.cols(1, N - 1);

    // ---- first stage: p on [1, Z] --------------------------------------
    arma::mat X1 = arma::join_rows(arma::ones<arma::mat>(T, 1), Z_sub);
    arma::vec fs_beta = arma::solve(X1, p);
    arma::vec p_hat   = X1 * fs_beta;

    const arma::vec  fs_res = p - p_hat;
    const double ss_res = arma::dot(fs_res, fs_res);
    const double ss_tot = arma::dot(p - arma::mean(p), p - arma::mean(p));
    const double r2     = (ss_tot > 0.0) ? 1.0 - ss_res / ss_tot : 0.0;
    const double df     = static_cast<double>(T) - static_cast<double>(k) - 1.0;
    const double Fstat  = (r2 < 1.0 && df > 0.0)
                            ? (r2 / (1.0 - r2)) * df / static_cast<double>(k)
                            : arma::datum::nan;

    // ---- second stage: each remaining residual on the fitted p ----------
    arma::vec b1(N, arma::fill::zeros);
    b1(0) = 1.0;
    arma::vec s21s11(N - 1, arma::fill::zeros);
    arma::mat X2 = arma::join_rows(arma::ones<arma::mat>(T, 1), p_hat);
    // One factorisation serves all N-1 equations: the regressor block is common.
    arma::mat coef = arma::solve(X2, q);          // 2 x (N-1)
    for (arma::uword i = 1; i < N; ++i) {
        b1(i)       = coef(1, i - 1);
        s21s11(i-1) = coef(1, i - 1);
    }

    // ---- shock-size normalisation (Gertler-Karadi 2015, fn. 4) ---------
    arma::mat pq_dm = resid_sub.each_row() - arma::mean(resid_sub, 0);
    arma::mat sigma_b = (pq_dm.t() * pq_dm) /
                        (static_cast<double>(T) - static_cast<double>(ntotcoeff));

    const double    S11 = sigma_b(0, 0);
    const arma::vec S21 = sigma_b.submat(1, 0, N - 1, 0);
    const arma::mat S22 = sigma_b.submat(1, 1, N - 1, N - 1);

    arma::mat Qm = s21s11 * S11 * s21s11.t()
                 - (S21 * s21s11.t() + s21s11 * S21.t()) + S22;
    arma::vec d  = S21 - s21s11 * S11;
    const double inner = S11 - arma::as_scalar(d.t() * arma::solve(Qm, d));
    if (!(inner > 0.0)) {
        Rcpp::stop("The IV shock-size normalisation is not positive; the "
                   "instrument is likely too weak for this system.");
    }
    const double sp = std::sqrt(inner);

    b1 *= sp * ((fs_beta(1) < 0.0) ? -1.0 : 1.0);

    IVColumnResult out;
    out.b1       = b1;
    out.B        = fCompleteB_cpp(b1, sigma);
    out.sigma_b  = sigma_b;
    out.fs_beta  = fs_beta;
    out.fs_F     = Fstat;
    out.fs_r2    = r2;
    out.shock_sd = sp;
    out.n_iv     = static_cast<int>(T);
    return out;
}

//' Proxy-SVAR Impact Column from an External Instrument
//'
//' Recovers the structural impact column of the instrumented variable by the
//' Mertens-Ravn / Gertler-Karadi two-stage procedure, and completes it to a
//' full invertible impact matrix.
//'
//' @param resid_sub T x N matrix of reduced-form VAR residuals restricted to
//'   the rows overlapping the instrument. The instrumented variable must be
//'   the first column.
//' @param Z_sub T x k matrix of instruments, row-aligned with \code{resid_sub}.
//' @param sigma N x N full-sample reduced-form residual covariance matrix.
//' @param ntotcoeff Integer number of coefficients per VAR equation
//'   (\code{N * p + c + n_exog}), used in the degrees-of-freedom correction of
//'   the instrument-subsample covariance.
//'
//' @return A list with \code{b1} (identified impact column), \code{B}
//'   (completed N x N impact matrix with \code{B \%*\% t(B) = sigma} and
//'   \code{B[, 1] = b1}), \code{sigma_b} (instrument-subsample covariance),
//'   \code{fs_beta}, \code{fs_F} and \code{fs_r2} (first-stage coefficients,
//'   F statistic on the excluded instruments, and R-squared), \code{shock_sd}
//'   (the implied shock standard deviation) and \code{n_iv}.
//'
//' @details
//' Columns 2 to N of \code{B} are a Cholesky-QR completion with no economic
//' content; they exist so that \code{B} is invertible, which the historical
//' decomposition and the narrative-restriction check require. Because the
//' instrument typically spans a shorter sample than the VAR, the completion
//' renormalises \code{b1} in the Cholesky basis so that \code{B \%*\% t(B)}
//' equals the full-sample \code{sigma} exactly and FEVD shares still sum to one.
//'
//' @references
//' Mertens, K., & Ravn, M. O. (2013). The dynamic effects of personal and
//' corporate income tax changes in the United States. \emph{American Economic
//' Review}, 103(4), 1212--1247.
//'
//' Gertler, M., & Karadi, P. (2015). Monetary policy surprises, credit costs,
//' and economic activity. \emph{AEJ: Macroeconomics}, 7(1), 44--76.
//'
//' @seealso \code{\link{fSR_cpp}}, \code{\link{fVAR}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fRecoverBIV_cpp(const arma::mat& resid_sub,
                           const arma::mat& Z_sub,
                           const arma::mat& sigma,
                           int ntotcoeff) {

    IVColumnResult r = fIVColumn_cpp(resid_sub, Z_sub, sigma, ntotcoeff);

    return Rcpp::List::create(Rcpp::Named("b1")       = r.b1,
                              Rcpp::Named("B")        = r.B,
                              Rcpp::Named("sigma_b")  = r.sigma_b,
                              Rcpp::Named("fs_beta")  = r.fs_beta,
                              Rcpp::Named("fs_F")     = r.fs_F,
                              Rcpp::Named("fs_r2")    = r.fs_r2,
                              Rcpp::Named("shock_sd") = r.shock_sd,
                              Rcpp::Named("n_iv")     = r.n_iv);
}
