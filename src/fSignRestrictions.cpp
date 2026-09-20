// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_DONT_USE_OPENMP

#include "fSignRestrictions.h"
#include "fGenerateQ.h"
#include "fWoldIRF.h"
#include "fVAR.h"
#include <RcppArmadillo.h>
#include <algorithm>

void fSignRotationPrep_cpp(const arma::mat& sigma,
                           const arma::mat& Bfix,
                           tidymacro::RNG&  rng,
                           SignRotScratch&  scr) {

    const arma::uword N       = sigma.n_rows;
    const arma::uword n_fixed = Bfix.n_cols;

    if (!arma::chol(scr.C, sigma, "lower")) {
        Rcpp::stop("The reduced-form covariance matrix is not positive definite.");
    }

    scr.startingMat.set_size(N, N);

    if (n_fixed == 0u) {
        // No pre-determined column: the rotation basis is the Cholesky factor.
        scr.startingMat = scr.C;
        return;
    }

    // Express the fixed columns in the Cholesky basis and complete them to an
    // orthonormal basis by Gram-Schmidt on standard normal directions.  The
    // completion is drawn once per posterior draw rather than once per rotation:
    // the Haar rotation applied to the free block afterwards makes any
    // orthonormal completion of the same subspace distributionally equivalent.
    arma::mat q(N, N, arma::fill::none);
    for (arma::uword j = 0; j < n_fixed; ++j) {
        arma::vec qj = arma::solve(arma::trimatl(scr.C), Bfix.col(j));
        // C \ Bfix is only approximately unit-norm when Bfix was calibrated on a
        // different covariance (the instrument subsample); normalise defensively.
        for (arma::uword l = 0; l < j; ++l) qj -= arma::dot(q.col(l), qj) * q.col(l);
        const double nrm = arma::norm(qj, 2);
        if (nrm < 1e-12) Rcpp::stop("Pre-determined impact columns are collinear.");
        q.col(j) = qj / nrm;
    }
    for (arma::uword j = n_fixed; j < N; ++j) {
        arma::vec r(N, arma::fill::none);
        double nrm = 0.0;
        int guard = 0;
        do {
            for (arma::uword i = 0; i < N; ++i) r(i) = rng.norm();
            for (arma::uword l = 0; l < j; ++l) r -= arma::dot(q.col(l), r) * q.col(l);
            nrm = arma::norm(r, 2);
        } while (nrm < 1e-10 && ++guard < 100);
        q.col(j) = r / nrm;
    }

    scr.startingMat = scr.C * q;
}

bool fSignRotation_cpp(const arma::mat&  SIGN,
                       const arma::cube& wold,
                       int               sr_hor,
                       int               sr_rot,
                       int               n_fixed,
                       tidymacro::RNG&   rng,
                       SignRotScratch&   scr,
                       arma::mat&        B_out,
                       int&              n_tried) {

    const arma::uword N   = scr.startingMat.n_rows;
    const arma::uword ws  = static_cast<arma::uword>(n_fixed); // first free column
    const arma::uword m   = N - ws;                            // free columns
    const arma::uword ds  = SIGN.n_cols;                       // shocks to match
    const arma::uword H   = static_cast<arma::uword>(sr_hor);

    if (SIGN.n_rows != N)  Rcpp::stop("SIGN must have one row per variable.");
    if (ws + ds > N)       Rcpp::stop("SIGN has more shock columns than free columns of B.");

    scr.used.assign(N, 0);
    scr.order.resize(N);

    if (scr.termaa.n_rows != N) scr.termaa.set_size(N, N);
    if (H > 1 && (scr.irfchk.n_rows != N || scr.irfchk.n_cols != m ||
                  scr.irfchk.n_slices != H)) {
        scr.irfchk.set_size(N, m, H);
    }

    for (int attempt = 1; attempt <= sr_rot; ++attempt) {

        // --- rotate the free block -------------------------------------
        fGenerateQ_inplace(scr.Qs, scr.Rs, scr.Gs, m, rng);
        scr.rotated = scr.startingMat.cols(ws, N - 1) * scr.Qs;

        if (ws > 0) scr.termaa.cols(0, ws - 1) = scr.startingMat.cols(0, ws - 1);
        scr.termaa.cols(ws, N - 1) = scr.rotated;

        // --- responses used by the sign check --------------------------
        // sr_hor == 1 checks the impact matrix itself, so no IRF is needed.
        if (H > 1) {
            scr.irfchk.slice(0) = scr.rotated;   // wold.slice(0) is the identity
            for (arma::uword h = 1; h < H; ++h) {
                scr.irfchk.slice(h) = wold.slice(h) * scr.rotated;
            }
        }

        // --- greedily match each shock to a free column ----------------
        std::fill(scr.used.begin(), scr.used.end(), 0);
        for (arma::uword i = 0; i < N; ++i) scr.order[i] = i;

        arma::uword matched = 0;
        for (arma::uword ii = 0; ii < ds; ++ii) {
            for (arma::uword jj = ws; jj < N; ++jj) {
                if (scr.used[jj]) continue;

                // Mirrors the toolbox rule: only a strictly wrong-signed
                // response rejects, so exact zeros and unrestricted rows pass.
                bool ok_pos = true, ok_neg = true;
                for (arma::uword h = 0; h < H && (ok_pos || ok_neg); ++h) {
                    const double* col = (H > 1)
                        ? scr.irfchk.slice(h).colptr(jj - ws)
                        : scr.termaa.colptr(jj);
                    for (arma::uword i = 0; i < N; ++i) {
                        const double s = SIGN(i, ii);
                        if (s == 0.0) continue;
                        const double prod = s * col[i];
                        if (prod < 0.0) ok_pos = false;
                        else if (prod > 0.0) ok_neg = false;
                        if (!ok_pos && !ok_neg) break;
                    }
                }

                if (ok_pos) {
                    scr.used[jj] = 1;
                    scr.order[ws + ii] = jj;
                    ++matched;
                    break;
                }
                if (ok_neg) {
                    scr.used[jj] = 1;
                    scr.termaa.col(jj) *= -1.0;
                    scr.order[ws + ii] = jj;
                    ++matched;
                    break;
                }
            }
        }

        if (matched == ds) {
            B_out.set_size(N, N);
            for (arma::uword j = 0; j < N; ++j) B_out.col(j) = scr.termaa.col(scr.order[j]);
            n_tried = attempt;
            return true;
        }
    }

    n_tried = sr_rot;
    B_out.reset();
    return false;
}

//' Find a Rotation Satisfying Sign Restrictions
//'
//' Draws random orthonormal rotations of the reduced-form covariance matrix
//' until the implied structural impact matrix satisfies a sign-restriction
//' pattern, optionally holding pre-identified columns fixed.
//'
//' @param sigma N x N reduced-form residual covariance matrix.
//' @param SIGN N x ds sign-restriction matrix: \code{+1} the response must be
//'   non-negative, \code{-1} non-positive, \code{0} unrestricted. \code{ds} is
//'   the number of shocks to be matched and must equal \code{N} minus the
//'   number of columns supplied in \code{Bfix}.
//' @param sr_hor Integer. Restrictions are imposed at horizons
//'   \code{0, ..., sr_hor - 1}. \code{sr_hor = 1} (default) restricts the
//'   impact matrix only, in which case \code{beta} is not used.
//' @param sr_rot Integer maximum number of rotations to attempt (default 500).
//' @param Bfix Optional N x q matrix of impact columns identified by another
//'   scheme and held fixed (default NULL). Used for \code{sign+iv}.
//' @param beta Optional VAR coefficient matrix, required when
//'   \code{sr_hor > 1} to build the Wold multipliers (default NULL).
//' @param p Integer lag order, required when \code{sr_hor > 1} (default 1).
//' @param c Integer intercept indicator, required when \code{sr_hor > 1}
//'   (default 1).
//' @param seed Integer seed for the rotation draws (default 42).
//'
//' @return A list with \code{B} (N x N impact matrix, or a 0 x 0 matrix when no
//'   admissible rotation was found), \code{n_tried} (rotations attempted) and
//'   \code{found} (logical).
//'
//' @details
//' Shocks are matched to columns greedily: shock \code{ii} takes the first
//' unmatched column whose responses satisfy \code{SIGN[, ii]}, with the sign of
//' that column flipped if the restrictions hold in reverse. This reproduces the
//' matching rule of \code{SignRestrictions.m} in the VAR Toolbox.
//'
//' @references
//' Rubio-Ramirez, J. F., Waggoner, D. F., & Zha, T. (2010). Structural vector
//' autoregressions. \emph{Review of Economic Studies}, 77(2), 665--696.
//'
//' @seealso \code{\link{fSR_cpp}}, \code{\link{fGenerateQ}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fSignRestrictions_cpp(const arma::mat& sigma,
                                 const arma::mat& SIGN,
                                 int sr_hor = 1,
                                 int sr_rot = 500,
                                 Rcpp::Nullable<arma::mat> Bfix = R_NilValue,
                                 Rcpp::Nullable<arma::mat> beta = R_NilValue,
                                 int p = 1,
                                 int c = 1,
                                 int seed = 42) {

    if (sr_hor < 1)  Rcpp::stop("sr_hor must be at least 1.");
    if (sr_rot < 1)  Rcpp::stop("sr_rot must be at least 1.");

    arma::mat Bfix_mat;
    if (Bfix.isNotNull()) Bfix_mat = Rcpp::as<arma::mat>(Bfix);

    arma::cube wold;
    if (sr_hor > 1) {
        if (beta.isNull()) {
            Rcpp::stop("beta is required when sr_hor > 1.");
        }
        VARResult vr;
        vr.beta   = Rcpp::as<arma::mat>(beta);
        vr.sigma  = sigma;
        vr.p      = p;
        vr.c      = c;
        vr.n_exog = 0;
        wold = fWoldIRF_cpp(vr, sr_hor - 1).irfwold;
    }

    tidymacro::RNG rng(static_cast<std::uint64_t>(seed));
    SignRotScratch scr;
    fSignRotationPrep_cpp(sigma, Bfix_mat, rng, scr);

    arma::mat B;
    int n_tried = 0;
    const bool found = fSignRotation_cpp(SIGN, wold, sr_hor, sr_rot,
                                         static_cast<int>(Bfix_mat.n_cols),
                                         rng, scr, B, n_tried);

    return Rcpp::List::create(Rcpp::Named("B")       = B,
                              Rcpp::Named("n_tried") = n_tried,
                              Rcpp::Named("found")   = found);
}
