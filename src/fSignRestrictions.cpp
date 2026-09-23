// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_DONT_USE_OPENMP

#include "fSignRestrictions.h"
#include "fGenerateQ.h"
#include "fWoldIRF.h"
#include "fVAR.h"
#include <RcppArmadillo.h>
#include <algorithm>
#include <stdexcept>

namespace {
// A conservative Gordan/Farkas certificate for the small free subspaces in
// the two-IV replication. If m+1 signed rows of full rank have a strictly
// positive null combination, A*q >= 0 implies q=0: no unit rotation can pass.
// Only well-conditioned certificates, far from zero weights, are used.
bool infeasible_impact_cpp(const arma::mat& checks, arma::uword nrows) {
    const arma::uword m = checks.n_rows;
    if (m == 0 || m > 3 || nrows < m + 1) return false;
    arma::mat rows = checks.cols(0, nrows - 1);
    for (arma::uword j = 0; j < nrows; ++j) {
        const double norm = arma::norm(rows.col(j), 2);
        if (!(norm > 0.0) || !std::isfinite(norm)) return false;
        rows.col(j) /= norm;
    }
    std::vector<arma::uword> subset(m + 1);
    for (arma::uword i = 0; i <= m; ++i) subset[i] = i;
    arma::mat system(m + 1, m + 1, arma::fill::ones);
    arma::vec target(m + 1, arma::fill::zeros), weights;
    target(m) = 1.0;
    for (int budget = 0; budget < 64; ++budget) {
        for (arma::uword j = 0; j <= m; ++j)
            system.submat(0, j, m - 1, j) = rows.col(subset[j]);
        if (arma::rcond(system) > 1e-6 &&
            arma::solve(weights, system, target,
                        arma::solve_opts::fast + arma::solve_opts::no_approx) &&
            weights.min() > 1e-8 &&
            arma::norm(system * weights - target, "inf") < 1e-12) return true;
        int i = static_cast<int>(m);
        while (i >= 0 && subset[i] == nrows - (m + 1) + i) --i;
        if (i < 0) break;
        ++subset[i];
        for (arma::uword j = i + 1; j <= m; ++j) subset[j] = subset[j - 1] + 1;
    }
    return false; // Absence of a certificate says nothing about feasibility.
}
} // namespace

void fSignRotationPrep_cpp(const arma::mat& sigma,
                           const arma::mat& Bfix,
                           tidymacro::RNG&  rng,
                           SignRotScratch&  scr) {

    const arma::uword N       = sigma.n_rows;
    const arma::uword n_fixed = Bfix.n_cols;
    scr.restrictions_ready = false;
    scr.infeasible = false;
    if (N == 0 || sigma.n_cols != N || !sigma.is_finite() ||
        n_fixed >= N || (n_fixed > 0 && (Bfix.n_rows != N || !Bfix.is_finite())))
        throw std::invalid_argument("Invalid covariance or fixed-column dimensions/values.");

    if (!arma::chol(scr.C, sigma, "lower")) {
        throw std::runtime_error("The reduced-form covariance matrix is not positive definite.");
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
        if (!std::isfinite(nrm) || nrm < 1e-12)
            throw std::runtime_error("Pre-determined impact columns are collinear.");
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
        if (!std::isfinite(nrm) || nrm < 1e-10)
            throw std::runtime_error("Failed to complete the fixed-column basis.");
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

    if (!scr.restrictions_ready) {
        if (SIGN.n_rows != N || ws + ds > N || sr_hor < 1 || sr_rot < 1 ||
            !SIGN.is_finite() ||
            arma::any(arma::vectorise((SIGN != -1.0) % (SIGN != 0.0) % (SIGN != 1.0))))
            throw std::invalid_argument("Invalid sign restrictions.");
        if (H > 1 && (wold.n_rows != N || wold.n_cols != N || wold.n_slices < H))
            throw std::invalid_argument("Too few Wold horizons for the restrictions.");
        scr.used.resize(N);
        scr.order.resize(N);
        scr.orientation.resize(N);
        scr.n_restr.resize(ds);
        scr.signed_basis.resize(ds);
        // Factor Psi_h * startingMat out of the rotation loop. Only constrained
        // rows are needed; each candidate sign is a short dot product with Q.
        arma::cube basis(N, m, H, arma::fill::none);
        basis.slice(0) = scr.startingMat.cols(ws, N - 1);
        for (arma::uword h = 1; h < H; ++h)
            basis.slice(h) = wold.slice(h) * basis.slice(0);
        for (arma::uword ii = 0; ii < ds; ++ii) {
            const arma::uvec rows = arma::find(SIGN.col(ii) != 0.0);
            scr.n_restr[ii] = rows.n_elem;
            arma::mat& checks = scr.signed_basis[ii];
            checks.set_size(m, rows.n_elem * H);
            for (arma::uword h = 0; h < H; ++h)
                for (arma::uword r = 0; r < rows.n_elem; ++r)
                    checks.col(h * rows.n_elem + r) =
                        SIGN(rows(r), ii) * basis.slice(h).row(rows(r)).t();
        }
        if (scr.impact_match_first && ws > 0)
            for (arma::uword ii = 0; ii < ds && !scr.infeasible; ++ii)
                scr.infeasible = infeasible_impact_cpp(scr.signed_basis[ii], scr.n_restr[ii]);
        scr.restrictions_ready = true;
    }

    if (scr.infeasible) {
        n_tried = 0;
        B_out.reset();
        return false;
    }
    for (int attempt = 1; attempt <= sr_rot; ++attempt) {
        fGenerateQ_inplace(scr.Qs, scr.Rs, scr.Gs, m, rng);
        std::fill(scr.used.begin(), scr.used.end(), 0);
        std::fill(scr.orientation.begin(), scr.orientation.end(), 1.0);
        for (arma::uword i = 0; i < N; ++i) scr.order[i] = i;

        bool matched = true;
        for (arma::uword ii = 0; ii < ds; ++ii) {
            const arma::mat& checks = scr.signed_basis[ii];
            const arma::uword count = scr.n_restr[ii] * (scr.impact_match_first ? 1 : H);
            bool found = false;
            for (arma::uword jj = ws; jj < N; ++jj) {
                if (scr.used[jj]) continue;
                bool ok_pos = true, ok_neg = true;
                const double* q = scr.Qs.colptr(jj - ws);
                for (arma::uword r = 0; r < count; ++r) {
                    const double* b = checks.colptr(r);
                    double value = 0.0;
                    for (arma::uword l = 0; l < m; ++l) value += b[l] * q[l];
                    if (!std::isfinite(value)) { ok_pos = ok_neg = false; break; }
                    if (value < 0.0) ok_pos = false;
                    else if (value > 0.0) ok_neg = false;
                    if (!ok_pos && !ok_neg) break;
                }
                if (ok_pos || ok_neg) {
                    scr.used[jj] = 1;
                    scr.order[ws + ii] = jj;
                    scr.orientation[ws + ii] = ok_pos ? 1.0 : -1.0;
                    found = true;
                    break;
                }
            }
            if (!found) { matched = false; break; }
        }
        if (!matched) continue;

        // The paper's collector fixes the impact assignment before checking
        // later horizons. It does not try another column after a later failure.
        if (scr.impact_match_first && H > 1) {
            for (arma::uword ii = 0; ii < ds && matched; ++ii) {
                const arma::mat& checks = scr.signed_basis[ii];
                const double* q = scr.Qs.colptr(scr.order[ws + ii] - ws);
                for (arma::uword r = scr.n_restr[ii]; r < checks.n_cols; ++r) {
                    const double* b = checks.colptr(r);
                    double value = 0.0;
                    for (arma::uword l = 0; l < m; ++l) value += b[l] * q[l];
                    if (!(scr.orientation[ws + ii] * value > 0.0)) {
                        matched = false; break;
                    }
                }
            }
            if (!matched) continue;
        }

        // Complete an optional partial SIGN without duplicating matched columns.
        arma::uword next = ws + ds;
        for (arma::uword jj = ws; jj < N; ++jj)
            if (!scr.used[jj]) scr.order[next++] = jj;
        scr.rotated = scr.startingMat.cols(ws, N - 1) * scr.Qs;
        B_out.set_size(N, N);
        if (ws > 0) B_out.cols(0, ws - 1) = scr.startingMat.cols(0, ws - 1);
        for (arma::uword j = ws; j < N; ++j)
            B_out.col(j) = scr.orientation[j] * scr.rotated.col(scr.order[j] - ws);
        n_tried = attempt;
        return true;
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

    if (sigma.n_rows == 0 || sigma.n_rows != sigma.n_cols ||
        SIGN.n_rows != sigma.n_rows || !sigma.is_finite() || !SIGN.is_finite())
        Rcpp::stop("Invalid covariance or sign matrix dimensions/values.");
    if (Bfix_mat.n_cols >= sigma.n_rows ||
        (Bfix_mat.n_cols > 0 && (Bfix_mat.n_rows != sigma.n_rows || !Bfix_mat.is_finite())))
        Rcpp::stop("Invalid fixed-column dimensions/values.");
    if (SIGN.n_cols + Bfix_mat.n_cols > sigma.n_rows)
        Rcpp::stop("Too many restricted and fixed columns.");
    if (p < 1 || (c != 0 && c != 1)) Rcpp::stop("Invalid lag order or intercept.");
    arma::cube wold;
    if (sr_hor > 1) {
        if (beta.isNull()) {
            Rcpp::stop("beta is required when sr_hor > 1.");
        }
        VARResult vr;
        vr.beta   = Rcpp::as<arma::mat>(beta);
        if (vr.beta.n_cols != sigma.n_rows ||
            vr.beta.n_rows < sigma.n_rows * p + c || !vr.beta.is_finite())
            Rcpp::stop("Invalid beta dimensions/values for the requested VAR.");
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

//' Collect All Rotations Satisfying Sign Restrictions
//'
//' Runs a fixed number of Haar rotations and retains every impact matrix that
//' satisfies the requested sign restrictions. This is the storage convention
//' used by Cesa-Bianchi and Sokol's \code{signRestrictions.m}; most callers
//' should use \code{\link{fSignRestrictions_cpp}}, which stops at the first
//' admissible rotation.
//'
//' @inheritParams fSignRestrictions_cpp
//'
//' @return A list with \code{Ball}, an N x N x \code{n_found} array of
//'   admissible impact matrices, \code{n_found}, \code{n_tried}, and
//'   \code{found}.
//'
//' @seealso \code{\link{fSignRestrictions_cpp}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fSignRestrictionsAll_cpp(const arma::mat& sigma,
                                    const arma::mat& SIGN,
                                    int sr_hor = 1,
                                    int sr_rot = 500,
                                    Rcpp::Nullable<arma::mat> Bfix = R_NilValue,
                                    Rcpp::Nullable<arma::mat> beta = R_NilValue,
                                    int p = 1,
                                    int c = 1,
                                    int seed = 42) {

    if (sr_hor < 1) Rcpp::stop("sr_hor must be at least 1.");
    if (sr_rot < 1) Rcpp::stop("sr_rot must be at least 1.");

    const arma::uword N = sigma.n_rows;
    if (sigma.n_cols != N) Rcpp::stop("sigma must be square.");
    if (SIGN.n_rows != N) Rcpp::stop("SIGN must have one row per variable.");

    arma::mat Bfix_mat;
    if (Bfix.isNotNull()) Bfix_mat = Rcpp::as<arma::mat>(Bfix);
    if (Bfix_mat.n_cols > 0 && Bfix_mat.n_rows != N) {
        Rcpp::stop("Bfix must have one row per variable.");
    }
    if (SIGN.n_cols + Bfix_mat.n_cols != N) {
        Rcpp::stop("SIGN must have N - ncol(Bfix) columns.");
    }

    if (sigma.n_rows == 0 || sigma.n_rows != sigma.n_cols ||
        SIGN.n_rows != sigma.n_rows || !sigma.is_finite() || !SIGN.is_finite())
        Rcpp::stop("Invalid covariance or sign matrix dimensions/values.");
    if (Bfix_mat.n_cols >= sigma.n_rows ||
        (Bfix_mat.n_cols > 0 && (Bfix_mat.n_rows != sigma.n_rows || !Bfix_mat.is_finite())))
        Rcpp::stop("Invalid fixed-column dimensions/values.");
    if (SIGN.n_cols + Bfix_mat.n_cols > sigma.n_rows)
        Rcpp::stop("Too many restricted and fixed columns.");
    if (p < 1 || (c != 0 && c != 1)) Rcpp::stop("Invalid lag order or intercept.");
    arma::cube wold;
    if (sr_hor > 1) {
        if (beta.isNull()) Rcpp::stop("beta is required when sr_hor > 1.");
        VARResult vr;
        vr.beta   = Rcpp::as<arma::mat>(beta);
        if (vr.beta.n_cols != sigma.n_rows ||
            vr.beta.n_rows < sigma.n_rows * p + c || !vr.beta.is_finite())
            Rcpp::stop("Invalid beta dimensions/values for the requested VAR.");
        vr.sigma  = sigma;
        vr.p      = p;
        vr.c      = c;
        vr.n_exog = 0;
        wold = fWoldIRF_cpp(vr, sr_hor - 1).irfwold;
    }

    tidymacro::RNG rng(static_cast<std::uint64_t>(seed));
    SignRotScratch scr;
    scr.impact_match_first = true;
    fSignRotationPrep_cpp(sigma, Bfix_mat, rng, scr);

    std::vector<double> accepted;
    const arma::uword block = N * N;
    accepted.reserve(block * static_cast<std::size_t>(std::min(sr_rot, 1024)));

    arma::mat B;
    for (int attempt = 0; attempt < sr_rot; ++attempt) {
        int n_tried = 0;
        const bool found = fSignRotation_cpp(
            SIGN, wold, sr_hor, 1, static_cast<int>(Bfix_mat.n_cols),
            rng, scr, B, n_tried);
        if (scr.infeasible) break;
        if (found) accepted.insert(accepted.end(), B.begin(), B.end());
        if (attempt % 4096 == 0) Rcpp::checkUserInterrupt();
    }

    const arma::uword n_found = accepted.size() / block;
    arma::cube Ball(accepted.data(), N, N, n_found, false);

    return Rcpp::List::create(
        Rcpp::Named("Ball")    = Ball,
        Rcpp::Named("n_found") = static_cast<int>(n_found),
        Rcpp::Named("n_tried") = scr.infeasible ? 0 : sr_rot,
        Rcpp::Named("n_ruled_out") = scr.infeasible ? sr_rot : 0,
        Rcpp::Named("infeasible") = scr.infeasible,
        Rcpp::Named("found")   = !accepted.empty());
}
