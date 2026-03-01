#include <RcppArmadillo.h>
#include "fcompanionMatrix.h"
#include "fremove_bias.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Internal C++ function (called from other C++ code)
// Calls fcompanionMatrix_cpp directly - no list overhead!
RemoveBiasResult fremove_bias_cpp(const arma::mat&  beta,
                                  int               c,
                                  int               p,
                                  const arma::cube& boot_beta) {

    const int Np_c = beta.n_rows;  // Implicit conversion from arma::uword to int
    const int N = beta.n_cols;     // Implicit conversion from arma::uword to int
    const int nboot = boot_beta.n_slices;  // Implicit conversion from arma::uword to int

    // ------------------------------------------------------------------ //
    //  Compute bootstrap mean across slices                               //
    // ------------------------------------------------------------------ //
    arma::mat boot_mean(Np_c, N, arma::fill::zeros);
    for (int b = 0; b < nboot; ++b) {
        boot_mean += boot_beta.slice(b);
    }
    boot_mean /= nboot;

    // bias = E[beta_boot]' - beta'   (both are Np_c x N here)
    arma::mat bias = boot_mean - beta.t();

    // ------------------------------------------------------------------ //
    //  Check stability of the original estimates                          //
    // ------------------------------------------------------------------ //
    CompanionMatrixResult comp0 = fcompanionMatrix_cpp(beta, c, p);
    arma::mat BigA0 = comp0.comp;
    arma::cx_vec ev0 = arma::eig_gen(BigA0);
    double max_ev = arma::max(arma::abs(ev0));

    arma::mat Beta;
    int corrections = 1;

    if (max_ev >= 1.0) {
        // Original is already explosive — return uncorrected
        RemoveBiasResult result;
        result.Beta = beta.t();
        result.corrections = corrections;
        return result;
    }

    // ------------------------------------------------------------------ //
    //  Apply full bias correction and verify stability                    //
    // ------------------------------------------------------------------ //
    Beta = beta.t() - bias;   // N x Np_c

    {
        CompanionMatrixResult comp1 = fcompanionMatrix_cpp(Beta.t(), c, p);
        arma::mat A1 = comp1.comp;
        arma::cx_vec ev1 = arma::eig_gen(A1);
        max_ev = arma::max(arma::abs(ev1));
    }

    // ------------------------------------------------------------------ //
    //  Iterative shrinkage if corrected system is explosive               //
    // ------------------------------------------------------------------ //
    double delta = 1.0;

    while (max_ev >= 1.0) {
        delta -= 0.01;
        corrections += 1;

        if (delta < 0.0 || corrections > 200) {
            // Give up: return the uncorrected estimates
            Beta = beta.t();
            break;
        }

        Beta = beta.t() - bias * delta;

        CompanionMatrixResult comp_loop = fcompanionMatrix_cpp(Beta.t(), c, p);
        arma::mat A_loop = comp_loop.comp;
        arma::cx_vec ev_loop = arma::eig_gen(A_loop);
        max_ev = arma::max(arma::abs(ev_loop));
    }

    // Return results as struct
    RemoveBiasResult result;
    result.Beta = Beta;
    result.corrections = corrections;
    return result;
}

//' Remove Small-Sample Bias from VAR Coefficient Estimates
//'
//' Computes the bootstrap mean of the VAR coefficients, subtracts the bias
//' from the point estimates, and iteratively shrinks the correction towards
//' zero if the bias-corrected system is explosive (all eigenvalues of the
//' companion matrix must be strictly less than one).
//'
//' @param beta (Np+c+M) x N matrix of OLS VAR coefficient estimates
//' @param c Integer intercept indicator (1 if intercept included, 0 otherwise)
//' @param p Integer VAR lag order
//' @param boot_beta (Np+c+M) x N x nboot cube of bootstrapped VAR coefficient
//'   estimates collected from a first-pass bootstrap
//'
//' @return A list with elements:
//'   Beta: N x (Np+c+M) bias-corrected coefficient matrix (transposed
//'         relative to the input so it matches the convention used by
//'         \code{fbootstrapCholCorrected})
//'   corrections: integer count of iterative shrinkage steps applied
//'                (1 means the first attempt was already stable)
//'
//' @details
//' The bias is estimated as
//'
//'   bias = mean(boot_beta, dim=3) - beta'
//'
//' If the original estimates are already explosive the function returns the
//' uncorrected estimates.  Otherwise the function applies the full bias
//' correction and checks stability.  If the corrected system is still
//' explosive (rare in practice) it reduces the correction factor by one
//' percent per iteration until stability is achieved or 200 attempts are
//' exhausted, at which point the uncorrected estimates are returned.
//'
//' This mirrors \code{remove_bias.m} from the MATLAB reference implementation.
//'
//' @seealso \code{\link{fbootstrapCholCorrected}}, \code{\link{fcompanionMatrix}}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fremove_bias(const arma::mat&  beta,
                        int               c,
                        int               p,
                        const arma::cube& boot_beta) {
    // Call the C++ function
    RemoveBiasResult result = fremove_bias_cpp(beta, c, p, boot_beta);
    
    // Return results as a list for R
    return Rcpp::List::create(
        Rcpp::Named("Beta")        = result.Beta,
        Rcpp::Named("corrections") = result.corrections
    );
}
