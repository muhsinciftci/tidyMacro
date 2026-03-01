// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// CRITICAL: Define this BEFORE including RcppArmadillo to prevent conflicts
#define ARMA_DONT_USE_OPENMP

#include "fbootstrapIV_mbb.h"
#include "fOLS.h"
#include "fVAR.h"
#include "fgenerateVARdata.h"
#include "fmbb_var.h"
#include "fwoldIRF.h"
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// Internal C++ function (called from other C++ code)
// Calls _cpp versions of internal functions - no list overhead!
BootstrapIVMBBResult
fbootstrapIV_mbb_cpp(const arma::mat &y, const VARResult &var_result,
                     const arma::mat &Z, int nboot, int blocksize,
                     const arma::ivec &adjustZ, const arma::ivec &adjustu,
                     int policyvar, int horizon, double prc,
                     Rcpp::Nullable<arma::mat> exog, int n_threads) {

  // Extract VAR components from struct (much faster than from Rcpp::List)
  const arma::mat &beta = var_result.beta;
  const arma::mat &residuals = var_result.residuals;
  const int p = var_result.p;
  const int c = var_result.c;
  const int n_exog = var_result.n_exog;

  // If original VAR had exogenous variables, they must be provided
  if (n_exog > 0 && exog.isNull()) {
    Rcpp::stop("Original VAR estimation used exogenous variables. You must "
               "provide the 'exog' parameter.");
  }

  // If original VAR had no exogenous variables, they should not be provided
  if (n_exog == 0 && exog.isNotNull()) {
    Rcpp::stop("Original VAR estimation did not use exogenous variables. Do "
               "not provide the 'exog' parameter.");
  }

  const int N = y.n_cols; // Implicit conversion from arma::uword to int
  const int H = horizon + 1;

  // Convert 1-based policyvar to 0-based
  const int policyvar_idx = policyvar - 1;

  // Pre-allocate storage for bootstrapped IRFs
  arma::cube ivirf_boot(N, H, nboot, arma::fill::none);

  // Balance IV and residuals to overlap (convert 1-based to 0-based indices)
  const int Z_start = adjustZ(0) - 1;
  const int Z_end = adjustZ(1) - 1;
  const int u_start = adjustu(0) - 1;
  const int u_end = adjustu(1) - 1;

  const arma::mat Zstar = Z.rows(Z_start, Z_end);
  const arma::mat Epsstar = residuals.rows(u_start, u_end);

// Determine and set number of threads for OpenMP
int actual_threads = 1;  // Default to 1 if no OpenMP

#ifdef _OPENMP
// If n_threads is 0, use all available cores minus 1
if (n_threads == 0) {
int max_threads = omp_get_max_threads();
actual_threads = std::max(1, max_threads - 1);  // Leave one core free
} else {
actual_threads = n_threads;
}
omp_set_num_threads(actual_threads);
Rprintf("Using %d thread(s) for parallel bootstrap computation...\n", actual_threads);
#else
Rprintf("OpenMP not available. Running in single-threaded mode.\n");
actual_threads = 1;
#endif

// Bootstrap loop - parallelized with OpenMP
// Each iteration is independent, making this perfectly parallelizable
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int b = 0; b < nboot; ++b) {
    // Each thread needs its own workspace
    arma::mat res_boot(residuals.n_rows, residuals.n_cols, arma::fill::none);

    // Generate new residuals and new instruments using direct matrix overload
    // No Rcpp::wrap overhead - thread-safe and faster!
    MBBVARResult mbb_result = fmbb_var_cpp(Epsstar, p, blocksize, Zstar);
    const arma::mat &eps_boot = mbb_result.eps_boot; // Use reference - no copy!
    const arma::mat &Znew = mbb_result.M_boot;       // Use reference - no copy!

    // Pad the data with the original residuals
    if (u_start > 0) {
      res_boot.rows(0, u_start - 1) = residuals.rows(0, u_start - 1);
    }
    res_boot.rows(u_start, u_end) = eps_boot;

    if (u_end < static_cast<int>(residuals.n_rows) - 1) {
      res_boot.rows(u_end + 1, residuals.n_rows - 1) =
          residuals.rows(u_end + 1, residuals.n_rows - 1);
    }

    // Bootstrap the VAR
    arma::mat varboot = fgenerateVARdata(y, p, c, beta, res_boot);

    // Estimate VAR on bootstrap sample using _cpp version
    VARResult var_loop;
    if (exog.isNotNull()) {
      var_loop = fVAR_cpp(varboot, p, c, exog);
    } else {
      var_loop = fVAR_cpp(varboot, p, c, R_NilValue);
    }

    const arma::mat &residuals_loop =
        var_loop.residuals; // Use reference - no copy!

    // Extract policy variable residuals
    arma::vec u_p_loop_final =
        residuals_loop.col(policyvar_idx).rows(u_start, u_end);

    // Extract non-policy residuals
    arma::mat u_q_final = residuals_loop.rows(u_start, u_end);
    u_q_final.shed_col(policyvar_idx);

    // Compute Wold IRFs using _cpp version - no list overhead!
    WoldIRFResult wold_result = fwoldIRF_cpp(var_loop, horizon);
    const arma::cube &wold_loop =
        wold_result.irfwold; // Use reference - no copy!

    // First-stage regression: regress policy residual on instrument using _cpp
    // version
    OLSResult ols1_result = fOLS_cpp(u_p_loop_final, Znew, 0);
    const arma::mat &uhat =
        ols1_result.fitted_partial; // Use reference - no copy!

    // Second-stage: compute impact on non-policy variables
    arma::mat sq_sp = arma::solve(uhat.t() * uhat, uhat.t() * u_q_final);

    // Build structural impact vector (normalized sp = 1)
    arma::vec s(N, arma::fill::zeros);
    int counter = 0;
    for (int i = 0; i < N; ++i) {
      if (i == policyvar_idx) {
        s(i) = 1.0;
      } else {
        s(i) = sq_sp(counter);
        counter++;
      }
    }

    // Compute structural IRFs: IRF_h = Wold_h * s
    for (int h = 0; h < H; ++h) {
      ivirf_boot.slice(b).col(h) = wold_loop.slice(h) * s;
    }
  }

  // Compute percentiles
  const double up_pct = 50.0 + prc * 0.5;
  const double low_pct = 50.0 - prc * 0.5;

  // Pre-allocate output matrices
  arma::mat upper(N, H, arma::fill::none);
  arma::mat lower(N, H, arma::fill::none);
  arma::mat medianirf(N, H, arma::fill::none);
  arma::mat meanirf(N, H, arma::fill::none);

// Compute percentiles for each IRF element - also parallelized
// Each (i,h) combination is independent
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
  for (int i = 0; i < N; ++i) {
    for (int h = 0; h < H; ++h) {
      // Each thread needs its own boot_vals vector
      arma::vec boot_vals(nboot, arma::fill::none);
      // Extract bootstrap values for this element
      for (int b = 0; b < nboot; ++b) {
        boot_vals(b) = ivirf_boot(i, h, b);
      }

      arma::vec sorted = arma::sort(boot_vals);

      // Linear interpolation for percentiles
      const double idx_up = (up_pct / 100.0) * (nboot - 1);
      const double idx_low = (low_pct / 100.0) * (nboot - 1);
      const double idx_med = 0.5 * (nboot - 1);

      const int floor_up = static_cast<int>(idx_up);
      const int floor_low = static_cast<int>(idx_low);
      const int floor_med = static_cast<int>(idx_med);

      const double frac_up = idx_up - floor_up;
      const double frac_low = idx_low - floor_low;
      const double frac_med = idx_med - floor_med;

      // Compute percentiles with bounds checking
      if (floor_up < nboot - 1) {
        upper(i, h) =
            sorted(floor_up) * (1.0 - frac_up) + sorted(floor_up + 1) * frac_up;
      } else {
        upper(i, h) = sorted(nboot - 1);
      }

      if (floor_low < nboot - 1) {
        lower(i, h) = sorted(floor_low) * (1.0 - frac_low) +
                      sorted(floor_low + 1) * frac_low;
      } else {
        lower(i, h) = sorted(nboot - 1);
      }

      if (floor_med < nboot - 1) {
        medianirf(i, h) = sorted(floor_med) * (1.0 - frac_med) +
                          sorted(floor_med + 1) * frac_med;
      } else {
        medianirf(i, h) = sorted(nboot - 1);
      }

      // Mean
      meanirf(i, h) = arma::mean(boot_vals);
    }
  }

  // Return results as struct
  BootstrapIVMBBResult result;
  result.upper = upper;
  result.lower = lower;
  result.meanirf = meanirf;
  result.medianirf = medianirf;
  return result;
}

//' Bootstrap IV Impulse Response Functions using Moving Block Bootstrap
//'
//' Computes bootstrap confidence bands for instrumental variable identified
//' impulse response functions using the moving block bootstrap to preserve
//' temporal dependence in residuals and instruments.
//'
//' @param y T x N matrix of endogenous variables
//' @param var_result List output from fVAR containing beta, residuals, p, c
//' @param Z (T-p) x K matrix of instrumental variables
//' @param nboot Integer number of bootstrap replications
//' @param blocksize Integer block size for moving block bootstrap
//' @param adjustZ Integer vector of length 2: [start, end] indices for Z
// alignment ' @param adjustu Integer vector of length 2: [start, end] indices
// for residuals alignment ' @param policyvar Integer index (1-based) of the
// policy variable ' @param horizon Integer maximum impulse response horizon '
//@param prc Double percentile for confidence bands (e.g., 68 for 68% CI) '
//@param exog Optional matrix of exogenous variables (T x M). Default is NULL.
//' @param n_threads Integer number of threads for parallel computation.
//'   Default is 0 (uses all available cores). Set to 1 for single-threaded
// execution. '   If OpenMP is not available, automatically falls back to
// single-threaded.
//'
//' @return A list containing:
//'   \itemize{
//'     \item upper: N x (horizon+1) matrix of upper confidence bands
//'     \item lower: N x (horizon+1) matrix of lower confidence bands
//'     \item meanirf: N x (horizon+1) matrix of mean impulse responses
//'     \item medianirf: N x (horizon+1) matrix of median impulse responses
//'   }
//'
//' @details
//' This function implements the moving block bootstrap for IV-identified SVARs.
//' The first stage regresses the policy variable residual on the instrument(s),
//' and the second stage recovers the structural impact matrix. The
// normalization ' sets the policy variable shock to have unit impact on itself.
//'
//' The function uses OpenMP for parallel computation when available,
// significantly ' speeding up bootstrap iterations. Each bootstrap replication
// and percentile ' computation is independent and can be parallelized. If
// OpenMP is not available, ' the function automatically falls back to
// single-threaded execution.
//'
//' @examples
//' \dontrun{
//' var_result <- fVAR(y, p = 2, c = 1)
//' result <- fbootstrapIV_mbb(y, var_result, Z,
//'                             nboot = 1000, blocksize = 10,
//'                             adjustZ = c(1, 100), adjustu = c(1, 100),
//'                             policyvar = 1, horizon = 20, prc = 68)
//'
//' # Use 4 threads for parallel computation
//' result <- fbootstrapIV_mbb(y, var_result, Z,
//'                             nboot = 1000, blocksize = 10,
//'                             adjustZ = c(1, 100), adjustu = c(1, 100),
//'                             policyvar = 1, horizon = 20, prc = 68,
//'                             n_threads = 4)
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List fbootstrapIV_mbb(const arma::mat &y, const Rcpp::List &var_result,
                            const arma::mat &Z, int nboot, int blocksize,
                            const arma::ivec &adjustZ,
                            const arma::ivec &adjustu, int policyvar,
                            int horizon, double prc,
                            Rcpp::Nullable<arma::mat> exog = R_NilValue,
                            int n_threads = 0) {

  // Extract elements from R list and construct VARResult struct
  VARResult var_result_struct;
  var_result_struct.beta = Rcpp::as<arma::mat>(var_result["beta"]);
  var_result_struct.residuals = Rcpp::as<arma::mat>(var_result["residuals"]);
  var_result_struct.sigma_full = Rcpp::as<arma::mat>(var_result["sigma_full"]);
  var_result_struct.p = Rcpp::as<int>(var_result["p"]);
  var_result_struct.c = Rcpp::as<int>(var_result["c"]);

  // Check if n_exog exists in the list
  if (var_result.containsElementNamed("n_exog")) {
    var_result_struct.n_exog = Rcpp::as<int>(var_result["n_exog"]);
  } else {
    var_result_struct.n_exog = 0;
  }

  // Call the C++ function with the struct
  BootstrapIVMBBResult result =
      fbootstrapIV_mbb_cpp(y, var_result_struct, Z, nboot, blocksize, adjustZ,
                           adjustu, policyvar, horizon, prc, exog, n_threads);

  // Return results as a list for R
  return Rcpp::List::create(Rcpp::Named("upper") = result.upper,
                            Rcpp::Named("lower") = result.lower,
                            Rcpp::Named("meanirf") = result.meanirf,
                            Rcpp::Named("medianirf") = result.medianirf);
}
