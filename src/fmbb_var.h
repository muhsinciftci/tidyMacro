#ifndef FMBB_VAR_H
#define FMBB_VAR_H

#include <RcppArmadillo.h>

// Struct to hold MBB VAR results
struct MBBVARResult {
  arma::mat eps_boot;  // (T-p) x N matrix of bootstrapped residuals (centered)
  arma::mat M_boot;    // (T-p) x K matrix of bootstrapped instruments (centered)
};

// Precomputed state: build once before the bootstrap loop, reuse across reps.
struct MBBSamplerState {
  arma::cube Blocks;           // (blockSize x N x nBlocks_avail)
  arma::cube MBlocks;          // (blockSize x K x nBlocks_avail), empty if !has_M
  arma::mat  centering_full;   // (final_size x N)
  arma::mat  Mcentering_full;  // (final_size x K), empty if !has_M
  int nBlock;
  int nBlocks_avail;
  int final_size;
  int actualBlockSize;
  int N;
  int k;
  bool has_M;
};

// Build the sampler state once from fixed data (eps, lags, BlockSize, M).
MBBSamplerState fmbb_precompute_cpp(const arma::mat& eps, int lags, int BlockSize,
                                    const arma::mat& M);
MBBSamplerState fmbb_precompute_cpp(const arma::mat& eps, int lags, int BlockSize);

// Draw one bootstrap sample from the precomputed state (call inside the loop).
MBBVARResult fmbb_sample_cpp(const MBBSamplerState& state);

// Internal C++ function - Nullable version (for use when M might be null)
MBBVARResult fmbb_var_cpp(const arma::mat& eps,
                          int lags,
                          int BlockSize,
                          Rcpp::Nullable<arma::mat> M);

// Internal C++ function - Direct matrix version (for use in loops, thread-safe)
MBBVARResult fmbb_var_cpp(const arma::mat& eps,
                          int lags,
                          int BlockSize,
                          const arma::mat& M);

// R wrapper function (for calling from R)
Rcpp::List fmbb_var(const arma::mat& eps,
                    int lags,
                    int BlockSize,
                    Rcpp::Nullable<arma::mat> M);

#endif // FMBB_VAR_H
