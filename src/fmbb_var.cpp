#include <RcppArmadillo.h>
#include "fmbb_var.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Internal C++ function (called from other C++ code)
MBBVARResult fmbb_var_cpp(const arma::mat& eps, 
                          int lags, 
                          int BlockSize, 
                          Rcpp::Nullable<arma::mat> M) {
         
  const int Tp = eps.n_rows;  // Implicit conversion from arma::uword to int
  const int N = eps.n_cols;   // Implicit conversion from arma::uword to int
  const int T = Tp + lags;
  
  // Compute block size if not provided or invalid
  int actualBlockSize = BlockSize;
  if (actualBlockSize <= 0) {
    actualBlockSize = static_cast<int>(std::round(5.03 * std::pow(Tp, 0.25)));
  }
  
  // Number of blocks needed
  const int nBlock = static_cast<int>(std::ceil(static_cast<double>(T - lags) / actualBlockSize));
  const int nBlocks_avail = T - lags - actualBlockSize + 1;
  const int final_size = T - lags;
  
  // Check if M is provided
  const bool has_M = M.isNotNull();
  arma::mat M_mat;
  int k = 0;
  if (has_M) {
    M_mat = Rcpp::as<arma::mat>(M);
    k = M_mat.n_cols;  // Implicit conversion from arma::uword to int
  }
  
  // Create blocks for residuals
  arma::cube Blocks(actualBlockSize, N, nBlocks_avail, arma::fill::none);
  for (int j = 0; j < nBlocks_avail; ++j) {
    Blocks.slice(j) = eps.rows(j, j + actualBlockSize - 1);
  }
  
  // Create blocks for instruments if provided
  arma::cube MBlocks;
  if (has_M) {
    MBlocks.set_size(actualBlockSize, k, nBlocks_avail);
    for (int j = 0; j < nBlocks_avail; ++j) {
      MBlocks.slice(j) = M_mat.rows(j, j + actualBlockSize - 1);
    }
  }
  
  // Pre-compute centering for residuals
  arma::mat centering(actualBlockSize, N, arma::fill::none);
  const int center_end = T - lags - actualBlockSize;
  for (int j = 0; j < actualBlockSize; ++j) {
    centering.row(j) = arma::mean(eps.rows(j, j + center_end), 0);
  }
  
  // Replicate and trim centering matrix
  arma::mat centering_full(final_size, N, arma::fill::none);
  for (int i = 0; i < nBlock; ++i) {
    int start_row = i * actualBlockSize;
    int end_row = std::min(start_row + actualBlockSize - 1, final_size - 1);
    int rows_to_copy = end_row - start_row + 1;
    centering_full.rows(start_row, end_row) = centering.rows(0, rows_to_copy - 1);
  }
  
  // Pre-compute centering for instruments if provided
  arma::mat Mcentering_full;
  if (has_M) {
    arma::mat Mcentering(actualBlockSize, k, arma::fill::zeros);
    for (int j = 0; j < actualBlockSize; ++j) {
      arma::mat subM = M_mat.rows(j, j + center_end);
      for (int col = 0; col < k; ++col) {
        arma::vec col_vec = subM.col(col);
        arma::uvec nonzero_idx = arma::find(col_vec != 0.0);
        if (nonzero_idx.n_elem > 0) {
          Mcentering(j, col) = arma::mean(col_vec.elem(nonzero_idx));
        }
      }
    }
    
    // Replicate and trim
    Mcentering_full.set_size(final_size, k);
    for (int i = 0; i < nBlock; ++i) {
      int start_row = i * actualBlockSize;
      int end_row = std::min(start_row + actualBlockSize - 1, final_size - 1);
      int rows_to_copy = end_row - start_row + 1;
      Mcentering_full.rows(start_row, end_row) = Mcentering.rows(0, rows_to_copy - 1);
    }
  }
  
  // Draw random block indices once
  arma::ivec index = arma::randi<arma::ivec>(nBlock, arma::distr_param(0, nBlocks_avail - 1));
  
  // Pre-allocate bootstrap matrices
  arma::mat U_boot(final_size, N, arma::fill::none);
  arma::mat M_boot;
  if (has_M) {
    M_boot.set_size(final_size, k);
  }
  
  // Fill bootstrap matrices by concatenating selected blocks
  for (int j = 0; j < nBlock; ++j) {
    int start_row = j * actualBlockSize;
    int end_row = std::min(start_row + actualBlockSize - 1, final_size - 1);
    int rows_to_copy = end_row - start_row + 1;
    
    U_boot.rows(start_row, end_row) = Blocks.slice(index(j)).rows(0, rows_to_copy - 1);
    
    if (has_M) {
      M_boot.rows(start_row, end_row) = MBlocks.slice(index(j)).rows(0, rows_to_copy - 1);
    }
  }
  
  // Center the bootstrapped residuals
  arma::mat eps_boot = U_boot - centering_full;
  
  // Center the bootstrapped instruments (only non-zero values)
  if (has_M) {
    for (int col = 0; col < k; ++col) {
      for (int i = 0; i < final_size; ++i) {
        if (M_boot(i, col) != 0.0) {
          M_boot(i, col) -= Mcentering_full(i, col);
        }
      }
    }
  } else {
    M_boot.set_size(0, 0);  // Empty matrix
  }
  
  // Return results as struct
  MBBVARResult result;
  result.eps_boot = eps_boot;
  result.M_boot = M_boot;
  return result;
}

// Overload: Direct matrix version (thread-safe, no Rcpp::wrap/as overhead)
// This version assumes M is always provided (for use in bootstrap loops)
MBBVARResult fmbb_var_cpp(const arma::mat& eps, 
                          int lags, 
                          int BlockSize, 
                          const arma::mat& M) {
         
  const int Tp = eps.n_rows;  // Implicit conversion from arma::uword to int
  const int N = eps.n_cols;   // Implicit conversion from arma::uword to int
  const int T = Tp + lags;
  
  // Compute block size if not provided or invalid
  int actualBlockSize = BlockSize;
  if (actualBlockSize <= 0) {
    actualBlockSize = static_cast<int>(std::round(5.03 * std::pow(Tp, 0.25)));
  }
  
  // Number of blocks needed
  const int nBlock = static_cast<int>(std::ceil(static_cast<double>(T - lags) / actualBlockSize));
  const int nBlocks_avail = T - lags - actualBlockSize + 1;
  const int final_size = T - lags;
  
  // M is always provided in this overload
  const int k = M.n_cols;  // Implicit conversion from arma::uword to int
  
  // Create blocks for residuals
  arma::cube Blocks(actualBlockSize, N, nBlocks_avail, arma::fill::none);
  for (int j = 0; j < nBlocks_avail; ++j) {
    Blocks.slice(j) = eps.rows(j, j + actualBlockSize - 1);
  }
  
  // Create blocks for instruments
  arma::cube MBlocks(actualBlockSize, k, nBlocks_avail, arma::fill::none);
  for (int j = 0; j < nBlocks_avail; ++j) {
    MBlocks.slice(j) = M.rows(j, j + actualBlockSize - 1);
  }
  
  // Pre-compute centering for residuals
  arma::mat centering(actualBlockSize, N, arma::fill::none);
  const int center_end = T - lags - actualBlockSize;
  for (int j = 0; j < actualBlockSize; ++j) {
    centering.row(j) = arma::mean(eps.rows(j, j + center_end), 0);
  }
  
  // Replicate and trim centering matrix
  arma::mat centering_full(final_size, N, arma::fill::none);
  for (int i = 0; i < nBlock; ++i) {
    int start_row = i * actualBlockSize;
    int end_row = std::min(start_row + actualBlockSize - 1, final_size - 1);
    int rows_to_copy = end_row - start_row + 1;
    centering_full.rows(start_row, end_row) = centering.rows(0, rows_to_copy - 1);
  }
  
  // Pre-compute centering for instruments
  arma::mat Mcentering(actualBlockSize, k, arma::fill::zeros);
  for (int j = 0; j < actualBlockSize; ++j) {
    arma::mat subM = M.rows(j, j + center_end);
    for (int col = 0; col < k; ++col) {
      arma::vec col_vec = subM.col(col);
      arma::uvec nonzero_idx = arma::find(col_vec != 0.0);
      if (nonzero_idx.n_elem > 0) {
        Mcentering(j, col) = arma::mean(col_vec.elem(nonzero_idx));
      }
    }
  }
  
  // Replicate and trim
  arma::mat Mcentering_full(final_size, k, arma::fill::none);
  for (int i = 0; i < nBlock; ++i) {
    int start_row = i * actualBlockSize;
    int end_row = std::min(start_row + actualBlockSize - 1, final_size - 1);
    int rows_to_copy = end_row - start_row + 1;
    Mcentering_full.rows(start_row, end_row) = Mcentering.rows(0, rows_to_copy - 1);
  }
  
  // Draw random block indices once
  arma::ivec index = arma::randi<arma::ivec>(nBlock, arma::distr_param(0, nBlocks_avail - 1));
  
  // Pre-allocate bootstrap matrices
  arma::mat U_boot(final_size, N, arma::fill::none);
  arma::mat M_boot(final_size, k, arma::fill::none);
  
  // Fill bootstrap matrices by concatenating selected blocks
  for (int j = 0; j < nBlock; ++j) {
    int start_row = j * actualBlockSize;
    int end_row = std::min(start_row + actualBlockSize - 1, final_size - 1);
    int rows_to_copy = end_row - start_row + 1;
    
    U_boot.rows(start_row, end_row) = Blocks.slice(index(j)).rows(0, rows_to_copy - 1);
    M_boot.rows(start_row, end_row) = MBlocks.slice(index(j)).rows(0, rows_to_copy - 1);
  }
  
  // Center the bootstrapped residuals
  arma::mat eps_boot = U_boot - centering_full;
  
  // Center the bootstrapped instruments (only non-zero values)
  for (int col = 0; col < k; ++col) {
    for (int i = 0; i < final_size; ++i) {
      if (M_boot(i, col) != 0.0) {
        M_boot(i, col) -= Mcentering_full(i, col);
      }
    }
  }
  
  // Return results as struct
  MBBVARResult result;
  result.eps_boot = eps_boot;
  result.M_boot = M_boot;
  return result;
}

//' Moving Block Bootstrap for VAR Residuals and Instruments
//' 
//' Computes a moving-block bootstrap sample for the reduced form errors
//' of a VAR(p) model. Optionally bootstraps a matrix of instruments M.
//' 
//' @param eps (T-p) x N matrix of VAR residuals
//' @param lags Integer lag order p of the VAR model
//' @param BlockSize Integer block size for moving block bootstrap. 
//'   If 0 or negative, automatically computed as round(5.03 * T^0.25)
//' @param M Optional (T-p) x K matrix of instruments to bootstrap (default = NULL)
//' 
//' @return A list containing:
//'   \itemize{
//'     \item eps_boot: (T-p) x N matrix of bootstrapped residuals (centered)
//'     \item M_boot: (T-p) x K matrix of bootstrapped instruments (centered), 
//'                   or empty matrix if M not provided
//'   }
//' 
//' @details
//' The moving block bootstrap preserves the temporal dependence structure
//' in the residuals. The default block size follows the recommendation from
//' "Proxy SVARs: Asymptotic Theory, Bootstrap Inference, and the Effects 
//' of Income Tax Changes in the United States".
//' 
//' Both residuals and instruments are centered to preserve their sample properties.
//' For instruments, only non-zero values are centered (to handle sparse instruments).
//' 
//' @examples
//' \dontrun{
//' # Bootstrap residuals only
//' eps <- matrix(rnorm(200), 100, 2)
//' result <- fmbb_var(eps, lags = 2, BlockSize = 0)
//' eps_boot <- result$eps_boot
//' 
//' # Bootstrap residuals and instruments
//' M <- matrix(rnorm(100), 100, 1)
//' result <- fmbb_var(eps, lags = 2, BlockSize = 0, M)
//' M_boot <- result$M_boot
//' }
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List fmbb_var(const arma::mat& eps, 
                    int lags, 
                    int BlockSize, 
                    Rcpp::Nullable<arma::mat> M = R_NilValue) {
  // Call the C++ function
  MBBVARResult result = fmbb_var_cpp(eps, lags, BlockSize, M);
  
  // Return results as a list for R
  return Rcpp::List::create(
    Rcpp::Named("eps_boot") = result.eps_boot,
    Rcpp::Named("M_boot") = result.M_boot
  );
}
