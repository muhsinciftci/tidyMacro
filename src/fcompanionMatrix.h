#ifndef FCOMPANIONMATRIX_H
#define FCOMPANIONMATRIX_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Struct to hold companion matrix results
struct CompanionMatrixResult {
  arma::mat comp;  // Companion matrix (Np x Np)
  int N;           // Number of variables in the VAR system
};

// Internal C++ function (for use in other C++ code)
CompanionMatrixResult fcompanionMatrix_cpp(arma::mat beta, int c, int p);

// R wrapper function (for calling from R)
Rcpp::List fcompanionMatrix(arma::mat beta, int c, int p);

#endif
