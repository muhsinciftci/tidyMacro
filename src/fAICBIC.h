#ifndef FAICBIC_H
#define FAICBIC_H

#include <RcppArmadillo.h>

// Struct to hold AIC/BIC results
struct AICBICResult {
  int aic; // Optimal lag length by Akaike IC
  int bic; // Optimal lag length by Bayesian IC
  int hq;  // Optimal lag length by Hannan-Quinn IC
};

// Internal C++ function (for use in other C++ code)
AICBICResult fAICBIC_cpp(const arma::mat& y, int pmax, int c, 
                         Rcpp::Nullable<arma::mat> exog);

// R wrapper function (for calling from R)
Rcpp::List fAICBIC(const arma::mat& y, int pmax, int c, 
                   Rcpp::Nullable<arma::mat> exog);

#endif // FAICBIC_H
