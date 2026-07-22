#ifndef FCHOLESKY_IRF_H
#define FCHOLESKY_IRF_H
 
#include <RcppArmadillo.h>
 
 // Function declaration
 arma::cube fCholeskyIRF(const arma::cube& wold, const arma::mat& S);
 
#endif // FCHOLESKY_IRF_H
