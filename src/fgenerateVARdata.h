#ifndef FGENERATEVARDATA_H
#define FGENERATEVARDATA_H

#include <RcppArmadillo.h>

// Function declaration
arma::mat fgenerateVARdata(const arma::mat& y, 
                           int p, 
                           int c, 
                           const arma::mat& beta, 
                           const arma::mat& residuals);

#endif // FGENERATEVARDATA_H
