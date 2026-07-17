#ifndef FGENERATEQ_H
#define FGENERATEQ_H

#include <RcppArmadillo.h>

// Fill Q with an N x N random orthonormal matrix using QR decomposition.
// R and G are caller-owned scratch buffers reused by hot loops.
void fGenerateQ_inplace(arma::mat& Q,
                        arma::mat& R,
                        arma::mat& G,
                        arma::uword N);

// Generate an N x N random orthonormal matrix using QR decomposition.
// Sign normalisation: diagonal of R forced positive for a unique decomposition
// (Rubio-Ramirez, Waggoner & Zha 2010).
arma::mat fGenerateQ_cpp(int N);

// R wrapper
arma::mat fGenerateQ(int N);

#endif // FGENERATEQ_H
