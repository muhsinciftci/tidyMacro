#ifndef FGETSHOCK_H
#define FGETSHOCK_H

#include <RcppArmadillo.h>
#include <string>

// Internal C++ function (for use in other C++ code)
// normalize: "unit" = shockSize * (U*Sig^{-1}*s) / (s'*Sig^{-1}*s)
//            "sd"   = (U*Sig^{-1}*s) / sqrt(s'*Sig^{-1}*s)  → unit-variance shock
arma::vec fGetShock_cpp(const arma::mat& residuals, const arma::mat& sigma_full,
                        const arma::vec& s, double shockSize,
                        std::string normalize);

// R wrapper function (for calling from R)
arma::vec fGetShock(const arma::mat& residuals, const arma::mat& sigma_full,
                    const arma::vec& s, double shockSize,
                    std::string normalize);

#endif // FGETSHOCK_H
