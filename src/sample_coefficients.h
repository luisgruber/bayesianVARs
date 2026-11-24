#ifndef SAMPLE_COEFFICIENTS_H
#define SAMPLE_COEFFICIENTS_H

#include <RcppArmadillo.h>

void sample_PHI(arma::mat& PHI, const arma::mat& PHI_prior,
                const arma::mat& Y, const arma::mat& X,
                const arma::mat& V_prior, const double floor, const bool stable, 
                const double stable_lim, const std::string sigma_type,
                const arma::mat& logvar, const arma::mat& facload, 
                const arma::mat& fac, const bool& huge, 
                const arma::mat& U, const int M);

void sample_U(arma::mat& U, const arma::mat& Ytilde, const arma::vec& V_i, const arma::mat& logvar,
              const double floor);

#endif
