#ifndef UTILITIES_CPP_H
#define UTILITIES_CPP_H

#include <RcppArmadillo.h>

void build_sigma(arma::mat& Sigma, arma::mat& Sigma_chol, const bool& factor,
                 const arma::mat& facload, const arma::rowvec& logvar_t,
                 const int& factors, const int& m, const arma::vec u,
                 const bool& return_chol);

void abs_floor(double& x, const double floor);

void abs_floor2(double& x, const double y, const double floor);

bool is_stable(const arma::mat& PHI, const int M, const double lim);

void stable_PHI(arma::mat& PHI, const arma::mat& PHI_old, const int M, const double lim);

#endif