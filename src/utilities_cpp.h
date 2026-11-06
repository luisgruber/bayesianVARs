#ifndef UTILITIES_CPP_H
#define UTILITIES_CPP_H

#include <RcppArmadillo.h>

void build_sigma(arma::mat& Sigma, arma::mat& Sigma_chol, const bool& factor,
                 const arma::mat& facload, const arma::rowvec& logvar_t,
                 const int& factors, const int& m, const arma::vec u,
                 const bool& return_chol);

void abs_floor(double& x, const double floor);

void abs_floor2(double& x, const double y, const double floor);

#endif