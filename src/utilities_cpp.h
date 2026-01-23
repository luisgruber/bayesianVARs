#pragma once
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

void build_sigma(arma::mat& Sigma, arma::mat& Sigma_chol, const bool& factor,
                 const arma::mat& facload, const arma::rowvec& logvar_t,
                 const int& factors, const int& m, const arma::vec u,
                 const bool& return_chol);
