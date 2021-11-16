#ifndef SAMPLE_COEFFICIENTS_H
#define SAMPLE_COEFFICIENTS_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

void sample_PHI(arma::mat& PHI, const arma::mat& PHI_prior, const arma::mat& Y,
                const arma::mat& X, const arma::mat& L, const arma::mat& d,
                const arma::mat& V_prior, const int& K, const int& M);

void sample_L(arma::mat& L, arma::mat& Ytilde, arma::vec& V_i, arma::mat& d);

void sample_V_i_DL(arma::vec& V_i, const arma::vec& coefs, const double& a ,
                   double& zeta, arma::vec& psi, arma::vec&theta);

double do_rgig1(double lambda, double chi, double psi);

#endif
