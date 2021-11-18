#ifndef SAMPLE_COEFFICIENTS_H
#define SAMPLE_COEFFICIENTS_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

void sample_PHI(arma::mat& PHI, const arma::mat& PHI_prior, const arma::mat& Y,
                const arma::mat& X, const arma::mat& L, const arma::mat& d,
                const arma::mat& V_prior, const int& K, const int& M);

void sample_L(arma::mat& L, arma::mat& Ytilde, arma::vec& V_i, arma::mat& d);

void sample_V_i_DL(arma::vec& V_i, const arma::vec coefs, const double a ,
                   double& zeta, arma::vec& psi, arma::vec& theta);

arma::colvec ddir_prep(arma::colvec x, arma::mat prep1, arma::rowvec prep2);

void sample_DL_hyper(double& a, const arma::vec& theta, const arma::mat& prep1,
                     const arma::rowvec& prep2, const double& zeta,
                     arma::vec& a_vec);

arma::colvec ddir_prep(arma::colvec x, arma::mat prep1, arma::rowvec prep2);

void sample_DL_hyper(double& a, const arma::vec& theta, const arma::mat& prep1,
                     const arma::rowvec& prep2, const double& zeta,
                     arma::vec& a_vec);

double do_rgig1(double lambda, double chi, double psi);

void sample_V_i_SSVS(arma::vec& V_i, arma::vec& gammas, arma::vec& p_i,
                     const arma::vec coeffs, const arma::vec tau0,
                     const arma::vec tau_1, const double s_a, const double s_b);

void sample_V_i_HMP(double& lambda_1, double& lambda_2, arma::vec& V_i, const double s1,
                    const double r1, const double s2, const double r2,
                    const arma::vec PHI_diff, const arma::vec V_i_prep,
                    const int n_ol, const int n_cl, const arma::uvec i_ol,
                    const arma::uvec i_cl);

#endif
