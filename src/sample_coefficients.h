#ifndef SAMPLE_COEFFICIENTS_H
#define SAMPLE_COEFFICIENTS_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

arma::vec mvrnorm1(arma::vec& mu, arma::mat& Sigma, double tol = 1e-06);

void sample_PHI(arma::mat& PHI, const arma::mat& PHI_prior, const arma::mat& Y,
                const arma::mat& X, const arma::mat& L, const arma::mat& d_sqrt,
                const arma::mat& V_prior, const int& K, const int& M, bool subs);

void sample_L(arma::mat& L, arma::mat& Ytilde, const arma::vec& V_i, const arma::mat& d_sqrt);

void sample_V_i_DL(arma::vec& V_i, const arma::vec& coefs, const double& a ,
                   double& zeta, arma::vec& psi, arma::vec& theta); //, bool hyper

void sample_V_i_DL_new(arma::vec& V_i, const arma::vec& coefs, const arma::vec& a ,
                       arma::vec& zeta, arma::vec& psi, arma::vec& theta,
                       arma::ivec& groups, const arma::ivec& i_vec );

arma::colvec ddir_prep(const arma::colvec& x, const arma::vec& prep1,
                       const arma::rowvec& prep2);


void sample_DL_hyper(double& a, const arma::vec& theta, const arma::mat& prep1,
                     const arma::rowvec& prep2, const double& zeta,
                     const arma::vec& a_vec);

void sample_DL_hyper_new(arma::vec& a, const arma::colvec& theta, const arma::vec& prep1,
                         const arma::mat& prep2, const arma::vec& zeta,
                         const arma::vec& a_vec, arma::ivec& groups,
                         const arma::ivec& i_vec );

double do_rgig1(double lambda, double chi, double psi);

void sample_V_i_SSVS(arma::vec& V_i, arma::vec& gammas, arma::vec& p_i,
                     const arma::vec& coeffs, const arma::vec& tau_0,
                     const arma::vec& tau_1, const double& s_a, const double& s_b);

void sample_V_i_HMP(double& lambda_1, double& lambda_2, arma::vec& V_i, const double& s1,
                    const double& r1, const double& s2, const double& r2,
                    const arma::vec& PHI_diff, const arma::vec& V_i_prep,
                    const int& n_ol, const int& n_cl, const arma::uvec& i_ol,
                    const arma::uvec& i_cl);

void sample_V_i_L_HMP(double& lambda_3, arma::vec& V_i_L, const double& s1,
                      const double& r1, const arma::vec& l);

void sample_V_i_R2D2(arma::vec& V_i, const arma::vec& coefs,  const double& api,
                     double& zeta, arma::vec& psi, arma::vec& theta, double& xi,
                     const double& a , const double& b);

void sample_V_i_R2D2_new(arma::vec& V_i, const arma::vec& coefs,  const arma::vec& api,
                         arma::vec& zeta, arma::vec& psi, arma::vec& theta, arma::vec& xi,
                         const arma::vec& a , arma::vec& b,
                         arma::ivec& groups, const arma::ivec& i_vec, const bool b_hyper );

void sample_PHI_SL(arma::mat& PHI, const arma::mat& PHI_prior, const arma::mat& Y,
                   const arma::mat& X, const arma::mat& L, const arma::mat& d_sqrt,
                   arma::mat& Gamma, const int& K, const int& M,
                   const double& nu_a, const double& nu_b);

void sample_L_SL(arma::mat& L, arma::mat& Ytilde, const arma::mat& d_sqrt,
              vec& omega, const double& nu_a, const double& nu_b);

void sample_prior_mean(vec& m_i ,const vec& coefs, const vec& v_i, const double& mu0,
                       const double& b0);

#endif
