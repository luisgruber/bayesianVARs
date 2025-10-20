#ifndef SAMPLE_COEFFICIENTS_H
#define SAMPLE_COEFFICIENTS_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

void sample_PHI_factor(arma::mat& PHI, const arma::mat& PHI_prior,
                       const arma::mat& Y, const arma::mat& X,
                       const arma::mat& logvaridi, const arma::mat& V_prior,
                       const arma::mat& facload, const arma::mat& fac,
                       const bool& huge);

void sample_PHI(arma::mat& PHI, const arma::mat PHI_prior, const arma::mat Y,
                const arma::mat X, const arma::mat U, const arma::mat d_sqrt,
                const arma::mat V_prior, const int M);

void sample_U(arma::mat& U, const arma::mat& Ytilde, const arma::vec& V_i,
              const arma::mat& d_sqrt);

void sample_V_i_DL(arma::vec& V_i, const arma::vec coefs, double& a ,
                   const double b, double& c,
                   const arma::vec a_vec, const arma::vec a_weight,
                   arma::vec& psi, arma::vec& lambda, double& xi, arma::uvec ind,
                   const bool hyper,const arma::vec norm_consts,
                   const double tol, const bool DL_plus,
                   const arma::vec c_vec, const bool c_rel_a);

double do_rgig(double lambda, double chi, double psi);

void sample_V_i_HS(arma::vec& V_i, const arma::vec coefs, arma::vec& theta,
                   double& zeta, arma::vec& nu, double& varpi ,arma::uvec ind);

void sample_V_i_GT(arma::vec& V_i, const arma::vec coefs, arma::vec& psi,
                   arma::vec& lambda, double& xi, double& a, const double b,
                   double& c, arma::uvec ind, const double tol,
                   const std::string priorkernel,
                   const double vs, const arma::vec norm_consts,
                   const arma::vec a_vec, const arma::vec a_weight,
                   const arma::vec c_vec, const bool hyper, const bool c_rel_a);

void sample_V_i_SSVS_beta(arma::vec& V_i, arma::vec& gammas, arma::vec& p_i,
                          const arma::vec coeffs, const arma::vec tau_0,
                          const arma::vec tau_1, const double s_a, const double s_b,
                          const bool hyper, arma::uvec ind);

void sample_V_i_HMP(double& lambda_1, double& lambda_2, arma::vec& V_i, const double& s1,
                    const double& r1, const double& s2, const double& r2,
                    const arma::vec& PHI_diff, const arma::vec& V_i_prep,
                    const int& n_ol, const int& n_cl, const arma::uvec& i_ol,
                    const arma::uvec& i_cl);

void sample_V_i_U_HMP(double& lambda_3, arma::vec& V_i_U, const double& s1,
                      const double& r1, const arma::vec& u);

#endif
