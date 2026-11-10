#ifndef SAMPLE_COEFFICIENTS_H
#define SAMPLE_COEFFICIENTS_H

#include <RcppArmadillo.h>

void sample_PHI_factor(arma::mat& PHI, const arma::mat& PHI_prior,
                       const arma::mat& Y, const arma::mat& X,
                       const arma::mat& logvaridi, const arma::mat& V_prior,
                       const arma::mat& facload, const arma::mat& fac,
                       const bool& huge, const double floor);

void sample_PHI(arma::mat& PHI, const arma::mat& PHI_prior, const arma::mat& Y,
                const arma::mat& X, const arma::mat& U, const arma::mat& d_sqrt,
                const arma::mat& V_prior, const int M, const double floor);

void sample_U(arma::mat& U, const arma::mat& Ytilde, const arma::vec& V_i,
              const arma::mat& d_sqrt, const double floor);

void sample_V_i_DL(arma::mat& V_i, const arma::mat& coefs, double& a ,
                   const double b, double& c,
                   const arma::vec& a_vec, const arma::vec& a_weight,
                   arma::mat& psi, arma::mat& lambda, double& xi, const arma::uvec& ind,
                   const bool hyper,const arma::vec& norm_consts,
                   const bool DL_plus,
                   const arma::vec& c_vec, const bool c_rel_a);

double do_rgig(double lambda, double chi, double psi);

void sample_V_i_HS(arma::mat& V_i, const arma::mat& coefs, arma::mat& lambda,
                   double& xi, arma::mat& psi, double& varpi, const arma::uvec& ind);

void sample_V_i_GT(arma::mat& V_i, const arma::mat& coefs, arma::mat& psi,
                   arma::mat& lambda, double& xi, double& a, const double b,
                   double& c, const arma::uvec& ind,
                   const std::string priorkernel,
                   const double vs, const arma::vec& norm_consts,
                   const arma::vec& a_vec, const arma::vec& a_weight,
                   const arma::vec& c_vec, const bool hyper, const bool c_rel_a);

void sample_V_i_SSVS_beta(arma::mat& V_i, arma::mat& lambda, double& a,
                          const arma::mat& coefs, const arma::mat& tau_0,
                          const arma::mat& tau_1, const double s_a, const double s_b,
                          const bool hyper, const arma::uvec& ind);

void sample_V_i_HMP(double& lambda_1, double& lambda_2, arma::mat& V_i, const double& s1,
                    const double& r1, const double& s2, const double& r2,
                    const arma::mat& PHI_diff, const arma::mat& V_i_prep,
                    const int& n_ol, const int& n_cl, const arma::uvec& i_ol,
                    const arma::uvec& i_cl);

void sample_V_i_U_HMP(double& lambda_3, arma::vec& V_i_U, const double& s1,
                      const double& r1, const arma::vec& u);

#endif
