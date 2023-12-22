#ifndef SAMPLE_COEFFICIENTS_H
#define SAMPLE_COEFFICIENTS_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

arma::vec mvrnorm1(arma::vec& mu, arma::mat& Sigma, double tol = 1e-06);

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

// void sample_V_i_DL_deprecated(arma::vec& V_i, const arma::vec coefs, double& a ,
//                    const arma::vec a_vec, const arma::vec prep1,
//                    const arma::vec prep2, double& zeta, arma::vec& psi,
//                    arma::vec& theta, arma::uvec ind, const bool hyper,
//                    const int method, const double tol);

// double do_rgig1(double lambda, double chi, double psi);
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

// void sample_V_i_NG(arma::vec& V_i, const arma::vec coefs, arma::vec& theta_tilde,
//                    double& zeta, double& a , const arma::vec a_vec,
//                    const double varrho0, const double varrho1, arma::uvec ind,
//                    const bool hyper, const double tol);

void sample_V_i_SSVS_beta(arma::vec& V_i, arma::vec& gammas, arma::vec& p_i,
                          const arma::vec coeffs, const arma::vec tau_0,
                          const arma::vec tau_1, const double s_a, const double s_b,
                          const bool hyper, arma::uvec ind);

// void sample_V_i_SSVS(arma::vec& V_i, arma::vec& gammas, arma::vec& p_i,
//                      const arma::vec& coeffs, const arma::vec& tau_0,
//                      const arma::vec& tau_1, const double& s_a, const double& s_b,
//                      const std::string type);

void sample_V_i_HMP(double& lambda_1, double& lambda_2, arma::vec& V_i, const double& s1,
                    const double& r1, const double& s2, const double& r2,
                    const arma::vec& PHI_diff, const arma::vec& V_i_prep,
                    const int& n_ol, const int& n_cl, const arma::uvec& i_ol,
                    const arma::uvec& i_cl);

void sample_V_i_U_HMP(double& lambda_3, arma::vec& V_i_U, const double& s1,
                      const double& r1, const arma::vec& u);

// void sample_V_i_R2D2(arma::vec& V_i, const arma::vec coefs, double& api,
//                      const arma::vec api_vec, double& zeta, arma::vec& psi,
//                      arma::vec& theta, double& xi, double& b,
//                      const arma::vec b_vec, arma::uvec ind, const bool hyper,
//                      const int method, const std::string kernel);

// void sample_PHI_SL(arma::mat& PHI, const arma::mat& PHI_prior, const arma::mat& Y,
//                    const arma::mat& X, const arma::mat& U, const arma::mat& d_sqrt,
//                    arma::mat& Gamma, const int& K, const int& M,
//                    const double& nu_a, const double& nu_b);
//
// void sample_L_SL(arma::mat& U, arma::mat& Ytilde, const arma::mat& d_sqrt,
//               vec& omega, const double& nu_a, const double& nu_b);

#endif
