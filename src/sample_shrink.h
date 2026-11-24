#ifndef SAMPLE_SHRINK_H
#define SAMPLE_SHRINK_H

#include <RcppArmadillo.h>

void sample_V_i_U(arma::mat& V_i, double&a, double& c, arma::mat& lambda,
                  arma::mat& psi, double& xi, double& varpi, const std::string prior,
                  const arma::mat& coefs, const double b, const arma::vec& a_vec,
                  const arma::vec& a_weight, const arma::vec& c_vec,
                  const bool c_rel_a, const bool hyper, const arma::vec& norm_consts,
                  const double vs, const std::string priorkernel,
                  const arma::mat& tau_0, const arma::mat& tau_1, const double s_a,
                  const double s_b, const arma::vec& s_r, const bool DL_plus,
                  const arma::uvec& ind);

void sample_V_i_PHI(arma::mat& V_i, arma::vec&a, arma::vec& c, arma::mat& lambda,
                    arma::mat& psi, arma::vec& xi, arma::vec& varpi, const std::string prior,
                    const arma::mat& coefs, const arma::vec& b, const arma::vec& a_vec,
                    const arma::vec& a_weight, const arma::vec& c_vec,
                    const bool c_rel_a, const bool a_hyper, const arma::vec& norm_consts,
                    const double GT_vs, const std::string GT_priorkernel,
                    const arma::mat& tau_0, const arma::mat& tau_1, const double SSVS_s_a,
                    const double SSVS_s_b, const arma::mat& HMP_srmat,
                    const arma::mat& V_i_prep, const bool DL_plus,
                    const Rcpp::IntegerVector groups, const arma::imat& i_mat);

#endif