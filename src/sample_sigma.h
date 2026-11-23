#ifndef SAMPLE_SIGMA_H
#define SAMPLE_SIGMA_H

#include <RcppArmadillo.h>
#include <stochvol.h>

void sample_sigma(arma::mat& facload, arma::mat& fac, arma::mat& logvar,
                  arma::vec& logvar0, arma::mat& sv_para, arma::mat& factor_tau2,
                  arma::vec& factor_lambda2, arma::umat& mixind, arma::mat&U,
                  const std::string sigma_type, const arma::mat& resid, const double facloadtol,
                  const arma::imat& factor_restr, const arma::uvec& factor_unrestrictedelements,
                  const arma::irowvec& factor_nonzerospercol, const arma::icolvec& factor_nonzerosperrow,
                  const Rcpp::NumericVector factor_sv_priorh0, const bool factor_ngprior, const bool factor_columnwise,
                  const Rcpp::NumericVector factor_aShrink, const Rcpp::NumericVector factor_cShrink, 
                  const Rcpp::NumericVector factor_dShrink, const Rcpp::NumericMatrix factor_homoskedastic, 
                  const double factor_offset, const Rcpp::NumericVector heteroscedastic, const int factor_interweaving,
                  const stochvol::ExpertSpec_FastSV expert_sv, const stochvol::ExpertSpec_FastSV expert_fac,
                  const std::vector<stochvol::PriorSpec> prior_specs, const double factor_B011inv,
                  const Rcpp::NumericVector chol_sv_offset, const arma::vec& chol_U_priorvariance,
                  const Rcpp::NumericMatrix chol_priorHomoscedastic, const double chol_abs_floor, 
                  const int T, const int M,
                  const int rep);

#endif