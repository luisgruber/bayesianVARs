#include <RcppArmadillo.h>
#include <factorstochvol.h>
#include "sample_coefficients.h"

void sample_sigma_chol(arma::mat& U, arma::mat& logvar, arma::vec& logvar0,
                       arma::umat& mixind,
                       arma::mat& sv_para,
                       const Rcpp::NumericVector sv_offset,
                       const arma::mat& resid, 
                       const arma::vec& U_priorvariance, 
                       const int T, const int M,
                       const Rcpp::NumericVector heteroscedastic,
                       const Rcpp::NumericMatrix priorHomoscedastic,
                       const std::vector<stochvol::PriorSpec> prior_specs,
                       const stochvol::ExpertSpec_FastSV expert_sv,
                       const double abs_floor) {
  sample_U(U, resid, U_priorvariance, logvar, abs_floor);
  const arma::mat str_resid = resid*U; // structural (orthogonalized) residuals
  
  for(int j =0; j<M; j++){
    if(heteroscedastic(j) == false){
      double s_p = priorHomoscedastic(j,1) + 0.5*arma::accu(arma::square(str_resid.col(j)));
      double d_i = 1. / R::rgamma(priorHomoscedastic(j,0)+.5*T, 1./s_p);
      logvar.col(j).fill(std::log(d_i));
    }else if(heteroscedastic(j) == true){
      const arma::mat resid_norm = arma::log(arma::square(str_resid.col(j)) + sv_offset(j));
      arma::vec h_j  = logvar.unsafe_col(j);  // unsafe_col reuses memory, logvar will automatically be overwritten
      arma::uvec mixind_j = mixind.unsafe_col(j);
      double mu = sv_para(0,j),
        phi = sv_para(1,j),
        sigma = sv_para(2,j),
        h0_j = logvar0(j);
      stochvol::update_fast_sv(resid_norm, mu, phi, sigma, h0_j, h_j, mixind_j, prior_specs[j], expert_sv); 
      sv_para.col(j) = arma::colvec({mu, phi, sigma}); 
      logvar0(j) = h0_j;
    }
  }
}

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
                  const int rep) {
  if(sigma_type == "factor") {
    const arma::mat residt = resid.t();
    try{
      factorstochvol::update_fsv(facload, fac, logvar, logvar0, sv_para,
                                 factor_tau2, factor_lambda2, mixind,
                                 residt, facloadtol,  factor_restr,
                                 factor_unrestrictedelements,
                                 factor_nonzerospercol, factor_nonzerosperrow,
                                 factor_sv_priorh0, factor_ngprior,
                                 factor_columnwise, factor_aShrink, factor_cShrink,
                                 factor_dShrink, factor_homoskedastic,
                                 factor_offset, heteroscedastic,
                                 factor_interweaving, expert_sv,  expert_fac,
                                 prior_specs, factor_B011inv, true,  false, rep); 
    } catch(const std::exception &e){
      Rcpp::stop("factorstochvol::update_fsv() failed: %s", e.what());
    } catch (...) {
      Rprintf("\nfactorstochvol::update_fsv() failed. Rethrowing exception:");
      throw;
    }
  } else if (sigma_type == "cholesky") {
    try{
      sample_sigma_chol(U, logvar, logvar0, mixind, sv_para, chol_sv_offset,
                        resid, chol_U_priorvariance, T, M, heteroscedastic,
                        chol_priorHomoscedastic, prior_specs, expert_sv,
                        chol_abs_floor);
    } catch(const std::exception &e){
      Rcpp::stop("sample_sigma_chol() failed: %s", e.what());
    } catch (...) {
      Rprintf("\nsample_sigma_chol() failed. Rethrowing exception:");
      throw;
    }
  } else {
    Rcpp::stop("sample_sigma() unaware of 'sigma_type'" + sigma_type);
  }
}
