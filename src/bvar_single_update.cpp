#include <RcppArmadillo.h>
#include <factorstochvol.h>
#include "sample_coefficients.h"
#include "sample_sigma.h"
#include "sample_shrink.h"

void single_bvar(arma::mat& PHI, arma::mat& PHI_priorvariance, arma::vec& a,  
                 arma::vec& c,arma::mat& psi, arma::mat& lambda, arma::vec& xi, 
                 arma::vec& varpi, arma::mat& U, arma::vec& u, 
                 arma::vec& U_priorvariance, double& a_U, const double b_U, 
                 double& c_U, arma::vec& psi_U, arma::vec& lambda_U, 
                 double& xi_U, double& varpi_U, arma::mat& armafacload,  
                 arma::mat& armafac, arma::mat& armatau2, 
                 arma::vec& armalambda2, arma::mat& logvar, arma::vec& logvar0, 
                 arma::mat& sv_para, arma::umat& mixind,  
                 const int rep, const int burnin, const arma::mat& Y, 
                 const arma::mat& X, const int K, const int T, const int M, 
                 const arma::imat& i_mat, const Rcpp::IntegerVector groups, 
                 const arma::mat& PHI_priormean, const std::string priorPHI,
                 const bool isIntercept, const arma::rowvec& priorIntercept,
                 const bool a_hyper, const bool DL_plus, const arma::vec& b, 
                 const arma::vec& a_vec, const arma::vec& a_weight, 
                 const arma::vec& norm_consts, const arma::vec& c_vec, 
                 const bool c_rel_a, const std::string GT_priorkernel, 
                 const double GT_vs, const arma::mat& tau_0, 
                 const arma::mat& tau_1, const double SSVS_s_a, 
                 const double SSVS_s_b, const arma::mat& HMP_srmat, 
                 const arma::mat& V_i_prep, const bool huge, 
                 const std::string sigma_type, const arma::uvec& U_upper_indices,
                 const arma::uvec& ind_U, const std::string priorU, 
                 const bool a_hyper_U, const bool DL_plus_U,
                 const arma::vec& a_vec_U, const arma::vec& a_weight_U,
                 const arma::vec& norm_consts_U, const arma::vec& c_vec_U, 
                 const bool c_rel_a_U, const std::string GT_priorkernel_U, 
                 const double GT_vs_U,   
                 const arma::vec& tau_0_U, const arma::vec& tau_1_U,
                 const double SSVS_s_a_U, const double SSVS_s_b_U,
                 const arma::vec& s_r_3, const Rcpp::NumericVector aShrink, 
                 const Rcpp::NumericVector cShrink, const Rcpp::NumericVector dShrink,
                 const double factor_facloadtol, 
                 const arma::imat& factor_armarestr, 
                 const arma::uvec& armafacloadtunrestrictedelements,
                 const arma::irowvec& nonzerospercol, 
                 const arma::icolvec& nonzerosperrow,
                 const Rcpp::NumericVector sv_priorh0, const bool ngprior, 
                 const bool columnwise,
                 const Rcpp::NumericMatrix factor_homoskedastic, 
                 const double factor_offset,
                 const Rcpp::NumericVector heteroscedastic, 
                 const int factor_interweaving,
                 const stochvol::ExpertSpec_FastSV expert_sv, 
                 const stochvol::ExpertSpec_FastSV expert_fac,
                 const std::vector<stochvol::PriorSpec> prior_specs, 
                 const double B011inv,
                 const Rcpp::NumericMatrix priorHomoscedastic, 
                 const Rcpp::NumericVector sv_offset,
                 const double abs_floor, const bool stable, 
                 const double stable_lim
){
  
  //----0) Sample hyperparameters of hierarchical priors (prior variances V_i)
  
  arma::mat PHI_diff = PHI.head_rows(K) - PHI_priormean.head_rows(K);  
  if(priorPHI == "DL" || priorPHI== "R2D2" || priorPHI=="SSVS" ||
     priorPHI =="HS" || priorPHI== "NG" || priorPHI == "GT" || priorPHI == "HMP"){
    if(!(priorPHI == "HMP" || (priorPHI == "SSVS" && !a_hyper )) || rep >= 0.1*burnin) {
      sample_V_i_PHI(PHI_priorvariance, a, c, lambda, psi, xi, varpi, 
                     priorPHI, PHI_diff, b, a_vec, a_weight, c_vec, c_rel_a, 
                     a_hyper, norm_consts, GT_vs, GT_priorkernel, tau_0, tau_1, 
                     SSVS_s_a, SSVS_s_b, HMP_srmat, V_i_prep,  DL_plus, groups, i_mat);
    }
  }
  
  //----1) Draw PHI (reduced form VAR coefficients)
  arma::mat coef_priorvariance(K+isIntercept, M);
  coef_priorvariance.head_rows(K) = PHI_priorvariance;
  if(isIntercept) {
    coef_priorvariance.row(K) = priorIntercept; // priorIntercept must be a rowvec!
  }
  if(K>0){
    try{
      sample_PHI(PHI, PHI_priormean, Y, X, coef_priorvariance, abs_floor, 
                 stable, stable_lim, sigma_type, logvar, 
                 armafacload, armafac, huge, U, M);
    } catch (const std::exception &e) {
      Rcpp::stop("sample_PHI() failed in run %i: %s", rep+1, e.what());
    } catch (...) {
      Rprintf("\nsample_PHI() failed in run %i: Rethrowing exception:", rep+1);
      throw;
    }
    if(!PHI.is_finite()){
      Rcpp::stop("non-finite PHI in rep %i.", rep+1);
    }
  }
  
  const arma::mat resid = Y - X*PHI;
  
  //----2) Draw Sigma_t
  if(sigma_type == "cholesky") {
    // Draw hyperparameters of hierarchical priors
    u = U(U_upper_indices);
    if(priorU == "DL" || priorU == "GT" || priorU == "SSVS" || priorU == "HS" ||
       priorU == "HMP") {
      sample_V_i_U(U_priorvariance, a_U, c_U, lambda_U, psi_U, xi_U, varpi_U, 
                   priorU, u, b_U, a_vec_U, a_weight_U, c_vec_U, c_rel_a_U, 
                   a_hyper_U, norm_consts_U, GT_vs_U, GT_priorkernel_U, tau_0_U, 
                   tau_1_U, SSVS_s_a_U, SSVS_s_b_U, s_r_3, DL_plus_U,
                   ind_U); 
    }
  }
  
  try{
    sample_sigma(armafacload, armafac, logvar, logvar0, sv_para, armatau2,
                 armalambda2, mixind, U, sigma_type, resid, factor_facloadtol,
                 factor_armarestr, armafacloadtunrestrictedelements,
                 nonzerospercol, nonzerosperrow, sv_priorh0, ngprior, columnwise,
                 aShrink, cShrink, dShrink, factor_homoskedastic, factor_offset, 
                 heteroscedastic, factor_interweaving, expert_sv, expert_fac,
                 prior_specs, B011inv, sv_offset, U_priorvariance, 
                 priorHomoscedastic, abs_floor, T, M, rep);
    
  } catch(const std::exception &e){
    Rcpp::stop("sample_sigma() failed in run %i: %s", rep+1, e.what());
  } catch (...) {
    Rprintf("\nsample_sigma() failed in run %i. Rethrowing exception:", rep+1);
    throw;
  }
}