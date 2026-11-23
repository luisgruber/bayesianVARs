#include <RcppArmadillo.h>
#include <factorstochvol.h>
#include <stochvol.h>
#include "sample_coefficients.h"

void single_bvar(arma::mat& PHI,
                 arma::mat& PHI_priorvariance, arma::vec& a, arma::vec& c, 
                 arma::mat& psi, arma::mat& lambda, arma::vec& xi, 
                 arma::vec& varpi, 
                 double& lambda_1, double& lambda_2,
                 arma::mat& U, arma::vec& u, arma::vec& U_priorvariance,
                 double& a_U, const double b_U, double& c_U, arma::vec& psi_U, 
                 arma::vec& lambda_U, double& xi_U, 
                 double& lambda_3,
                 double& varpi_U, 
                 arma::mat& armafacload, arma::mat& armafac, arma::mat& armatau2, 
                 arma::vec& armalambda2, arma::mat& logvar, arma::vec& logvar0, 
                 arma::mat& sv_para, arma::umat& mixind,
                 const int rep, const int burnin, const arma::mat& Y, 
                 const arma::mat& X, const int K, const int T, const int M, 
                 const int n_ol, const int n_cl,
                 const arma::imat& i_mat,
                 const arma::uvec& i_ol, const arma::uvec& i_cl, const arma::ivec& groups, 
                 const arma::mat& PHI_priormean, const std::string priorPHI,
                 const bool isIntercept, const arma::rowvec& priorIntercept,
                 const bool a_hyper, 
                 const bool DL_plus, const arma::vec& b, 
                 const arma::vec& a_vec, const arma::vec& a_weight, 
                 const arma::vec& norm_consts, const arma::vec& c_vec, 
                 const bool c_rel_a, const std::string GT_priorkernel, 
                 const double GT_vs, const arma::mat& tau_0, const arma::mat& tau_1,
                 const double SSVS_s_a, const double SSVS_s_b, 
                 const arma::vec& s_r_1,
                 const arma::vec& s_r_2, const arma::mat& V_i_prep,
                 const bool huge,   
                 const std::string sigma_type, const arma::uvec& U_upper_indices,
                 const arma::uvec& ind_U,
                 const std::string priorU, const bool a_hyper_U,
                 const bool DL_plus_U,
                 const arma::vec& a_vec_U, const arma::vec& a_weight_U,
                 const arma::vec& norm_consts_U, const arma::vec& c_vec_U, 
                 const bool c_rel_a_U, const std::string GT_priorkernel_U, 
                 const double GT_vs_U,   
                 const arma::vec& tau_0_U, const arma::vec& tau_1_U,
                 const double SSVS_s_a_U, const double SSVS_s_b_U,
                 const arma::vec& s_r_3, const Rcpp::NumericVector aShrink, 
                 const Rcpp::NumericVector cShrink, const Rcpp::NumericVector dShrink,
                 const double factor_facloadtol, 
                 const arma::imat& factor_armarestr, const arma::uvec& armafacloadtunrestrictedelements,
                 const arma::irowvec& nonzerospercol, const arma::icolvec& nonzerosperrow,
                 const Rcpp::NumericVector sv_priorh0, const bool ngprior, const bool columnwise,
                 const Rcpp::NumericMatrix factor_homoskedastic, const double factor_offset,
                 const Rcpp::NumericVector heteroscedastic, const int factor_interweaving,
                 const stochvol::ExpertSpec_FastSV expert_sv, const stochvol::ExpertSpec_FastSV expert_fac,
                 const std::vector<stochvol::PriorSpec> prior_specs, const double B011inv,
                 const Rcpp::NumericMatrix priorHomoscedastic, const Rcpp::NumericVector sv_offset,
                 const double abs_floor, const bool stable, const double stable_lim
                 ){
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

  //----2) Sample hyperparameters of hierarchical priors (prior variances V_i)

  arma::mat PHI_diff = PHI.head_rows(K) - PHI_priormean.head_rows(K);  
  if(priorPHI == "DL" || priorPHI== "R2D2" || priorPHI=="SSVS" ||
     priorPHI =="HS" || priorPHI== "NG" || priorPHI == "GT" ){
    
    int j=0;
    for(arma::ivec::const_iterator g=groups.begin(); g!=groups.end(); ++g){
      arma::uvec ind = arma::find(i_mat == *g);
      
      if(priorPHI=="DL"){
        
        sample_V_i_DL(PHI_priorvariance, PHI_diff, a(j),b(j),c(j), a_vec, a_weight, psi,
                      lambda, xi(j), ind, a_hyper, norm_consts, DL_plus,
                      c_vec, c_rel_a);
        
      }else if(priorPHI == "GT"){ // R2D2 is GT with exponential kernel, NG is GT with normal kernel
        sample_V_i_GT(PHI_priorvariance, PHI_diff, psi, lambda, xi(j), a(j), b(j), c(j),
                      ind, GT_priorkernel, GT_vs, norm_consts,
                      a_vec, a_weight, c_vec, a_hyper, c_rel_a);
      }else if(priorPHI == "SSVS"){
        
        if(rep >= 0.1*burnin || a_hyper){
          
          sample_V_i_SSVS_beta(PHI_priorvariance, lambda, a(j), PHI_diff, tau_0,
                               tau_1, SSVS_s_a, SSVS_s_b, a_hyper, ind);
          
        }
        
      }else if(priorPHI == "HS"){
        
        sample_V_i_HS(PHI_priorvariance, PHI_diff, lambda, xi(j), psi, varpi(j), ind);
        
      }
      j++;
    }
    
  }else if(priorPHI == "HMP"){
    
    if(rep >= 0.1*burnin){
      sample_V_i_HMP(lambda_1, lambda_2, PHI_priorvariance, s_r_1(0), s_r_1(1), s_r_2(0),
                     s_r_2(1), PHI_diff, V_i_prep, n_ol, n_cl,
                     i_ol, i_cl);
    }
  }
  
  arma::mat resid = Y - X*PHI;
  
  //----3) Draw Sigma_t
  
  if(sigma_type == "factor"){
    try{
      factorstochvol::update_fsv(armafacload,
                                 armafac,
                                 logvar,
                                 logvar0,
                                 sv_para,
                                 armatau2,
                                 armalambda2,
                                 mixind,
                                 resid.t(),
                                 factor_facloadtol,
                                 factor_armarestr,
                                 armafacloadtunrestrictedelements,
                                 nonzerospercol,
                                 nonzerosperrow,
                                 sv_priorh0,
                                 ngprior,
                                 columnwise,
                                 aShrink,
                                 cShrink,
                                 dShrink,
                                 factor_homoskedastic,
                                 factor_offset,
                                 heteroscedastic,
                                 factor_interweaving,
                                 expert_sv, // aka expert_idi
                                 expert_fac,
                                 prior_specs,
                                 B011inv,
                                 true, //const bool& samplefac,
                                 false,//const bool& signswitch,
                                 rep);
    } catch(const std::exception &e){
      Rcpp::stop("factorstochvol::update_fsv() failed in run %i: %s", rep+1, e.what());
    } catch (...) {
      Rprintf("\nfactorstochvol::update_fsv() failed in run %i. Rethrowing exception:", rep+1);
      throw;
    }


  }else if(sigma_type == "cholesky"){

    try{
      sample_U(U, resid, U_priorvariance, logvar, abs_floor);
    } catch (const std::exception& e) {
      Rcpp::stop("sample_U() failed in run %i: %s", rep+1, e.what());
    } catch(...){
      Rprintf("\nsample_U() failed in run %i. Rethrowing exception:", rep+1);
      throw;
    }

    //----3b) Draw hyperparameters of hierarchical priors
    u = U(U_upper_indices);
    if(priorU == "DL" ){


      try{
        sample_V_i_DL(U_priorvariance, u, a_U, b_U, c_U, a_vec_U, a_weight_U,
                      psi_U, lambda_U, xi_U, ind_U, a_hyper_U, norm_consts_U, 
                      DL_plus_U, c_vec_U, c_rel_a_U);
      } catch(const std::exception& e) {
        Rcpp::stop("sample_V_i_DL() failed in run %i: %s", rep+1, e.what());
      } catch (...) {
        Rprintf("\nsample_V_i_DL() failed in run %i. Rethrowing exception:", rep+1);
        throw;
      }

    }else if (priorU == "GT"){

      sample_V_i_GT(U_priorvariance, u, psi_U, lambda_U, xi_U, a_U, b_U, c_U, ind_U,
                    GT_priorkernel_U, GT_vs_U, norm_consts_U,
                    a_vec_U, a_weight_U, c_vec_U, a_hyper_U, c_rel_a_U);

    }else if(priorU == "SSVS"){

      sample_V_i_SSVS_beta(U_priorvariance, lambda_U, a_U, u, tau_0_U, tau_1_U, SSVS_s_a_U,
                           SSVS_s_b_U, a_hyper_U, ind_U);

    }else if(priorU == "HMP"){

      sample_V_i_U_HMP(lambda_3, U_priorvariance, s_r_3(0), s_r_3(1), u);
    }else if(priorU == "HS"){

      sample_V_i_HS(U_priorvariance, u, lambda_U, xi_U, psi_U, varpi_U, ind_U);

    }

    //----3c) Draw elements of D_t
    //        in case of SV use package stochvol
    arma::mat str_resid = resid*U; // structural (orthogonalized) residuals

    for(int j =0; j<M; j++){
      if(heteroscedastic(j) == false){
        double s_p = priorHomoscedastic(j,1) + 0.5*accu(square(str_resid.col(j)));
        double d_i = 1. / R::rgamma(priorHomoscedastic(j,0)+.5*T, 1./s_p);
        logvar.col(j).fill(std::log(d_i));
      }else if(heteroscedastic(j) == true){
        //const arma::mat resid_norm = log(square(str_resid) + sv_offset); // + 1e-40offset??
        const arma::mat resid_norm = log(square(str_resid.col(j)) + sv_offset(j));
        arma::vec h_j  = logvar.unsafe_col(j);  // unsafe_col reuses memory, logvar will automatically be overwritten
        arma::uvec mixind_j = mixind.unsafe_col(j);
        double mu = sv_para(0,j),
          phi = sv_para(1,j),
          sigma = sv_para(2,j),
          h0_j = logvar0(j);
        stochvol::update_fast_sv(resid_norm, mu, phi, sigma, h0_j, h_j, mixind_j, prior_specs[j], expert_sv); //resid_norm.col(j)
        sv_para.col(j) = arma::colvec({mu, phi, sigma}); //, h0_j
        logvar0(j) = h0_j;
      }
    }
  }
}