#include <RcppArmadillo.h>
#include <stochvol.h>
#include <progress.hpp>
#include <Rcpp/Benchmark/Timer.h>
#include "sample_coefficients.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
List bvar_cpp(const arma::mat Y,
              const arma::mat X,
              const int M,
              const int T,
              const int K,
              const int draws,
              const int burnin,
              const int intercept,
              const arma::vec priorIntercept,
              arma::mat PHI,
              arma::mat PHI0, //??? no const
              const List priorPHI_in,
              const List priorL_in,
              arma::mat L,
              const bool SV,
              const arma::vec priorHomoscedastic,
              const List sv_spec,
              arma::mat h,
              arma::mat sv_para,
              const arma::imat i_mat,
              const arma::ivec i_vec,
              const bool progressbar,
              const double PHI_tol,
              const double L_tol){

//-------------------------Preliminaries--------------------------------------//

  const double n = K*M; // number of VAR coefficients without intercept
  arma::mat PHI_diff(K+intercept,M); // will hold PHI - PHI0
  const arma::uvec i_ocl= arma::find(i_vec != 0.); // indicator for all coefficients except intercept
  const arma::uvec i_i= arma::find(i_vec == 0.); // indicator for intercepts

  //////////////////
  // indicators for coefficients without intercept
  const arma::ivec i_vec_small=i_vec(i_ocl);
  // indicators for cross-/own-lags (relevant for HM prior)
  const arma::uvec i_ol = arma::find(i_vec_small > 0);
  const arma::uvec i_cl = arma::find(i_vec_small < 0);
  const int n_ol = i_ol.size(); // nr. of ownlags
  const int n_cl = i_cl.size(); // nr. of crosslags
  /////////////////

//--------------------Initialization of hyperparameters-----------------------//

//---- PHI
  std::string priorPHI = priorPHI_in["prior"];
  // V_i holds prior variances (without intercepts)
  arma::vec V_i(n);

  arma::vec V_i_long(n+M*intercept); // ??? V_i plus intercept prior variances
  V_i_long(i_i) = priorIntercept;

  if(priorPHI == "normal"){
    arma::vec V_i_in = priorPHI_in["V_i"];
    //in case of 'normal' V_i is fixed at user specified values
    V_i = V_i_in;
  }

  //---- GL prior on PHI

  //  sub-groups for semi-global-local GL priors
  const int n_groups = priorPHI_in["n_groups"];
  arma::ivec groups = priorPHI_in["groups"];

  /////////// GL initilization

  arma::vec lambda(n, fill::value(1/static_cast<double>(n)));
  arma::vec psi(n, fill::ones);
  arma::vec xi(n_groups, fill::value(0.5));
  arma::vec a = priorPHI_in["a"];
  arma::vec b = priorPHI_in["b"];
  arma::vec c = priorPHI_in["c"];
  const double GL_tol = priorPHI_in["GL_tol"];

  ///////////

  //---- DL prior on PHI
//?  const double DL_tol = priorPHI_in["DL_tol"];
  const int DL_method = priorPHI_in["DL_method"];
  const bool DL_hyper = priorPHI_in["DL_hyper"];
//?  arma::vec DL_a = priorPHI_in["DL_a"];
//?  arma::vec DL_b = priorPHI_in["DL_b"];
//?  arma::vec DL_c = priorPHI_in["DL_c"];
  const bool DL_plus = priorPHI_in["DL_plus"];
  const arma::vec a_vec = priorPHI_in["a_vec"];
  const arma::vec a_weight = priorPHI_in["a_weight"];
  const arma::mat prep2 = priorPHI_in["prep2"]; // DL_deprecated
  const arma::vec prep1 = priorPHI_in["prep1"]; // DL_deprecated
  const arma::vec norm_consts = priorPHI_in["norm_consts"];

  // initialization of scaling, global, and local parameters
//?  arma::vec psi(n);psi.fill(1.0);
  arma::vec zeta(n_groups); zeta.fill(10);
//?  arma::vec lambda(n);lambda.fill(1/static_cast<double>(n));
//?  arma::vec DL_xi(n_groups, fill::value(0.5));
  if(priorPHI == "DL" ){
    V_i = psi % lambda % lambda * zeta(0) * zeta(0);
  }

  //----GT on PHI (GT refers to gamma type priors a la normal gamma or r2d2)
  const double GT_vs = priorPHI_in["GT_vs"];
  const std::string GT_priorkernel = priorPHI_in["GT_priorkernel"];
  if(priorPHI == "GT"){
    if(GT_priorkernel == "exponential"){
      V_i = psi%lambda/2;
    }else if(GT_priorkernel == "normal"){
      V_i = lambda;
    }
  }

  //----R2D2 on PHI
  // initialization of scaling, global, and local parameters
//?  vec xi(n_groups, fill::ones);
//?  arma::vec theta_r2d2(n); theta_r2d2.fill(1/static_cast<double>(n));
//?  vec zeta_r2d2(n_groups); zeta_r2d2.fill(10);
//?  arma::vec psi_r2d2(n); psi_r2d2.fill(1/static_cast<double>(n));

  //  bool R2D2_hyper;
//?  bool R2D2_hyper = priorPHI_in["R2D2_hyper"];
//?  arma::vec api = priorPHI_in["R2D2_api"];
//?  arma::vec b_r2d2 = priorPHI_in["R2D2_b"];
//?  const arma::vec api_vec = priorPHI_in["api_vec"];
//?  const arma::vec b_r2d2_vec = priorPHI_in["b_vec"];
//?  const int R2D2_method = priorPHI_in["R2D2_method"];
//?  std::string R2D2_kernel = priorPHI_in["R2D2_kernel"];
//?  if(priorPHI == "R2D2"){
    // to do: c = 1/2*api, vs
//?    if(R2D2_kernel=="laplace"){
//?      V_i = psi_r2d2%theta_r2d2*zeta_r2d2(0)/2;
//?    }else if(R2D2_kernel == "normal"){
//?      psi_r2d2.fill(1.0); // needed for sample_V_i_R2D2
//?      V_i = theta_r2d2*zeta_r2d2(0);
//?    }
//?    V_i = psi % lambda;
//?  }

  //---- Horseshoe on PHI
  // initialization of global, local and auxiliary scaling parameters
  arma::vec theta_hs(n); theta_hs.fill(1/static_cast<double>(n));
  arma::vec zeta_hs(n_groups); zeta_hs.fill(10);
  arma::vec nu(n); nu.fill(1);
  arma::vec varpi(n_groups); varpi.fill(1);
  if(priorPHI == "HS"){
    V_i = theta_hs*zeta_hs;
  }

  //---- NG on PHI
  //theta_ng, zeta_ng(j), NG_a, varrho0, varrho1
//?  arma::vec theta_ng(n); theta_ng.fill(0.1);
//?  arma::vec zeta_ng(n_groups); zeta_ng.fill(10);
//?  arma::vec NG_a = priorPHI_in["NG_a"];
//?  const double varrho0 = priorPHI_in["NG_varrho0"];
//?  const double varrho1 = priorPHI_in["NG_varrho1"];
//?  const arma::vec a_ng_vec = priorPHI_in["NG_a_vec"];
//?  const bool NG_hyper = priorPHI_in["NG_hyper"];
//?  if(priorPHI == "NG"){
//?    //psi =1;
//?    V_i = lambda;
//?  }

  //---- SSVS on PHI
  arma::vec tau_0 = priorPHI_in["SSVS_tau0"];
  arma::vec tau_1 = priorPHI_in["SSVS_tau1"];
  double SSVS_s_a = priorPHI_in["SSVS_s_a"];
  double SSVS_s_b = priorPHI_in["SSVS_s_b"];
  bool SSVS_hyper = priorPHI_in["SSVS_hyper"];
  arma::vec p_i = priorPHI_in["SSVS_p"];

  if(priorPHI == "SSVS"){
    if(SSVS_hyper){
      V_i = tau_1 % tau_1;
    }else{
      V_i = tau_0 % tau_0; // !!! tau_0 % tau_0
    }
  }
  arma::vec gammas(n, fill::zeros);

  //---- HMP on PHI
  double lambda_1=0.04;
  double lambda_2=0.0016;
  arma::vec V_i_prep = priorPHI_in["V_i_prep"];
  arma::vec s_r_1 = priorPHI_in["lambda_1"];
  arma::vec s_r_2 = priorPHI_in["lambda_2"];
  if(priorPHI == "HMP"){
    //V_i_long(i_i) = priorIntercept % V_i_prep(i_i);
    arma::vec V_i_prep_small = V_i_prep(i_ocl);
    V_i(i_ol) = lambda_1*V_i_prep_small(i_ol);
    V_i(i_cl) = lambda_2*V_i_prep_small(i_cl);
  }

  //Fill V_i_long with remaining prior variances
  V_i_long(i_ocl) = V_i; // ???
  arma::mat V_prior = arma::reshape(V_i_long, K + intercept, M);

//-------------- L
  std::string priorL = priorL_in["prior"];

  uvec L_upper_indices = trimatu_ind( size(L),  1);
  arma::vec l = L(L_upper_indices);
  const double n_L = l.size();
  arma::vec V_i_L(n_L);

  if(priorL == "normal"){
    arma::vec V_i_L_in = priorL_in["V_i"];
    V_i_L = V_i_L_in;
  }

  /// Initialize GL parameters
  arma::vec lambda_L(n_L, fill::value(1/static_cast<double>(n_L)));
  arma::vec psi_L(n_L, fill::ones);
  double xi_L = 0.5;
  double a_L = priorL_in["a"];
  double b_L = priorL_in["b"];
  double c_L = priorL_in["c"];
  const double GL_tol_L = priorL_in["GL_tol"];
  ///

  //---- DL prior on L
//?  arma::vec psi_L(n_L); psi_L.fill(1.0);
//?  double zeta_L=10;
//?  arma::vec theta_L(n_L); theta_L.fill(1/static_cast<double>(n_L));

  bool DL_hyper_L = priorL_in["DL_hyper"];
//?  double DL_b_L = priorL_in["DL_b"];

//?  arma::vec prep2_L = priorL_in["prep2"];;
//?  arma::vec prep1_L = priorL_in["prep1"];
  const bool DL_plus_L = priorL_in["DL_plus"];
  const arma::vec a_vec_L = priorL_in["a_vec"];
  const arma::vec a_weight_L = priorL_in["a_weight"];
  const arma::vec norm_consts_L = priorL_in["norm_consts"];
  if(priorL == "DL"){
    V_i_L= psi_L % lambda_L;
  }

  //----GT on L (GT refers to gamma type priors a la normal gamma or r2d2)
  const double GT_vs_L = priorL_in["GT_vs"];
  const std::string GT_priorkernel_L = priorL_in["GT_priorkernel"];
  if(priorL == "GT"){
    if(GT_priorkernel_L == "exponential"){
      V_i_L = psi_L%lambda_L/2;
    }else if(GT_priorkernel_L == "normal"){
      V_i_L = lambda_L;
    }
  }

  //---- Horseshoe on L
  // initialization of global, local and auxiliary scaling parameters
  arma::vec theta_hs_L(n_L); theta_hs_L.fill(1/static_cast<double>(n_L));
  double zeta_hs_L = 10;
  arma::vec nu_L(n_L, fill::ones);
  double varpi_L = 1;
  if(priorL == "HS"){
    V_i_L = theta_hs_L*zeta_hs_L;
  }
  //----R2D2 on L

//?  double xi_L = 1;
//?  arma::vec theta_L_r2d2(n_L); theta_L_r2d2.fill(1/static_cast<double>(n_L));
//?  double zeta_L_r2d2 = 10;
//?  arma::vec psi_L_r2d2(n_L); psi_L_r2d2.fill(1/static_cast<double>(n_L));

//?  bool R2D2_L_hyper = priorL_in["R2D2_hyper"];
//?  arma::vec api_vec_L = priorL_in["api_vec"];
//?  NumericVector b_vec_L_r2d2 = priorL_in["b_vec"];
//?  double b_L_r2d2 = priorL_in["R2D2_b"];
//?  double api_L = priorL_in["R2D2_api"];
//?  if(priorL == "R2D2"){
//?    V_i_L = psi_L_r2d2 % theta_L_r2d2 * zeta_L_r2d2 /2.;
//?  }

  //---- SSVS on L
  bool SSVS_hyper_L = priorL_in["SSVS_hyper"];
  arma::vec p_i_L = priorL_in["SSVS_p"];
  arma::vec tau_0_L = priorL_in["SSVS_tau0"];
  arma::vec tau_1_L = priorL_in["SSVS_tau1"];
  double SSVS_s_a_L = priorL_in["SSVS_s_a"];
  double SSVS_s_b_L = priorL_in["SSVS_s_b"];
  if(priorL == "SSVS"){
    V_i_L = tau_0_L % tau_0_L;
  }
  arma::vec gammas_L(n_L, fill::zeros);

  //---- HMP on L
  double lambda_3 = 0.001;
  NumericVector s_r_3 = priorL_in["lambda_3"];
  if(priorL == "HMP"){
    arma::vec V_i_L_tmp(n_L); V_i_L_tmp.fill(1.0);
    V_i_L= lambda_3*V_i_L_tmp;
  }

  //-----------------------------SV-settings----------------------------------//
  // Import sv_spec
  NumericVector sv_priormu = sv_spec["priormu"] ;
  NumericVector sv_priorphi= sv_spec["priorphi"];
  NumericVector sv_priorsigma2 = sv_spec["priorsigma2"] ;
  NumericVector sv_offset = sv_spec["sv_offset"];

  // create prior specification object for the update_*_sv functions
  using stochvol::PriorSpec;
  const PriorSpec prior_spec = {
    PriorSpec::Latent0(),  // stationary prior distribution on h0
    PriorSpec::Mu(PriorSpec::Normal(sv_priormu[0], sv_priormu[1])),  // N(mu,sd)
    PriorSpec::Phi(PriorSpec::Beta(sv_priorphi[0], sv_priorphi[1])),  //
    PriorSpec::Sigma2(PriorSpec::Gamma(sv_priorsigma2[0], sv_priorsigma2[1]))  //
  };  // stochvol would even offen more :heavy-tailed, leverage, regression

  //expert settings: these are the same settings as the default settings for the
  //idiosyncratic variances from package factorstochvol
  //(https://github.com/gregorkastner/factorstochvol, accessed 2021-11-12)
  using stochvol::ExpertSpec_FastSV;
  const ExpertSpec_FastSV expert_sv {
    true,  // interweave
    stochvol::Parameterization::CENTERED,  // centered baseline always
    1e-8,  // B011inv,
    1e-12,  //B022inv,
    2,  // MHsteps,
    ExpertSpec_FastSV::ProposalSigma2::INDEPENDENCE,  // independece proposal for sigma
    -1,  // MHcontrol unused for independence prior for sigma
    ExpertSpec_FastSV::ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL  // immediately reject (mu,phi,sigma) if proposed phi is outside (-1, 1)
  };

  // initialize mixing indicators
  arma::umat mixind(T,M); mixind.fill(5);
  //initialize d_sqrt
  arma::mat d_sqrt = arma::exp(h/2);


  //------------------------------------STORAGE-------------------------------//

  arma::cube sv_latent_draws(draws, T, M);
  arma::cube sv_para_draws(draws, 4,M);
  arma::cube PHI_draws(draws, (K+intercept), M); // ??? K + intercept
  arma::cube L_draws(draws, M, M);
  arma::mat V_prior_draws(draws, n+M*intercept);

  int phi_hyperparameter_size(0);
  if(priorPHI == "DL" ){
    phi_hyperparameter_size += 2*n_groups + 2*n; // a + xi + n(lambda + psi)
  }else if(priorPHI == "GT"){
    phi_hyperparameter_size += n_groups + 2*n; // xi + n(lambda + psi)
  }else if(priorPHI == "R2D2"){
    phi_hyperparameter_size += 4*n_groups + 2*n; // (b+)xi + zeta + n(theta + psi)
  }else if(priorPHI == "HS"){
    phi_hyperparameter_size += 2*n_groups + 2*n; // zeta + varpi + n(theta + nu)
  }else if (priorPHI == "NG"){
    phi_hyperparameter_size += n + 2*n_groups;
  }else if(priorPHI == "SSVS"){
    phi_hyperparameter_size += 2*n; // n(gammas + p_i)
  }else if(priorPHI == "HMP"){
    phi_hyperparameter_size += 2; // lambda_1 + lambda_2
  }
  arma::mat phi_hyperparameter_draws(draws, phi_hyperparameter_size);

  int l_hyperparameter_size(0);
  if(priorL == "DL" ){
    l_hyperparameter_size += 2. + 2*n_L;
  }else if(priorL == "GT"){
    l_hyperparameter_size += 1. + 2*n_L; // xi + n(lambda + psi)
  }else if(priorL == "R2D2"){
    l_hyperparameter_size += 4. + 2*n_L; // b + api +xi + zeta + n(theta + psi)
  }else if(priorL == "SSVS"){
    l_hyperparameter_size += 2*n_L;
  }else if(priorL == "HMP"){
    l_hyperparameter_size += 1;
  }else if(priorL == "HS"){
    l_hyperparameter_size += 2.+2*n_L;
  }
  arma::mat l_hyperparameter_draws(draws, l_hyperparameter_size);

  //indicator vector needed for DL and R2D2 on L
  arma::uvec ind_L(n_L);
  for(int i=0; i<n_L; ++i){
    ind_L(i)=i;
  }
  //-----------------------------------SAMPLER--------------------------------//

  const int tot = draws + burnin;
  // Initialize progressbar
  //Progress p(tot, progressbar);
  Timer timer;
  timer.step("start");
  for(int rep = 0; rep < tot; rep++){

    //----1) Draw PHI (reduced form VAR coefficients)
    arma::mat PHI_old = PHI;
    try{
      sample_PHI(PHI, PHI0, Y, X, L, d_sqrt, V_prior, M, true);
    } catch(...){
      ::Rf_error("Couldn't sample PHI in rep %i.", rep);
    }

    if(!PHI.is_finite()){
      ::Rf_error("non-finite PHI in rep %i.", rep);
    }

    if(priorPHI == "DL" || priorPHI == "R2D2" || priorPHI == "NG" || priorPHI == "GT"){
      // to do: as function argument
      PHI_diff = PHI - PHI0;
      for (int ii = 0; ii<K; ii++) {
        for (int jj = 0; jj<M; jj++){
          if(PHI_diff(ii,jj) == 0) {

            if(R::rbinom( 1, 0.5 )==0){
              PHI(ii,jj) = PHI0(ii,jj) + PHI_tol;//eps(ii,jj);//PHI_tol;
              PHI_diff(ii,jj) = PHI_tol;
            }else{
              PHI(ii,jj) = PHI0(ii,jj) - PHI_tol;//eps(ii,jj);//PHI_tol;
              PHI_diff(ii,jj) = -PHI_tol;
            }
          }else if(PHI_diff(ii,jj) < PHI_tol && PHI_diff(ii,jj) > 0){ //eps(ii,jj)
            PHI(ii,jj) = PHI0(ii,jj) + PHI_tol;//eps(ii,jj);
            PHI_diff(ii,jj) = PHI_tol;
          }else if (PHI_diff(ii,jj) > -PHI_tol && PHI_diff(ii,jj) < 0){ //-eps(ii,jj)
            PHI(ii,jj) = PHI0(ii,jj) - PHI_tol;//eps(ii,jj);
            PHI_diff(ii,jj) = -PHI_tol;
          }
        }
      }
    }else{
      PHI_diff = PHI - PHI0;
    }

    //----2) Sample hyperparameters of hierarchical priors (prior variances V_i)

    if(priorPHI == "DL" || priorPHI== "R2D2" || priorPHI=="SSVS" ||
       priorPHI =="HS" || priorPHI== "NG" || priorPHI == "GT" ){

      arma::ivec::iterator g;
      int j=0;

      for(g=groups.begin(); g!=groups.end(); ++g){
        arma::uvec ind = arma::find(i_vec_small == *g);
        arma::uvec indplus = arma::find(i_vec == *g);

        if(priorPHI=="DL"){
          if(DL_method==2.){
            sample_V_i_DL(V_i, PHI_diff(i_ocl), a(j),b(j),c(j), a_vec, a_weight, psi,
                          lambda, xi(j), ind, DL_hyper, norm_consts, GL_tol, DL_plus);
          }else if(DL_method==1.){
            sample_V_i_DL_deprecated(V_i, PHI_diff(i_ocl), a(j) , a_vec,
                                     prep1, prep2, zeta(j), psi, lambda, ind,
                                     DL_hyper,1., GL_tol);
          }

        }else if(priorPHI == "GT"){
          sample_V_i_GT(V_i, PHI_diff(i_ocl), psi, lambda, xi(j), a(j), b(j), c(j),
                        ind, GL_tol, GT_priorkernel, GT_vs);
        }else if(priorPHI == "SSVS"){

          if(rep > 0.1*burnin || SSVS_hyper){

            sample_V_i_SSVS_beta(V_i, gammas, p_i, PHI_diff(i_ocl), tau_0,
                                 tau_1, SSVS_s_a, SSVS_s_b, SSVS_hyper, ind);

          }

        }else if(priorPHI == "HS"){

          sample_V_i_HS(V_i, PHI_diff(i_ocl), theta_hs, zeta_hs(j), nu, varpi(j), ind);

        }
//?        else if(priorPHI == "NG"){
//?
//?          sample_V_i_NG(V_i, PHI_diff(i_ocl), theta_ng, zeta_ng(j), NG_a(j),
//?                        a_ng_vec, varrho0, varrho1, ind, NG_hyper,0);
//?          if(!V_i.is_finite()){
//?            ::Rf_error("non-finite V_i in rep %i after group %i.", rep, j+1);
//?          }
//?
//?        }else if(priorPHI=="R2D2"){
//?          sample_V_i_R2D2(V_i, PHI_diff(indplus), api(j), api_vec, zeta_r2d2(j),
//?                        psi_r2d2, theta_r2d2, xi(j), b_r2d2(j), b_r2d2_vec, ind,
//?                        R2D2_hyper, R2D2_method, R2D2_kernel);
//?      }
        j += 1;
      }

    }else if(priorPHI == "HMP"){

      if(rep > 0.1*burnin){
        sample_V_i_HMP(lambda_1, lambda_2, V_i, s_r_1(0), s_r_1(1), s_r_2(0),
                       s_r_2(1), PHI_diff(i_ocl), V_i_prep(i_ocl), n_ol, n_cl,
                       i_ol, i_cl);
      }
    }

    V_i_long(i_ocl) = V_i;

    V_prior = reshape(V_i_long, K+intercept, M);


    //----3) Draw Sigma_t = inv(L)'*D_t*inv(L), where L is upper unitriangular
    //       and D diagonal
    //----3a) Draw free off-diagonal elements in L

    arma::mat resid = Y - X*PHI;

    try{
      sample_L(L, resid, V_i_L, d_sqrt);
    }
    catch(...){
      ::Rf_error("Couldn't sample L in rep %i.", rep);
    }

    if(priorL == "DL" || priorL == "R2D2" || priorL == "GT"){

      for (int jj = 0; jj<M; jj++){
        for (int ii = 0; ii<jj; ii++) {
          if(L(ii,jj) == 0) {
            if(R::rbinom( 1, 0.5 )==0){
              L(ii,jj) = L_tol;
            }else{
              L(ii,jj) = -L_tol;
            }
          }else
            if(L(ii,jj) < L_tol && L(ii,jj) > 0){
              L(ii,jj) = L_tol;
            }else if (L(ii,jj) > -L_tol && L(ii,jj) < 0){
              L(ii,jj) = -L_tol;
            }
        }
      }
    }

    //----3b) Draw hyperparameters of hierarchical priors
    l = L(L_upper_indices);
    if(priorL == "DL" ){


      try{
//?        sample_V_i_DL_deprecated(V_i_L, l, DL_b_L , b_vec, prep1_L, prep2_L, zeta_L, psi_L,
//?                      theta_L, ind_L, DL_hyper_L, 1.,0);
          sample_V_i_DL(V_i_L, l, a_L, b_L, c_L, a_vec_L, a_weight_L,
                        psi_L, lambda_L, xi_L, ind_L, DL_hyper_L, norm_consts_L,
                        GL_tol_L, DL_plus_L);
      } catch (...) {
        ::Rf_error("Couldn't sample V_i_L (DL prior)  in run %i", rep);

      }

    }else if (priorL == "GT"){

      sample_V_i_GT(V_i_L, l, psi_L, lambda_L, xi_L, a_L, b_L, c_L, ind_L,
                    GL_tol_L, GT_priorkernel_L, GT_vs_L);

    }else if(priorL == "SSVS"){

      sample_V_i_SSVS_beta(V_i_L, gammas_L, p_i_L, l, tau_0_L, tau_1_L, SSVS_s_a_L,
                      SSVS_s_b_L, SSVS_hyper_L, ind_L);

    }else if(priorL == "HMP"){

      sample_V_i_L_HMP(lambda_3, V_i_L, s_r_3(0), s_r_3(1), l);
    }else if(priorL == "HS"){

      sample_V_i_HS(V_i_L, l, theta_hs_L, zeta_hs_L, nu_L, varpi_L, ind_L);

    }
//?    else if(priorL == "R2D2"){
//?
//?      sample_V_i_R2D2(V_i_L, l, api_L, api_vec_L, zeta_L_r2d2, psi_L_r2d2,
//?                      theta_L_r2d2, xi_L, b_L_r2d2, b_vec_L_r2d2, ind_L, R2D2_L_hyper, 1., "laplace" );
//?
//?    }

    //----3c) Draw elements of D_t
    //        in case of SV use package stochvol
    arma::mat str_resid = resid*L; // structural (orthogonalized) residuals

    if(SV == false ){
      for(int i =0; i<M; i++){
        double s_p = priorHomoscedastic(i,1) + 0.5*accu(square(str_resid.col(i)));
        double d_i = 1. / R::rgamma(priorHomoscedastic(i,0)+T/2, 1./s_p);
        d_sqrt.col(i).fill(sqrt(d_i));
        h.col(i).fill(log(d_i));
      }
    }else if(SV == true){
      //const arma::mat resid_norm = log(square(str_resid) + sv_offset); // + 1e-40offset??
      for(int j=0; j < M; j++){
        const arma::mat resid_norm = log(square(str_resid.col(j)) + sv_offset(j));
        arma::vec h_j  = h.unsafe_col(j);  // unsafe_col reuses memory, h will automatically be overwritten
        arma::uvec mixind_j = mixind.unsafe_col(j);
        double mu = sv_para(0,j),
          phi = sv_para(1,j),
          sigma = sv_para(2,j),
          h0_j = sv_para(3,j);
        stochvol::update_fast_sv(resid_norm, mu, phi, sigma, h0_j, h_j, mixind_j, prior_spec, expert_sv); //resid_norm.col(j)
        sv_para.col(j) = arma::colvec({mu, phi, sigma, h0_j});
      }
      d_sqrt = exp(h/2);
    }

    //-------Store draws after burnin
    if(rep >= burnin){

      PHI_draws.row(rep-burnin) = PHI;
      L_draws.row(rep-burnin) = L;
      sv_latent_draws.row(rep-burnin) = h;
      sv_para_draws.row((rep-burnin)) = sv_para;
      V_prior_draws.row(rep-burnin) = V_i_long;

      if(priorPHI == "DL" ){

        phi_hyperparameter_draws(rep-burnin, span(0,(n_groups-1))) = a;
        phi_hyperparameter_draws(rep-burnin, span(n_groups,(n_groups+n-1))) = psi.as_row();
        phi_hyperparameter_draws(rep-burnin, span((n_groups+n),(n_groups+2*n-1))) = lambda.as_row();
        if(DL_method==2.){
        phi_hyperparameter_draws(rep-burnin, span((n_groups+2*n),(phi_hyperparameter_size-1.))) = xi;
        }
        if(DL_method==1.){
          phi_hyperparameter_draws(rep-burnin, span((n_groups+2*n),phi_hyperparameter_size-1.)) = zeta;
        }
        //phi_hyperparameter_draws(rep-burnin, span(0,(n_groups-1))) = zeta;
        //phi_hyperparameter_draws(rep-burnin, span(n_groups,(n_groups+n-1))) = psi.as_row();
        //phi_hyperparameter_draws(rep-burnin, span((n_groups+n),(n_groups+2*n-1))) = theta.as_row();
        //phi_hyperparameter_draws(rep-burnin, span((n_groups+2*n),phi_hyperparameter_size-1.)) = DL_a;

      }if(priorPHI == "GT" ){

        phi_hyperparameter_draws(rep-burnin, span(0,(n_groups-1))) = xi;
        phi_hyperparameter_draws(rep-burnin, span(n_groups,(n_groups+n-1))) = psi.as_row();
        phi_hyperparameter_draws(rep-burnin, span((n_groups+n),(phi_hyperparameter_size-1))) = lambda.as_row();

      }else if(priorPHI == "HS"){

        phi_hyperparameter_draws(rep-burnin, span(0, (n_groups -1))) = zeta_hs ;
        phi_hyperparameter_draws(rep-burnin, span(n_groups, (n_groups+n-1))) = theta_hs;
        phi_hyperparameter_draws(rep-burnin, span((n_groups+n), (n_groups+n+n_groups-1))) = varpi;
        phi_hyperparameter_draws(rep-burnin, span((n_groups+n+n_groups),phi_hyperparameter_size-1))= nu;

      }else if(priorPHI == "SSVS"){

        phi_hyperparameter_draws(rep-burnin, span(0, (n-1.))) = gammas;
        phi_hyperparameter_draws(rep-burnin, span(n, (phi_hyperparameter_size-1.))) = p_i;

      }else if(priorPHI == "HMP"){
        phi_hyperparameter_draws(rep-burnin, 0) = lambda_1;
        phi_hyperparameter_draws(rep-burnin, 1) = lambda_2;
      }
//?      else if(priorPHI == "R2D2"){
//?
//?        phi_hyperparameter_draws(rep-burnin, span(0,(n_groups-1))) = zeta_r2d2 ;
//?        phi_hyperparameter_draws(rep-burnin, span(n_groups,(n_groups+n-1))) = trans(psi_r2d2.as_col());
//?        phi_hyperparameter_draws(rep-burnin, span((n_groups+n),(n_groups+2*n-1))) = theta_r2d2.as_row();
//?        phi_hyperparameter_draws(rep-burnin, span((n_groups+2*n),phi_hyperparameter_size-1.-2*n_groups)) = xi ;
//?        phi_hyperparameter_draws(rep-burnin, span((2*n_groups+2*n ),phi_hyperparameter_size -1.-n_groups)) = b_r2d2;
//?        phi_hyperparameter_draws(rep-burnin, span((3*n_groups+2*n ),phi_hyperparameter_size-1.)) = api;
//?
//?      }else if(priorPHI == "NG"){
//?
//?        phi_hyperparameter_draws(rep-burnin, span(0, (n_groups -1))) = zeta_ng ;
//?        phi_hyperparameter_draws(rep-burnin, span(n_groups, (2*n_groups -1))) = NG_a ;
//?        phi_hyperparameter_draws(rep-burnin, span(2*n_groups, (phi_hyperparameter_size-1.))) = theta_ng;
//?
//?      }

      if(priorL == "DL"){

        l_hyperparameter_draws(rep-burnin, 0) = a_L;
        l_hyperparameter_draws(rep-burnin, span(1,n_L)) = psi_L.as_row();
        l_hyperparameter_draws(rep-burnin, span(n_L+1.,2*n_L)) = lambda_L.as_row();
        l_hyperparameter_draws(rep-burnin, l_hyperparameter_size-1.) = xi_L;

      }else if(priorL == "GT"){

        l_hyperparameter_draws(rep-burnin, 0) = xi_L;
        l_hyperparameter_draws(rep-burnin, span(1,n_L)) = psi_L.as_row();
        l_hyperparameter_draws(rep-burnin, span(n_L+1.,2*n_L)) = lambda_L.as_row();

      }else if(priorL == "SSVS"){

        l_hyperparameter_draws(rep-burnin, span(0, (n_L-1.))) = gammas_L;
        l_hyperparameter_draws(rep-burnin, span(n_L, (l_hyperparameter_size-1.))) = p_i_L;

      }else if(priorL == "HMP"){

        l_hyperparameter_draws(rep-burnin, 0) = lambda_3;
      }else if(priorL == "HS"){

        l_hyperparameter_draws(rep-burnin, 0) = zeta_hs_L ;
        l_hyperparameter_draws(rep-burnin, 1) = varpi_L;
        l_hyperparameter_draws(rep-burnin, span(2, n_L+1)) = theta_hs_L;
        l_hyperparameter_draws(rep-burnin, span((n_L+2),(l_hyperparameter_size-1)))= nu_L;

      }

//?      else if(priorL == "R2D2"){
//?
//?        l_hyperparameter_draws(rep-burnin, 0) = zeta_L_r2d2 ;
//?        l_hyperparameter_draws(rep-burnin, span(1,(n_L))) = psi_L_r2d2.as_row();
//?        l_hyperparameter_draws(rep-burnin, span(n_L+1.,(2*n_L))) = theta_L_r2d2.as_row();
//?        l_hyperparameter_draws(rep-burnin, l_hyperparameter_size-3.) = xi_L ;
//?        l_hyperparameter_draws(rep-burnin, l_hyperparameter_size-2.) = b_L_r2d2;
//?        l_hyperparameter_draws(rep-burnin, l_hyperparameter_size-1.) = api_L;
//?
//?      }
    }

    //p.increment();
    if(progressbar){
      //Rprintf("\r %i / %i",
      //        rep+1, tot);
      Rprintf("\r###  %i / %i ### (%3.0f%%) ###",
              rep+1, tot, 100.*(rep+1)/tot);
    }

  }
  timer.step("end");
  NumericVector time(timer);

  List out = List::create(
    Named("PHI") = PHI_draws,
    Named("L") = L_draws,
    Named("sv_latent") = sv_latent_draws,
    Named("sv_para") = sv_para_draws,
    Named("phi_hyperparameter") = phi_hyperparameter_draws,
    Named("l_hyperparameter") = l_hyperparameter_draws,
    Named("bench") = time,
    Named("V_prior") = V_prior_draws,
    Named("V_i") = V_i,
    Named("a") = a,
    Named("b") = b,
    Named("c") = c,
    Named("GT_vs") = GT_vs,
    Named("a_L") = a_L,
    Named("b_L") = b_L,
    Named("c_L") = c_L,
    Named("GT_vs_L") = GT_vs_L
  );

  return out;
}
