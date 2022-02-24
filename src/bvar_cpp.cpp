#include <RcppArmadillo.h>
#include <stochvol.h>
#include <progress.hpp>
#include <Rcpp/Benchmark/Timer.h>
#include "sample_coefficients.h"
//#include "SL.h"

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
              const List persistence_in){

//-------------------------Preliminaries--------------------------------------//

  const double n = K*M; // number of VAR coefficients // ??? without intercept
  arma::mat PHI_diff; // will hold PHI - PHI0
  const arma::uvec i_ol = arma::find(i_vec > 0); // indicator for ownlags
  const arma::uvec i_cl = arma::find(i_vec < 0); // indicator for crosslags
  const arma::uvec i_fol = arma::find(i_vec == 1); // indicator for first ownlags
  const int n_ol = i_ol.size(); // nr. of ownlags
  const int n_cl = i_cl.size(); // nr. of crosslags
  const arma::uvec i_ocl= arma::find(i_vec != 0.); // indicator for all coefficients except intercept
  const arma::uvec i_i= arma::find(i_vec == 0.); // indicator for intercepts

//--------------------Initialization of hyperparameters-----------------------//

//---- Hyperprior on mean of first ownlag coefficients //???
  bool priorPersistence = persistence_in["hyper"];
  double mu0;
  double b0;
  vec m_i;
  if(priorPersistence==true){
    mu0 =  persistence_in["persistence"];
    b0 =  persistence_in["priorPersistence"];
    m_i = PHI0(i_fol);
  }

//---- PHI
  std::string priorPHI = priorPHI_in["prior"];
  // V_i holds prior variances (without intercepts)
  arma::vec V_i(n);

  arma::vec V_i_long(n+M*intercept); // ??? V_i plus intercept prior variances
  V_i_long(i_i) = priorIntercept; // in case of HM prior, these will be scaled later (probably not???)

  if(priorPHI == "normal"){
    arma::vec V_i_in = priorPHI_in["V_i"];
    //in case of 'normal' V_i is fixed at user specified values
    V_i = V_i_in;
  }

  //---- DL prior on PHI
  arma::vec psi(n);psi.fill(1.0);
  double zeta=10;
  arma::vec theta(n);theta.fill(1/static_cast<double>(n));

   double DL_a;
   if(priorPHI == "DL" || priorPHI == "DL_h"){
    double DL_a_in = priorPHI_in["DL_a"];
    DL_a = DL_a_in;

    V_i = psi % theta % theta * zeta * zeta;
   }

  NumericVector a_vec_in;
  NumericVector prep2_in;
  NumericMatrix a_mat_in(1000,n);
  NumericMatrix prep1_in(n,1000);
  if(priorPHI == "DL_h"){
    a_vec_in = priorPHI_in["a_vec"];
    prep2_in = priorPHI_in["prep2"];
    a_mat_in = wrap(priorPHI_in["a_mat"]);
    prep1_in = wrap(priorPHI_in["prep1"]);

  }
  arma::vec a_vec(a_vec_in.begin(), a_vec_in.length(), false);
  arma::rowvec prep2(prep2_in.begin(), prep2_in.length(), false);
  arma::mat a_mat(a_mat_in.begin(), a_mat_in.nrow(), a_mat_in.ncol(), false);
  arma::mat prep1(prep1_in.begin(), prep1_in.nrow(), prep1_in.ncol(), false);

  //----R2D2 on PHI
  double api;
  double a_r2d2;
  double xi = 1;
  arma::vec theta_r2d2(n); theta_r2d2.fill(1/static_cast<double>(n));
  double zeta_r2d2 = 10;
  arma::vec psi_r2d2(n); psi_r2d2.fill(1/static_cast<double>(n));

  double b_r2d2;
  if(priorPHI== "R2D2"){
    double b_r2d2_in = priorPHI_in["R2D2_b"];
    b_r2d2= b_r2d2_in;

    api = 1/(pow(n,(b_r2d2/2)) * pow(T,(b_r2d2/2)) *log(T));
    a_r2d2 = n*api;
    V_i = psi_r2d2%theta_r2d2*zeta_r2d2/2;
  }

  //---- SSVS on PHI
  arma::vec tau_0;
  arma::vec tau_1;
  double SSVS_s_a;
  double SSVS_s_b;
  if(priorPHI == "SSVS"){
    arma::vec tau_0_in = priorPHI_in["SSVS_tau0"];
    tau_0 = tau_0_in;
    arma::vec tau_1_in = priorPHI_in["SSVS_tau1"];
    tau_1 = tau_1_in;
    double SSVS_s_a_in = priorPHI_in["SSVS_s_a"];
    SSVS_s_a = SSVS_s_a_in;
    double SSVS_s_b_in = priorPHI_in["SSVS_s_b"];
    SSVS_s_b = SSVS_s_b_in;

    V_i = tau_0 % tau_0;
  }
  arma::vec gammas(n, fill::zeros);
  arma::vec p_i(n); p_i.fill(0.5);

  //----------------------------------------------SL
  mat Gamma(K,M, fill::ones);
  double nu_a;
  double nu_b;
  if(priorPHI == "SL"){
    double nu_a_in = priorPHI_in["nu_a"];
    double nu_b_in = priorPHI_in["nu_b"];
    nu_a = nu_a_in;
    nu_b = nu_b_in;
  }

  //---- HMP on PHI
  double lambda_1=0.04;
  double lambda_2=0.0016;
  NumericVector V_i_prep_in;
  NumericVector s_r_1_in;
  NumericVector s_r_2_in;
  if(priorPHI == "HMP"){
     V_i_prep_in = priorPHI_in["V_i_prep"];
     s_r_1_in = priorPHI_in["lambda_1"];
     s_r_2_in = priorPHI_in["lambda_2"];
     }
  arma::vec V_i_prep(V_i_prep_in.begin(), V_i_prep_in.length(), false);
  if(priorPHI == "HMP"){
    //V_i_long(i_i) = priorIntercept % V_i_prep(i_i);
    V_i_long(i_ol) = lambda_1*V_i_prep(i_ol);
    V_i_long(i_cl) = lambda_2*V_i_prep(i_cl);
    V_i = V_i_long(i_ocl);
  }
  arma::vec s_r_1(s_r_1_in.begin(), s_r_1_in.length(), false);
  arma::vec s_r_2(s_r_2_in.begin(), s_r_2_in.length(), false);

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
  //---- DL prior on L
  arma::vec psi_L(n_L); psi_L.fill(1.0);
  double zeta_L=10;
  arma::vec theta_L(n_L); theta_L.fill(1/static_cast<double>(n_L));

  double DL_b;
  if(priorL == "DL" || priorL == "DL_h"){
    double DL_b_in = priorL_in["DL_b"];
    DL_b = DL_b_in;

    V_i_L= psi_L% theta_L % theta_L * zeta_L * zeta_L;
  }

  NumericVector b_vec_in;
  NumericVector prep2_L_in;
  NumericMatrix b_mat_in(1000,n_L);
  NumericMatrix prep1_L_in(n_L,1000);
  if(priorL == "DL_h"){
    b_vec_in = priorL_in["b_vec"];
    prep2_L_in = priorL_in["prep2"];
    b_mat_in = wrap(priorL_in["b_mat"]);
    prep1_L_in = wrap(priorL_in["prep1"]);
  }
  arma::vec b_vec(b_vec_in.begin(), b_vec_in.length(), false);
  arma::rowvec prep2_L(prep2_L_in.begin(), prep2_L_in.length(), false);
  arma::mat b_mat(b_mat_in.begin(), b_mat_in.nrow(), b_mat_in.ncol(), false);
  arma::mat prep1_L(prep1_L_in.begin(), prep1_L_in.nrow(), prep1_L_in.ncol(), false);

  //----R2D2 on L
  double api_L;
  double a_L_r2d2;
  double xi_L = 1;
  arma::vec theta_L_r2d2(n_L); theta_L_r2d2.fill(1/static_cast<double>(n_L));
  double zeta_L_r2d2 = 10;
  arma::vec psi_L_r2d2(n_L); psi_L_r2d2.fill(1/static_cast<double>(n_L));

  double b_L_r2d2;
  if(priorL == "R2D2"){
    double b_L_r2d2_in = priorL_in["R2D2_b"];
    b_L_r2d2 = b_L_r2d2_in;
    api_L = 1/(pow(n_L,(b_L_r2d2/2)) * pow(T,(b_L_r2d2/2)) *log(T));
    a_L_r2d2 = n_L*api_L;
    V_i_L = psi_L_r2d2 % theta_L_r2d2 * zeta_L_r2d2 /2.;
  }

  //---- SSVS on L
  arma::vec tau_0_L;
  arma::vec tau_1_L;
  double SSVS_s_a_L;
  double SSVS_s_b_L;
  if(priorL == "SSVS"){
    arma::vec tau_0_L_in = priorL_in["SSVS_tau0"];
    tau_0_L = tau_0_L_in;
    arma::vec tau_1_L_in = priorL_in["SSVS_tau1"];
    tau_1_L = tau_1_L_in;
    double SSVS_s_a_L_in = priorL_in["SSVS_s_a"];
    SSVS_s_a_L = SSVS_s_a_L_in;
    double SSVS_s_b_L_in = priorL_in["SSVS_s_b"];
    SSVS_s_b_L = SSVS_s_b_L_in;

    V_i_L = tau_0_L % tau_0_L;
  }
  arma::vec gammas_L(n_L, fill::zeros);
  arma::vec p_i_L(n_L); p_i_L.fill(0.5);

  //---- HMP on L
  double lambda_3 = 0.001;
  NumericVector s_r_3_in;
  if(priorL == "HMP"){
    s_r_3_in = priorL_in["lambda_3"];
    arma::vec V_i_L_tmp(n_L); V_i_L_tmp.fill(1.0);
    V_i_L= lambda_3*V_i_L_tmp;
  }
  arma::vec s_r_3(s_r_3_in.begin(), s_r_3_in.length(), false);


  //----------------------------------------------SL
  vec omega(n_L, fill::ones);
  double nu_a_L;
  double nu_b_L;
  if(priorPHI == "SL"){
    double nu_a_L_in = priorL_in["nu_a"];
    double nu_b_L_in = priorL_in["nu_b"];
    nu_a_L = nu_a_L_in;
    nu_b_L = nu_b_L_in;
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

  int mid = 0;
  int mim = 0;
  if(priorPersistence==true){
    mid +=draws;
    mim +=M;
  }
  mat m_i_draws(mid, mim);

  int phi_hyperparameter_size(0);
  if(priorPHI == "DL" || priorPHI == "DL_h"){
    phi_hyperparameter_size += 2. + 2*n; // a + zeta + n(theta + psi)
  }else if(priorPHI == "R2D2"){
    phi_hyperparameter_size += 2. + 2*n; // xi + zeta + n(theta + psi)
  }else if(priorPHI == "SSVS"){
    phi_hyperparameter_size += 2*n; // n(gammas + p_i)
  }else if(priorPHI == "HMP"){
    phi_hyperparameter_size += 2; // lambda_1 + lambda_2
  }else if(priorPHI == "SL"){
    phi_hyperparameter_size += n;
  }
  arma::mat phi_hyperparameter_draws(draws, phi_hyperparameter_size);

  int l_hyperparameter_size(0);
  if(priorL == "DL" || priorL == "DL_h"){
    l_hyperparameter_size += 2. + 2*n_L;
  }else if(priorL == "R2D2"){
    l_hyperparameter_size += 2. + 2*n_L; // xi + zeta + n(theta + psi)
  }else if(priorL == "SSVS"){
    l_hyperparameter_size += 2*n_L;
  }else if(priorL == "HMP"){
    l_hyperparameter_size += 1;
  }else if(priorL == "SL"){
    l_hyperparameter_size += n_L;
  }
  arma::mat l_hyperparameter_draws(draws, l_hyperparameter_size);

  //L = arma::eye(size(L));
  arma::mat eps(K,M);
  if(priorPHI == "DL" || priorPHI == "DL_h" || priorPHI == "R2D2"){
    for(int i=0; i<K; i++){
      for(int j=0; j<M; j++){
        if(PHI0(i,j)==0){
          eps(i,j)=1e-100;
        }else {
          eps(i,j)=1e-10;
        }
      }
    }
  }

  //-----------------------------------SAMPLER--------------------------------//

  const int tot = draws + burnin;
  // Initialize progressbar
  Progress p(tot, progressbar);
  Timer timer;
  timer.step("start");
  for(int rep = 0; rep < tot; rep++){

    //----1) Draw PHI (reduced form VAR coefficients)
    if(priorPHI == "SL"){
      sample_PHI_SL(PHI, PHI0, Y, X, L, d_sqrt, Gamma, K, M, nu_a, nu_b);
    }else{
    try{
      sample_PHI(PHI, PHI0, Y, X, L, d_sqrt, V_prior, K, M, false);
    } catch(...){
      ::Rf_error("Couldn't sample PHI in rep %i.", rep);
    }

    if(priorPHI == "DL" || priorPHI == "DL_h" || priorPHI == "R2D2"){
      PHI_diff = PHI - PHI0;
      // if regularization gets extreme, often there appear zeros (numerical issue)
      // coefficients must not be zero, otherwise problems with do_rgig1
      // anyhow, a realized value of a continous pd cannot be exactly zero
      for (int ii = 0; ii<K; ii++) {
        for (int jj = 0; jj<M; jj++){
          if(PHI_diff(ii,jj) == 0) {

            if(R::rbinom( 1, 0.5 )==0){
              PHI(ii,jj) = PHI0(ii,jj) + 1e-100;//eps(ii,jj);//1e-100;
              PHI_diff(ii,jj) = 1e-100;
            }else{
              PHI(ii,jj) = PHI0(ii,jj) - 1e-100;//eps(ii,jj);//1e-100;
              PHI_diff(ii,jj) = -1e-100;
            }
          }else if(PHI_diff(ii,jj) < 1e-100 && PHI_diff(ii,jj) > 0){ //eps(ii,jj)
              PHI(ii,jj) = PHI0(ii,jj) + 1e-100;//eps(ii,jj);
            PHI_diff(ii,jj) = 1e-100;
          }else if (PHI_diff(ii,jj) > -1e-100 && PHI_diff(ii,jj) < 0){ //-eps(ii,jj)
              PHI(ii,jj) = PHI0(ii,jj) - 1e-100;//eps(ii,jj);
            PHI_diff(ii,jj) = -1e-100;
          }
        }
      }
    }else{
      PHI_diff = PHI - PHI0;
      }

    //----2) Sample hyperparameters of hierarchical priors (prior variances V_i)

    if(priorPHI == "DL" || priorPHI == "DL_h"){

      if(priorPHI == "DL_h"){
        sample_DL_hyper(DL_a, theta, prep1, prep2, zeta, a_vec);
      }
      //try{
        sample_V_i_DL(V_i, PHI_diff(i_ocl), DL_a , zeta, psi, theta); //, priorPHI == "DL_h"
     // }catch(...){
      //  ::Rf_error("Couldn't sample V_i (DL prior) in run %i.",  rep);
      //}


    }else if(priorPHI == "R2D2"){


      //try{
      sample_V_i_R2D2(V_i, PHI_diff(i_ocl), api , zeta_r2d2, psi_r2d2,
                      theta_r2d2, xi, a_r2d2, b_r2d2 ); //, priorPHI == "DL_h"
      // }catch(...){
      //  ::Rf_error("Couldn't sample V_i (R2D2 prior) in run %i.",  rep);
      //}


    }else if(priorPHI == "SSVS"){

      if(rep > 0.1*burnin){
        sample_V_i_SSVS(V_i, gammas, p_i, PHI_diff(i_ocl), tau_0, tau_1, SSVS_s_a, SSVS_s_b);

      }
    }else if(priorPHI == "HMP"){

      if(rep > 0.1*burnin){
        sample_V_i_HMP(lambda_1, lambda_2, V_i_long, s_r_1(0), s_r_1(1), s_r_2(0),
                       s_r_2(1), PHI_diff, V_i_prep, n_ol, n_cl, i_ol, i_cl);
        V_i = V_i_long(i_ocl);
      }
    }

    V_i_long(i_ocl) = V_i;

    V_prior = reshape(V_i_long, K+intercept, M);

    //----(optional) Sample prior mean of first ownlags //???
    if(priorPersistence==true){
      sample_prior_mean(m_i, PHI(i_fol), V_i(i_fol), mu0, b0);
      PHI0(i_fol) = m_i;
    }
    }//end if SL

    //----3) Draw Sigma_t = inv(L)'*D_t*inv(L), where L is upper unitriangular
    //       and D diagonal
    //----3a) Draw free off-diagonal elements in L

    arma::mat resid = Y - X*PHI;

    if(priorL == "SL"){
      sample_L_SL(L, resid, d_sqrt, omega, nu_a_L, nu_b_L);
    }else{

      try{
        sample_L(L, resid, V_i_L, d_sqrt);
      }
      catch(...){
        ::Rf_error("Couldn't sample L in rep %i.", rep);
      }

      if(priorL == "DL" || priorL == "DL_h" || priorL == "R2D2"){

        for (int jj = 0; jj<M; jj++){
          for (int ii = 0; ii<jj; ii++) {
            if(L(ii,jj) == 0) {
              if(R::rbinom( 1, 0.5 )==0){
                L(ii,jj) = 1e-100;
              }else{
                L(ii,jj) = -1e-100;
              }
            }else
              if(L(ii,jj) < 1e-100 && L(ii,jj) > 0){
                L(ii,jj) = 1e-100;
              }else if (L(ii,jj) > -1e-100 && L(ii,jj) < 0){
                L(ii,jj) = -1e-100;
              }
          }
        }
      }

      //----3b) Draw hyperparameters of hierarchical priors
      l = L(L_upper_indices);
      if(priorL == "DL" || priorL == "DL_h"){

        if(priorL == "DL_h"){

          sample_DL_hyper(DL_b, theta_L, prep1_L, prep2_L, zeta_L, b_vec);
        }
        try{
          sample_V_i_DL(V_i_L, l, DL_b , zeta_L, psi_L, theta_L); //, priorL == "DL_h"
        } catch (...) {
          ::Rf_error("Couldn't sample V_i_L (DL prior)  in run %i", rep);

        }

      }else if(priorL == "R2D2"){

        sample_V_i_R2D2(V_i_L, l, api_L , zeta_L_r2d2, psi_L_r2d2,
                        theta_L_r2d2, xi_L, a_L_r2d2, b_L_r2d2 );

      }else if(priorL == "SSVS"){

        sample_V_i_SSVS(V_i_L, gammas_L, p_i_L, l, tau_0_L, tau_1_L, SSVS_s_a_L, SSVS_s_b_L);

      }else if(priorL == "HMP"){

        sample_V_i_L_HMP(lambda_3, V_i_L, s_r_3(0), s_r_3(1), l);
      }

    }//--------------end if SL

    //----3c) Draw elements of D_t
    //        in case of SV use package stochvol
    arma::mat str_resid = resid*L; // structural (orthogonalized) residuals

    if(SV == false || (priorPHI == "SL" && rep < 0.1*burnin)){ //|| (priorPHI == "SL" && rep < 0.1*burnin) //???
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

      if(priorPersistence==true){
        m_i_draws.row(rep-burnin) = m_i;
      }

      if(priorPHI == "DL" || priorPHI == "DL_h"){

        phi_hyperparameter_draws(rep-burnin, 0) = zeta;
        phi_hyperparameter_draws(rep-burnin, span(1,(n))) = trans(psi.as_col());
        phi_hyperparameter_draws(rep-burnin, span(n+1.,(2*n))) = trans(theta.as_col());
        phi_hyperparameter_draws(rep-burnin, phi_hyperparameter_size-1.) = DL_a;

      }else if(priorPHI == "R2D2"){

        phi_hyperparameter_draws(rep-burnin, 0) = zeta_r2d2 ;
        phi_hyperparameter_draws(rep-burnin, span(1,(n))) = trans(psi_r2d2.as_col());
        phi_hyperparameter_draws(rep-burnin, span(n+1.,(2*n))) = trans(theta_r2d2.as_col());
        phi_hyperparameter_draws(rep-burnin, phi_hyperparameter_size-1.) = xi ;

      }else if(priorPHI == "SSVS"){

        phi_hyperparameter_draws(rep-burnin, span(0, (n-1.))) = gammas;
        phi_hyperparameter_draws(rep-burnin, span(n, (phi_hyperparameter_size-1.))) = p_i;

      }else if(priorPHI == "HMP"){
        phi_hyperparameter_draws(rep-burnin, 0) = lambda_1;
        phi_hyperparameter_draws(rep-burnin, 1) = lambda_2;
      }else if(priorPHI == "SL"){

        phi_hyperparameter_draws(rep-burnin, span(0, (n-1.))) = vectorise(Gamma);

      }

      if(priorL == "DL" || priorL == "DL_h"){

        l_hyperparameter_draws(rep-burnin, 0) = zeta_L;
        l_hyperparameter_draws(rep-burnin, span(1,n_L)) = trans(psi_L.as_col());
        l_hyperparameter_draws(rep-burnin, span(n_L+1.,2*n_L)) = trans(theta_L.as_col());
        l_hyperparameter_draws(rep-burnin, l_hyperparameter_size-1.) = DL_b;

      }else if(priorL == "R2D2"){

        l_hyperparameter_draws(rep-burnin, 0) = zeta_L_r2d2 ;
        l_hyperparameter_draws(rep-burnin, span(1,(n_L))) = trans(psi_L_r2d2.as_col());
        l_hyperparameter_draws(rep-burnin, span(n_L+1.,(2*n_L))) = trans(theta_L_r2d2.as_col());
        l_hyperparameter_draws(rep-burnin, l_hyperparameter_size-1.) = xi_L ;

      }else if(priorL == "SSVS"){

        l_hyperparameter_draws(rep-burnin, span(0, (n_L-1.))) = gammas_L;
        l_hyperparameter_draws(rep-burnin, span(n_L, (l_hyperparameter_size-1.))) = p_i_L;

      }else if(priorL == "HMP"){

        l_hyperparameter_draws(rep-burnin, 0) = lambda_3;
      }else if(priorL == "SL"){
        l_hyperparameter_draws(rep-burnin, span(0, (n_L-1.))) = omega;
      }
    }

    p.increment();
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
    Named("i_ocl") = i_ocl, // ???
    Named("i_i") = i_i, // ???
    Named("s_r_3") = s_r_3, // ???
    Named("tau0") = tau_0, // ???
    Named("tau1") = tau_1,// ???
    Named("V_i_long") = V_i_long, // ???
    Named("V_i_L") = V_i_L,
    //Named("a_vec") = a_vec, // ???
    //Named("a_mat") = a_mat, // ???
    //Named("prep1") = prep1, // ???
    //Named("prep2") = prep2, // ???
    Named("Gamma") = Gamma, // ???
    Named("m_i") = m_i_draws,
    Named("priorPers") = priorPersistence,
    Named("PHI0")=PHI0
  );

  return out;
}
