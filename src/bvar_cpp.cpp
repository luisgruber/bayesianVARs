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
              const arma::mat PHI0,
              const List priorPHI_in,
              const List priorL_in,
              arma::mat L,
              const List sv_spec,
              arma::mat h,
              arma::mat sv_para,
              const arma::imat i_mat,
              const arma::ivec i_vec,
              const bool progressbar){

//-------------------------Preliminaries--------------------------------------//

  const double n = K*M; // number of VAR coefficients // ??? without intercept
  arma::mat PHI_diff; // will hold PHI - PHI0
  const arma::uvec i_ol = arma::find(i_vec > 0); // indicator for ownlags
  const arma::uvec i_cl = arma::find(i_vec < 0); // indicator for crosslags
  const int n_ol = i_ol.size(); // nr. of ownlags
  const int n_cl = i_cl.size(); // nr. of crosslags
  const arma::uvec i_ocl= arma::find(i_vec != 0.); // indicator for all coefficients except intercept
  const arma::uvec i_i= arma::find(i_vec == 0.); // indicator for intercepts

//--------------------Initialization of hyperparameters-----------------------//

//---- PHI
  std::string priorPHI = priorPHI_in["prior"];
  // V_i holds prior variances (without intercepts)
  arma::vec V_i(n);

  arma::vec V_i_long(n+M*intercept); // ??? V_i plus intercept prior variances


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

  // in case of hyperprior on a (discrete uniform)
  //arma::vec a_vec(1000); // will hold grid of possible a
  //arma::rowvec prep2(1000);
  //arma::mat a_mat(1000,n); // grid as matrix
  //arma::mat prep1(n,1000);
  NumericVector a_vec_in;
  NumericVector prep2_in;
  NumericMatrix a_mat_in(1000,n);
  NumericMatrix prep1_in(n,1000);
  if(priorPHI == "DL_h"){
    a_vec_in = priorPHI_in["a_vec"];
    prep2_in = priorPHI_in["prep2"];
    a_mat_in = wrap(priorPHI_in["a_mat"]);
    prep1_in = wrap(priorPHI_in["prep1"]);
    //double dist = 0.5 - 1/static_cast<double>(n); // grid should range from 1/n to 0.5
    //double stps = dist / (1000-1); // compute stepsize

    //a_vec(0) = 1/static_cast<double>(n);
    //for(int i=1; i<1000; ++i){
    //  a_vec(i) = a_vec(i-1) + stps;
    //}
    //a_mat.each_col() = a_vec;
    // some precalculations for the conditional posterior
    // for faster evaluation of Dirichlet density
    //prep1 = (a_mat).t() - 1;
    //prep2 = vectorise(arma::lgamma(arma::sum(a_mat.t(),0)) -
    //  arma::sum(arma::lgamma(a_mat.t()),0));

  }
  arma::vec a_vec(a_vec_in.begin(), a_vec_in.length(), false);
  arma::rowvec prep2(prep2_in.begin(), prep2_in.length(), false);
  arma::mat a_mat(a_mat_in.begin(), a_mat_in.nrow(), a_mat_in.ncol(), false);
  arma::mat prep1(prep1_in.begin(), prep1_in.nrow(), prep1_in.ncol(), false);

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
    V_i_long(i_i) = priorIntercept % V_i_prep(i_i);
    V_i_long(i_ol) = lambda_1*V_i_prep(i_ol);
    V_i_long(i_cl) = lambda_2*V_i_prep(i_cl);
    V_i = V_i_long(i_ocl);
  }
  arma::vec s_r_1(s_r_1_in.begin(), s_r_1_in.length(), false);
  arma::vec s_r_2(s_r_2_in.begin(), s_r_2_in.length(), false);


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

  // in case of hyperprior on b
  //arma::vec b_vec(1000);
  //arma::rowvec prep2_L(1000);
  //arma::mat b_mat(1000,n_L);
  //arma::mat prep1_L(n_L,1000);
  NumericVector b_vec_in;
  NumericVector prep2_L_in;
  NumericMatrix b_mat_in(1000,n_L);
  NumericMatrix prep1_L_in(n_L,1000);
  if(priorL == "DL_h"){
    b_vec_in = priorL_in["b_vec"];
    prep2_L_in = priorL_in["prep2"];
    b_mat_in = wrap(priorL_in["b_mat"]);
    prep1_L_in = wrap(priorL_in["prep1"]);
   // double dist0 = 0.5 - 1/static_cast<double>(n_L); // grid should range from 1/n to 0.5
   // double stps0 = dist0 / (1000-1); // compute stepsize

   // b_vec(0) = 1/static_cast<double>(n_L);
  //  for(int i=1; i<1000; ++i){
   //   b_vec(i) = b_vec(i-1) + stps0;
   // }
   // b_mat.each_col() = b_vec;
    // some precalculations for the conditional posterior
    // for faster evaluation of Dirichlet density
   // prep1_L = (b_mat).t() - 1;
   // prep2_L = vectorise(arma::lgamma(arma::sum(b_mat.t(),0)) -
   //   arma::sum(arma::lgamma(b_mat.t()),0));

  }
  arma::vec b_vec(b_vec_in.begin(), b_vec_in.length(), false);
  arma::rowvec prep2_L(prep2_L_in.begin(), prep2_L_in.length(), false);
  arma::mat b_mat(b_mat_in.begin(), b_mat_in.nrow(), b_mat_in.ncol(), false);
  arma::mat prep1_L(prep1_L_in.begin(), prep1_L_in.nrow(), prep1_L_in.ncol(), false);

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

  V_i_long(i_ocl) = V_i; // ???
  V_i_long(i_i) = priorIntercept;
  //arma::mat V_prior = arma::reshape(V_i, K, M);
  arma::mat V_prior = arma::reshape(V_i_long, K + intercept, M);

  //-----------------------------SV-settings----------------------------------//
  // Import sv_spec
  NumericVector sv_priormu = sv_spec["priormu"] ;
  NumericVector sv_priorphi= sv_spec["priorphi"];
  NumericVector sv_priorsigma2 = sv_spec["priorsigma2"] ;

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

  int phi_hyperparameter_size(0);
  if(priorPHI == "DL" || priorPHI == "DL_h"){
    phi_hyperparameter_size += 2. + 2*n; // a + zeta + n(theta + psi)
  }else if(priorPHI == "SSVS"){
    phi_hyperparameter_size += 2*n; // n(gammas + p_i)
  }else if(priorPHI == "HMP"){
    phi_hyperparameter_size += 2; // lambda_1 + lambda_2
  }
  arma::mat phi_hyperparameter_draws(draws, phi_hyperparameter_size);

  int l_hyperparameter_size(0);
  if(priorL == "DL" || priorL == "DL_h"){
    l_hyperparameter_size += 2. + 2*n_L;
  }else if(priorL == "SSVS"){
    l_hyperparameter_size += 2*n_L;
  }else if(priorL == "HMP"){
    l_hyperparameter_size += 1;
  }
  arma::mat l_hyperparameter_draws(draws, l_hyperparameter_size);

  //L = arma::eye(size(L));

  //-----------------------------------SAMPLER--------------------------------//

  const int tot = draws + burnin;
  // Initialize progressbar
  Progress p(tot, progressbar);
  Timer timer;
  timer.step("start");
  for(int rep = 0; rep < tot; rep++){

    //----1) Draw PHI (reduced form VAR coefficients)
    try{
      sample_PHI(PHI, PHI0, Y, X, L, d_sqrt, V_prior, K, M, false);
    } catch(...){
      ::Rf_error("Couldn't sample PHI in rep %i.", rep);
    }

    if(priorPHI == "DL" || priorPHI == "DL_h"){
      // if regularization gets extreme, often there appear zeros (numerical issue)
      // coefficients must not be zero, otherwise problems with do_rgig1
      // anyhow, a realized value of a continous pd cannot be exactly zero
      for (int ii = 0; ii<K; ii++) {
        for (int jj = 0; jj<M; jj++){
          if(PHI(ii,jj) == 0) {
            if(R::rbinom( 1, 0.5 )==0){
              PHI(ii,jj) = 1e-100;
            }else{
              PHI(ii,jj) = -1e-100;
            }
          }else
            if(PHI(ii,jj) < 1e-100 && PHI(ii,jj) > 0){
              PHI(ii,jj) = 1e-100;
            }else if (PHI(ii,jj) > -1e-100 && PHI(ii,jj) < 0){
              PHI(ii,jj) = -1e-100;
            }
        }
      }
    }

    PHI_diff = PHI - PHI0;

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


    }else if(priorPHI == "SSVS"){

      if(rep > 0.1*burnin){
        sample_V_i_SSVS(V_i, gammas, p_i, PHI_diff(i_ocl), tau_0, tau_1, SSVS_s_a, SSVS_s_b);

      }
    }else if(priorPHI == "HMP"){

      if(rep > 0.1*burnin){
        sample_V_i_HMP(lambda_1, lambda_2, V_i, s_r_1(0), s_r_1(1), s_r_2(0),
                       s_r_2(1), PHI_diff(i_ocl), V_i_prep, n_ol, n_cl, i_ol, i_cl);
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

    if(priorL == "DL" || priorL == "DL_h"){

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

    //----3b) Draw elements of D_t
    //        in case of SV use package stochvol
    arma::mat str_resid = resid*L; // structural (orthogonalized) residuals

    const arma::mat resid_norm = log(square(str_resid)); // + offset??
    for(int j=0; j < M; j++){
      arma::vec h_j  = h.unsafe_col(j);  // unsafe_col reuses memory, h will automatically be overwritten
      arma::uvec mixind_j = mixind.unsafe_col(j);
      double mu = sv_para(0,j),
        phi = sv_para(1,j),
        sigma = sv_para(2,j),
        h0_j = sv_para(3,j);
      stochvol::update_fast_sv(resid_norm.col(j), mu, phi, sigma, h0_j, h_j, mixind_j, prior_spec, expert_sv);
      sv_para.col(j) = arma::colvec({mu, phi, sigma, h0_j});
    }
    d_sqrt = exp(h/2);

    //----4) Draw hyperparameters of hierarchical priors
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

    }else if(priorL == "SSVS"){

      sample_V_i_SSVS(V_i_L, gammas_L, p_i_L, l, tau_0_L, tau_1_L, SSVS_s_a_L, SSVS_s_b_L);

    }else if(priorL == "HMP"){

      sample_V_i_L_HMP(lambda_3, V_i_L, s_r_3(0), s_r_3(1), l);
    }

    //-------Store draws after burnin
    if(rep >= burnin){

      PHI_draws.row(rep-burnin) = PHI;
      L_draws.row(rep-burnin) = L;
      sv_latent_draws.row(rep-burnin) = h;
      sv_para_draws.row((rep-burnin)) = sv_para;

      if(priorPHI == "DL" || priorPHI == "DL_h"){

        phi_hyperparameter_draws(rep-burnin, 0) = zeta;
        phi_hyperparameter_draws(rep-burnin, span(1,(n))) = trans(psi.as_col());
        phi_hyperparameter_draws(rep-burnin, span(n+1.,(2*n))) = trans(theta.as_col());
        phi_hyperparameter_draws(rep-burnin, phi_hyperparameter_size-1.) = DL_a;

      }else if(priorPHI == "SSVS"){

        phi_hyperparameter_draws(rep-burnin, span(0, (n-1.))) = gammas;
        phi_hyperparameter_draws(rep-burnin, span(n, (phi_hyperparameter_size-1.))) = p_i;

      }else if(priorPHI == "HMP"){
        phi_hyperparameter_draws(rep-burnin, 0) = lambda_1;
        phi_hyperparameter_draws(rep-burnin, 1) = lambda_2;
      }

      if(priorL == "DL" || priorL == "DL_h"){

        l_hyperparameter_draws(rep-burnin, 0) = zeta_L;
        l_hyperparameter_draws(rep-burnin, span(1,n_L)) = trans(psi_L.as_col());
        l_hyperparameter_draws(rep-burnin, span(n_L+1.,2*n_L)) = trans(theta_L.as_col());
        l_hyperparameter_draws(rep-burnin, l_hyperparameter_size-1.) = DL_b;

      }else if(priorL == "SSVS"){

        l_hyperparameter_draws(rep-burnin, span(0, (n_L-1.))) = gammas_L;
        l_hyperparameter_draws(rep-burnin, span(n_L, (l_hyperparameter_size-1.))) = p_i_L;

      }else if(priorL == "HMP"){

        l_hyperparameter_draws(rep-burnin, 0) = lambda_3;
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
    Named("a_vec") = a_vec,
    Named("a_mat") = a_mat,
    Named("prep1") = prep1,
    Named("prep2") = prep2
  );

  return out;
}
