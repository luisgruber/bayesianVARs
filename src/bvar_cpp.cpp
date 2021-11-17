#include <RcppArmadillo.h>
#include <stochvol.h>
#include <progress.hpp>
#include <Rcpp/Benchmark/Timer.h>
#include "sample_coefficients.h"

using namespace Rcpp;

//[[Rcpp::export]]
List bvar_cpp(const arma::mat Y,
              const arma::mat X,
              const int M,
              const int T,
              const int K,
              const int draws,
              const int burnin,
              arma::mat PHI,
              const arma::mat PHI0,
              const std::string priorPHI,
              double DL_a, // no const because of DL_hyper
              const std::string priorL,
              const double DL_b,
              arma::mat L,
              arma::vec V_i,
              arma::vec V_i_L,
              const List sv_spec,
              arma::mat h,
              arma::mat sv_para,
              const bool progressbar){

  //-----------------------Preliminaries--------------------------------------//
  const double n = K*M; // number of VAR coefficients
  arma::mat PHI_diff; // will hold PHI - PHI0
  arma::mat V_prior = arma::reshape(V_i, K, M);

  //----------------------Initialization of hyperpriors-----------------------//
  // DL prior on PHI
  arma::vec psi(n, arma::fill::value(1.0));
  double zeta=10;
  arma::vec theta(n, arma::fill::value(1.0));
  // with hyperprior on a
  arma::vec a_vec(1000);
  arma::rowvec prep2(1000);
  arma::mat a_mat(1000,n);
  arma::mat prep1(n,1000);
  if(priorPHI == "DL_h"){
    double dist = 0.5 - 1/static_cast<double>(n); // grid should range from 1/n to 0.5
    double stps = dist / (1000-1); // compute stepsize

    a_vec(0) = 1/static_cast<double>(n);
    for(int i=1; i<1000; ++i){
      a_vec(i) = a_vec(i-1) + stps;
    }
    a_mat.each_col() = a_vec;
    // some precalculations for the conditional posterior
    prep1 = (a_mat).t() - 1;
    prep2 = vectorise(arma::lgamma(arma::sum(a_mat.t(),0)) -
      arma::sum(arma::lgamma(a_mat.t()),0));

  }

  // DL prior on PHI
  uvec L_upper_indices = trimatu_ind( size(L),  1);
  arma::vec l = L(L_upper_indices);
  const double n_L = l.size();
  arma::vec psi_L(n_L, arma::fill::value(1.0));
  double zeta_L=10;
  arma::vec theta_L(n_L, arma::fill::value(1.0));

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
  arma::umat mixind(T,M, arma::fill::value(5));
  //initialize d
  arma::mat d = arma::exp(h);

  //------------------------------------STORAGE-------------------------------//

  arma::cube sv_latent_draws(draws, T, M);
  arma::cube sv_para_draws(draws, 4,M);
  arma::cube PHI_draws(draws, K, M);
  arma::cube L_draws(draws, M, M);
  int hyperparameter_size;
  if(priorPHI == "DL" || priorPHI == "DL_h"){
    hyperparameter_size = 2 + 2*n;
  }
  if(priorL == "DL"){
    hyperparameter_size += 1 + 2*n_L;
  }
  arma::mat hyperparameter_draws(draws, hyperparameter_size);


  //-----------------------------------SAMPLER--------------------------------//

  const int tot = draws + burnin;
  // Initialize progressbar
  Progress p(tot, progressbar);
  Timer timer;
  timer.step("start");
  for(int rep = 0; rep < tot; rep++){

    //----1) Draw PHI (reduced form VAR coefficients)

    sample_PHI(PHI, PHI0, Y, X, L, d, V_prior, K, M);

    if(priorPHI == "DL"){
      // coefficients must not be zero, otherwise problems with do_rgig1
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

    //----2) Sample hyperparameter of hierarchical priors

    if(priorPHI == "DL" || priorPHI == "DL_h"){

     if(priorPHI == "DL_h"){
       sample_DL_hyper(DL_a, theta, prep1, prep2, zeta, a_vec);
     }
     sample_V_i_DL(V_i, PHI_diff, DL_a , zeta, psi, theta);
    }
    V_prior = reshape(V_i, K, M);

    //----3) Draw Sigma_t = inv(L)'*D_t*inv(L), where L is upper unitriangular
    //       and D diagonal
    //----3a) Draw free off-diagonal elements in L


    arma::mat resid = Y - X*PHI;
    sample_L(L, resid, V_i_L, d);

    if(priorL == "DL"){

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
    d = exp(h);

    //----4) Draw hyperparameters of hierarchical priors
    if(priorL == "DL"){
      l = L(L_upper_indices);
      sample_V_i_DL(V_i_L, l, DL_b , zeta_L, psi_L, theta_L);
    }

    //-------Store draws after burnin
    if(rep >= burnin){

      PHI_draws.row(rep-burnin) = PHI;
      L_draws.row(rep-burnin) = L;
      sv_latent_draws.row(rep-burnin) = h;
      sv_para_draws.row((rep-burnin)) = sv_para;
      hyperparameter_draws(rep-burnin, 0) = zeta;
      hyperparameter_draws(rep-burnin, span(1,(n))) = trans(psi.as_col());
      hyperparameter_draws(rep-burnin, span(n+1,(2*n))) = trans(theta.as_col());
      hyperparameter_draws(rep-burnin, 2*n+1) = zeta_L;
      hyperparameter_draws(rep-burnin, span(2*n+2,2*n+1+n_L)) = trans(psi_L.as_col());
      hyperparameter_draws(rep-burnin, span(2*n+2+n_L,hyperparameter_size-2)) = trans(theta_L.as_col());
      hyperparameter_draws(rep-burnin, hyperparameter_size-1) = DL_a;
    }

    p.increment();
  }
  timer.step("end");
  NumericVector time(timer);

  List out = List::create(
    Named("PHI_draws") = PHI_draws,
    Named("L_draws") = L_draws,
    Named("sv_latent") = sv_latent_draws,
    Named("sv_para_draws") = sv_para_draws,
    Named("hyperparameter_draws") = hyperparameter_draws,
    Named("bench") = time
  );

  return out;
}
