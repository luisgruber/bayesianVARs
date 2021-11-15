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
              arma::mat PHI_prior,
              arma::mat L,
              arma::vec V_i,
              arma::vec V_i_L,
              const List sv_spec,
              arma::mat h,
              arma::mat sv_para,
              const bool progressbar){


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


  //-----------------------------------SAMPLER--------------------------------//

  const int tot = draws + burnin;
  // Initialize progressbar
  Progress p(tot, progressbar);
  Timer timer;
  timer.step("start");
  for(int rep = 0; rep < tot; rep++){

    //----1) Draw PHI (reduced form VAR coefficients)
    sample_PHI(PHI, PHI_prior,  Y,  X,L,  d, V_i, K,  M);

    //----2) Sample hyperparameter of hierarchical priors
    //----3) Draw Sigma_t = inv(L)'*D_t*inv(L), where L is upper unitriangular
    //       and D diagonal
    //----3a) Draw free off-diagonal elements in L
    arma::mat resid = Y - X*PHI;
    sample_L(L, resid, V_i_L, d);
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

    //-------Store draws after burnin
    if(rep >= burnin){

      PHI_draws.row(rep-burnin) = PHI;
      L_draws.row(rep-burnin) = L;
      sv_latent_draws.row(rep-burnin) = h;
      sv_para_draws.row((rep-burnin)) = sv_para;

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
    Named("bench") = time
  );

  return out;
}
