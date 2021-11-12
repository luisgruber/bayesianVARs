#include <RcppArmadillo.h>
#include <stochvol.h>

using namespace Rcpp;

//[[Rcpp::export]]
List bvar_cpp(arma::mat str_resid, //structural residuals
               const int draws,
               const int burnin,
               const int M,
               const int T,
               const List sv_spec,
               arma::mat h,
               arma::mat sv_para){

  //-----------------------------SV-settings-----------------------------------
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

  //expert settings: these are the settings as the default settings for the idiosyncratic variances from package factorstochvol (https://github.com/gregorkastner/factorstochvol, accessed 2021-11-12)
  using stochvol::ExpertSpec_FastSV;
  const ExpertSpec_FastSV expert_sv {
    true,  // interweave
    stochvol::Parameterization::CENTERED,  // centered baseline always
    1e-8,  // B011inv,
    1e-12,  //B022inv,
    2,  // MHsteps,
    ExpertSpec_FastSV::ProposalSigma2::INDEPENDENCE,  // independece proposal for sigma
    -1,  // MHcontorl unused for independence prior for sigma
    ExpertSpec_FastSV::ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL  // immediately reject (mu,phi,sigma) if proposed phi is outside (-1, 1)
  };

  // initialize mixing indicators
  arma::umat mixind(T,M, arma::fill::value(5));

  //------------------------------------STORAGE-------------------------------//

  arma::cube sv_latent_draws(draws, T, M);
  arma::cube sv_para_draws(draws, 4,M);


  //-----------------------------------SAMPLER--------------------------------//
  const int tot = draws + burnin;
  for(int rep = 0; rep < tot; rep++){

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

    if(rep >= burnin){

      sv_latent_draws.row(rep-burnin) = h;
      sv_para_draws.row((rep-burnin)) = sv_para;

    }

  }

  List out = List::create(
    Named("h") = sv_latent_draws,
    Named("para") = sv_para_draws
  );

  return out;
}
