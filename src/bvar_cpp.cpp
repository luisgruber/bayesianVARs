#include <RcppArmadillo.h>
#include <factorstochvol.h>
#include <stochvol.h>
#include "bvar_single_update.h"
#include "utilities_cpp.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List bvar_cpp(const arma::mat& Y,
              const arma::mat& X,
              const int& M,
              const int& T,
              const int& K,
              const int& draws,
              const int& burnin,
              const int& thin,
              const std::string& tvp_keep,
              const bool& intercept,
              const arma::rowvec& priorIntercept,
              arma::mat& PHI_priormean,
              const List priorPHI_in,
              const List priorSigma_in,
              const List Rstartvals_in,
              const arma::imat& i_mat,
              const bool& progressbar,
              const bool& huge,
              const double& abs_floor){
  
  // Rcpp objects reuse memory
  // arma make copy by default unless deactivated
  
  //--------------------Initialization -----------------------//
  
  const bool stable = priorPHI_in["stable"];
  const double stable_lim = priorPHI_in["stable_lim"];
  
  List startvals_in(clone(Rstartvals_in)); // clone makes copy
  //---- PHI
  
  Rcpp::List phi_startvals = startvals_in["phi"];
  NumericMatrix PHI_in = phi_startvals["PHI"];
  arma::mat PHI(PHI_in.begin(), PHI_in.nrow(), PHI_in.ncol(), false);
  std::string priorPHI = priorPHI_in["prior"];
  
  if(stable) {
    if(!is_stable(PHI, M, stable_lim)) {
      Rcpp::stop("Starting value for PHI does not fulfill stable criterion, although prior on PHI is specified to do so. See, ?specify_prior_phi and ?get_startvals!");
    }
  }
  
  // PHI_priorvariance holds prior variances of coefficicients without intercept
  // priorIntercept holds fixed prior variances of intercept coefficients
  arma::mat PHI_priorvariance(K, M);
  
  if(priorPHI == "normal"){
    NumericMatrix V_i_in = priorPHI_in["V_i"];
    arma::mat armaV_i_in(V_i_in.begin(), V_i_in.nrow(), V_i_in.ncol(), false);
    //in case of 'normal' PHI_priorvariance is fixed at user specified values
    PHI_priorvariance = armaV_i_in;
  }
  
  //  sub-groups for semi-global-local GL priors and SSVS prior
  const int n_groups = priorPHI_in["n_groups"]; // number of distinct groups
  IntegerVector groups_in = priorPHI_in["groups"];
  arma::ivec groups(groups_in.begin(), groups_in.size(),false);
  
  // Hyperparameters of DL, GT (NG or R2D2), HS and SSVS
  // lambda is (primary) local shrinkage scale parameter, in case of SSVS it is 0-1 indicator
  NumericMatrix lambda_in = phi_startvals["lambda"];
  arma::mat lambda(lambda_in.begin(), lambda_in.nrow(), lambda_in.ncol(), false);
  // psi is auxiliary local shrinkage scale
  NumericMatrix psi_in = phi_startvals["psi"];
  arma::mat psi(psi_in.begin(), psi_in.nrow(), psi_in.ncol(), false);
  // xi is (semi) global shrinkage scale parameter
  NumericVector xi_in = phi_startvals["xi"];
  arma::vec xi(xi_in.begin(), xi_in.length(), false);
  NumericVector a_in = phi_startvals["a"];
  // a is either (semi) global concentration parameter that determines degree of shrinkage, or
  // in case of SSVS prior inclusion probability
  arma::vec a(a_in.begin(), a_in.length(), false);
  //a_hyper is true if a discrete prior is placed on a
  const bool a_hyper = priorPHI_in["a_hyper"];
  // b is shape parameter of InverseGamma prior on (semi global) local scales (NG and R2D2 only)
  NumericVector b_in = priorPHI_in["b"];
  arma::vec b(b_in.begin(), b_in.length(), false);
  // c is scale parameter of InverseGamma prior on (semi global) local scales (NG only, for R2D2 c=0.5*a)
  NumericVector c_in = priorPHI_in["c"];
  arma::vec c(c_in.begin(), c_in.length(), false);
  
  // a_vec holds support points of discrete prior on a
  NumericVector a_vec_in = priorPHI_in["a_vec"]; 
  arma::vec a_vec(a_vec_in.begin(), a_vec_in.size(), false);
  // a_weight holds prior pobabilities of discrete prior
  NumericVector a_weight_in = priorPHI_in["a_weight"]; 
  arma::vec a_weight(a_weight_in.begin(), a_weight_in.size(), false);
  // norm_consts holds pre-computed evaluations of the Gamma function 
  // (needed for the computation of full condtional posterior of discrete prior)
  NumericVector norm_consts_in = priorPHI_in["norm_consts"];
  arma::vec norm_consts(norm_consts_in.begin(), norm_consts_in.size(), false);
  // c_vec holds support points of c if c is specified as a proportion of a
  NumericVector c_vec_in = priorPHI_in["c_vec"];
  arma::vec c_vec(c_vec_in.begin(), c_vec_in.size(), false); 
  // c_rel_a is true if c is specified as a proportion of a, e.g. R2D2 assumes c=0.5a
  const bool c_rel_a = priorPHI_in["c_rel_a"];
  
  //---- DL prior specific parameters
  // DL_plus (hidden functionality) is true, if DL prior is extended to a full (semi) global local prior
  const bool DL_plus = priorPHI_in["DL_plus"];
  if(priorPHI == "DL" ){
    PHI_priorvariance = psi % lambda % lambda;
  }
  
  //----GT prior specific parameters  (GT refers to gamma type priors a la NG or R2D2)
  // GT_vs is a scaling parameter of prior variance (for NG GT_vs=1; for R2D2 GT_vs = 0.5)
  const double GT_vs = priorPHI_in["GT_vs"];
  // for NG GT_priorkernel = "normal", for R2D2 GT_priorkernel = "exponential"
  const std::string GT_priorkernel = priorPHI_in["GT_priorkernel"];
  if(priorPHI == "GT"){
    PHI_priorvariance = lambda*GT_vs;
    if(GT_priorkernel == "exponential"){
      PHI_priorvariance %= psi;
    }
  }
  
  if(c_rel_a && a_hyper) {
    for(int g=0; g<n_groups; ++g){
      arma::vec tmp = c_vec.elem(arma::find(a_vec == a(g)));
      c(g) = tmp(0);
    }
  }
  
  //---- Horseshoe specifics
  // varpi is (semi) global auxiliary scale
  NumericVector varpi_in = phi_startvals["varpi"];
  arma::vec varpi(varpi_in.begin(), varpi_in.length(), false);
  if(priorPHI == "HS"){
    arma::ivec::iterator g;
    int j=0;
    for(g=groups.begin(); g!=groups.end(); ++g){
      arma::uvec ind = arma::find(i_mat == *g);
      PHI_priorvariance(ind) = lambda(ind) * xi(j);
      j++;
    }
  }
  
  //---- SSVS specifics
  Rcpp::NumericMatrix tau_0_rcpp = priorPHI_in["SSVS_tau0"];
  arma::mat tau_0(tau_0_rcpp.begin(), tau_0_rcpp.nrow(), tau_0_rcpp.ncol(), false);
  Rcpp::NumericMatrix tau_1_rcpp = priorPHI_in["SSVS_tau1"];
  arma::mat tau_1(tau_1_rcpp.begin(), tau_1_rcpp.nrow(), tau_1_rcpp.ncol(), false);
  double SSVS_s_a = priorPHI_in["SSVS_s_a"];
  double SSVS_s_b = priorPHI_in["SSVS_s_b"];
  if(priorPHI == "SSVS"){
    PHI_priorvariance = arma::square(tau_1%lambda + (lambda-1)%tau_0);
  }
  
  //---- HMP specifics
  // indicators for cross-/own-lags 
  const arma::uvec i_ol = arma::find(i_mat > 0);
  const arma::uvec i_cl = arma::find(i_mat < 0);
  const int n_ol = i_ol.size(); // nr. of ownlags
  const int n_cl = i_cl.size(); // nr. of crosslags
  double lambda_1= phi_startvals["lambda_1"];//0.04;
  double lambda_2= phi_startvals["lambda_2"];//0.0016;
  arma::mat V_i_prep = priorPHI_in["V_i_prep"];
  arma::vec s_r_1 = priorPHI_in["lambda_1"];
  arma::vec s_r_2 = priorPHI_in["lambda_2"];
  if(priorPHI == "HMP"){
    PHI_priorvariance(i_ol) = lambda_1*V_i_prep(i_ol);
    PHI_priorvariance(i_cl) = lambda_2*V_i_prep(i_cl);
  }
  
  // --------------------- Sigma_t
  const std::string sigma_type = priorSigma_in["type"];
  List sigma_startvals = startvals_in["sigma"];
  
  //-------- U
  NumericMatrix U_in = sigma_startvals["U"];
  arma::mat U(U_in.begin(), U_in.nrow(), U_in.ncol(), false);
  
  std::string priorU = priorSigma_in["cholesky_U_prior"];
  
  arma::uvec U_upper_indices(1);
  if(sigma_type == "cholesky"){
    U_upper_indices = trimatu_ind( size(U),  1);
  }
  arma::vec u = U(U_upper_indices);
  const arma::uword n_U = u.n_elem;
  arma::vec U_priorvariance(n_U);
  
  if(priorU == "normal"){
    Rcpp::NumericVector V_i_U_in_rcpp = priorSigma_in["cholesky_V_i"];
    arma::vec V_i_U_in(V_i_U_in_rcpp.begin(), V_i_U_in_rcpp.length(), false);
    U_priorvariance = V_i_U_in;
  }
  
  // DL, GT (NG & R2D2) and SSVS
  NumericVector lambda_U_in = sigma_startvals["lambda_U"];
  arma::vec lambda_U(lambda_U_in.begin(), lambda_U_in.length(), false);
  NumericVector psi_U_in = sigma_startvals["psi_U"];
  arma::vec psi_U(psi_U_in.begin(), psi_U_in.length(), false);
  double xi_U = sigma_startvals["xi_U"];
  double a_U = sigma_startvals["a_U"];
  const double b_U = priorSigma_in["cholesky_b"];
  double c_U = priorSigma_in["cholesky_c"];
  
  Rcpp::NumericVector a_vec_U_rcpp = priorSigma_in["cholesky_a_vec"];
  arma::vec a_vec_U(a_vec_U_rcpp.begin(), a_vec_U_rcpp.length(), false);
  Rcpp::NumericVector a_weight_U_rcpp = priorSigma_in["cholesky_a_weight"];
  arma::vec a_weight_U(a_weight_U_rcpp.begin(), a_weight_U_rcpp.length(), false);
  Rcpp::NumericVector norm_consts_U_rcpp = priorSigma_in["cholesky_norm_consts"];
  arma::vec norm_consts_U(norm_consts_U_rcpp.begin(), norm_consts_U_rcpp.length(), false);
  Rcpp::NumericVector c_vec_U_rcpp = priorSigma_in["cholesky_c_vec"];
  arma::vec c_vec_U(c_vec_U_rcpp.begin(), c_vec_U_rcpp.length(), false);
  const bool c_rel_a_U = priorSigma_in["cholesky_c_rel_a"];
  const bool a_hyper_U = priorSigma_in["cholesky_a_hyper"];
  
  //---- DL specifics
  const bool DL_plus_U = priorSigma_in["cholesky_DL_plus"];
  if(priorU == "DL"){
    U_priorvariance= psi_U % lambda_U % lambda_U;
  }
  
  //----GT specifics
  const double GT_vs_U = priorSigma_in["cholesky_GT_vs"];
  const std::string GT_priorkernel_U = priorSigma_in["cholesky_GT_priorkernel"];
  if(priorU == "GT"){
    U_priorvariance = lambda_U*GT_vs_U;
    if(GT_priorkernel_U == "exponential"){
      U_priorvariance %= psi_U;
    }
  }
  if(c_rel_a_U && a_hyper_U) {
    arma::vec tmp = c_vec_U.elem(arma::find(a_vec_U == a_U));
    c_U = tmp(0);
  }
  
  //---- Horseshoe specifics
  double varpi_U = sigma_startvals["varpi_U"];
  if(priorU == "HS"){
    U_priorvariance = lambda_U*xi_U;
  }
  
  //---- SSVS specifics
  Rcpp::NumericVector tau_0_U_rcpp = priorSigma_in["cholesky_SSVS_tau0"];
  arma::vec tau_0_U(tau_0_U_rcpp.begin(), tau_0_U_rcpp.length(), false);
  Rcpp::NumericVector tau_1_U_rcpp = priorSigma_in["cholesky_SSVS_tau1"];
  arma::vec tau_1_U(tau_1_U_rcpp.begin(), tau_1_U_rcpp.length(), false);
  double SSVS_s_a_U = priorSigma_in["cholesky_SSVS_s_a"];
  double SSVS_s_b_U = priorSigma_in["cholesky_SSVS_s_b"];
  if(priorU == "SSVS"){
    U_priorvariance = arma::square(tau_1_U%lambda_U + (lambda_U-1)%tau_0_U);
  }
  
  //---- HMP specifics
  double lambda_3 = sigma_startvals["lambda_3"];//0.001;
  NumericVector s_r_3 = priorSigma_in["cholesky_lambda_3"];
  if(priorU == "HMP"){
    U_priorvariance.fill(lambda_3);
  }
  //-----------------------factorstochvol----------------------//
  const int factors = priorSigma_in["factor_factors"];
  if(sigma_type == "cholesky" && factors != 0){
    Rcpp::stop("cholesky with factors,...,");
  }
  
  NumericMatrix facload = sigma_startvals["facload"];
  arma::mat armafacload(facload.begin(), facload.nrow(), facload.ncol(), false);
  NumericMatrix fac = sigma_startvals["fac"];
  arma::mat armafac(fac.begin(), fac.nrow(), fac.ncol(), false);
  NumericMatrix tau2 = sigma_startvals["tau2"];
  arma::mat armatau2(tau2.begin(), tau2.nrow(), tau2.ncol(), false);
  
  //----------  prior settings
  // shrinkage hyperparamters for factor loadings
  const List factor_shrinkagepriors = priorSigma_in["factor_shrinkagepriors"];
  const NumericVector aShrink = factor_shrinkagepriors["a"]; // a is hyperparameter of local prior
  const NumericVector cShrink = factor_shrinkagepriors["c"]; // c and d are hyperparameters of global prior
  const NumericVector dShrink = factor_shrinkagepriors["d"];
  const bool ngprior = priorSigma_in["factor_ngprior"];
  const bool columnwise = priorSigma_in["factor_columnwise"];
  arma::vec armalambda2((ngprior && columnwise) ? factors : ngprior ? M : 0 ); // conditional statement: if ngprior && columnwise, then factors, if ngprior (not columnwise), then M, else 0
  // underbound of absolute value of factorloadings for better preventing numerical issues
  const double factor_facloadtol = priorSigma_in["factor_facloadtol"];
  
  /* Could there be any situation where I need an offset? */
  const double factor_offset = 0; // unify with cholesky_sv ???
  
  // restriction on factor loadings
  IntegerMatrix factor_restriction = priorSigma_in["factor_restrinv"];
  arma::imat factor_armarestr(factor_restriction.begin(), factor_restriction.nrow(), factor_restriction.ncol(), false);
  arma::uvec armafacloadtunrestrictedelements = arma::find(factor_armarestr.t() != 0);
  arma::irowvec nonzerospercol = arma::sum(factor_armarestr, 0);
  arma::icolvec nonzerosperrow = arma::sum(factor_armarestr, 1);
  for (unsigned int i = 0; i < armafacload.n_rows; i++) {
    for (unsigned int j = 0; j < armafacload.n_cols; j++) {
      if (factor_armarestr(i, j) == 0) {
        armafacload(i,j) = 0.;
        armatau2(i,j) = 0.;
      }
    }
  }
  
  //-----------------------------SV-settings----------------------------------//
  // vector of length M+factors indicating whether time-varying or constant variance should be estimated
  NumericVector heteroscedastic = priorSigma_in["sv_heteroscedastic"];
  // in case of constant idiosyncratic variance, homoscedastic collects shape and scale parameters of inverse-gamma prior of idi variances
  NumericMatrix factor_homoskedastic = priorSigma_in["factor_priorhomoskedastic"]; // unify with cholesky_sv ???
  const int factor_interweaving = priorSigma_in["factor_interweaving"]; 
  
  NumericMatrix logvar_in = sigma_startvals["sv_logvar"];
  arma::mat logvar(logvar_in.begin(), logvar_in.nrow(), logvar_in.ncol(), false);
  NumericMatrix sv_para_in = sigma_startvals["sv_para"];
  arma::mat sv_para(sv_para_in.begin(), sv_para_in.nrow(), sv_para_in.ncol(), false);
  NumericVector logvar0_in = sigma_startvals["sv_logvar0"];
  arma::vec logvar0(logvar0_in.begin(), logvar0_in.length(), false);
  
  NumericMatrix priorHomoscedastic = priorSigma_in["cholesky_priorhomoscedastic"];
  
  // Import sv_spec
  NumericVector sv_priormu = priorSigma_in["sv_priormu"];
  NumericVector sv_priorphi= priorSigma_in["sv_priorphi"];
  NumericMatrix sv_priorsigma2 = priorSigma_in["sv_priorsigma2"];
  NumericVector sv_priorh0 = priorSigma_in["sv_priorh0"];
  NumericVector sv_offset = priorSigma_in["cholesky_sv_offset"];
  
  std::vector<stochvol::PriorSpec> prior_specs(M+factors);
  {
    //using stochvol::PriorSpec;
    for (int j = 0; j < M; j++) {
      if(heteroscedastic(j) == true){
        prior_specs[j] = {
          (sv_priorh0(j) <= 0) ? stochvol::PriorSpec::Latent0() : stochvol::PriorSpec::Latent0(stochvol::PriorSpec::Constant(sv_priorh0(j))), // ? : conditional statement, similar to if else
                          stochvol::PriorSpec::Mu(stochvol::PriorSpec::Normal(sv_priormu[0], sv_priormu[1])),
                          stochvol::PriorSpec::Phi(stochvol::PriorSpec::Beta(sv_priorphi[0], sv_priorphi[1])),
                          stochvol::PriorSpec::Sigma2(stochvol::PriorSpec::Gamma(sv_priorsigma2(j,0), sv_priorsigma2(j,1)))//,
          // stochvol would allow for more:
          //PriorSpec::Nu(PriorSpec::Infinity()), Nu: conditional t innovations with nu degrees of freedom  (Infinity indicates conditional standard normal innovations, the default)
          //PriorSpec::Rho(PriorSpec::Constant(0))//, Rho: leverage
          //PriorSpec::Covariates(PriorSpec::MultivariateNormal{{priorbeta[0]}, {std::pow(priorbeta[1], -2)}}) not needed, stochvol could model a constant mean
        };
      }
    }
    for (int j = M; j < (M+factors); j++) {
      if(heteroscedastic(j) == true){
        prior_specs[j] = {
          (sv_priorh0(j) <= 0) ? stochvol::PriorSpec::Latent0() : stochvol::PriorSpec::Latent0(stochvol::PriorSpec::Constant(sv_priorh0(j))),
                          stochvol::PriorSpec::Mu(stochvol::PriorSpec::Constant(0)),
                          stochvol::PriorSpec::Phi(stochvol::PriorSpec::Beta(sv_priorphi(2), sv_priorphi(3))),
                          stochvol::PriorSpec::Sigma2(stochvol::PriorSpec::Gamma(sv_priorsigma2(j,0), sv_priorsigma2(j,1)))
        };
      }
    }
  }
  
  //expert settings: these are the same settings as the default settings for the
  //idiosyncratic variances from package factorstochvol
  //(https://github.com/gregorkastner/factorstochvol, accessed 2021-11-12)
  const double B011inv = 1e-8;
  const double B022inv = 1e-12;
  const stochvol::ExpertSpec_FastSV expert_sv { // used for cholesky-sv and idiosyncratic variances of fsv
    true,  // interweave
    stochvol::Parameterization::CENTERED,  // centered baseline always
    B011inv,  // B011inv,
    B022inv,  //B022inv,
    2,  // MHsteps,
    stochvol::ExpertSpec_FastSV::ProposalSigma2::INDEPENDENCE,  // independece proposal for sigma
    -1,  // MHcontrol unused for independence prior for sigma
    stochvol::ExpertSpec_FastSV::ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL  // immediately reject (mu,phi,sigma) if proposed phi is outside (-1, 1)
  };
  const stochvol::ExpertSpec_FastSV expert_fac { // used for factor variances of fsv
    true,  // interweave
    stochvol::Parameterization::CENTERED,  // centered baseline always
    B011inv,  // B011inv,
    B022inv,  //B022inv,
    3,  // MHsteps,
    stochvol::ExpertSpec_FastSV::ProposalSigma2::INDEPENDENCE,  // independece proposal for sigma
    -1,  // MHcontrol unused for independence prior for sigma
    stochvol::ExpertSpec_FastSV::ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL  // immediately reject (mu,phi,sigma) if proposed phi is outside (-1, 1)
  };
  
  // initialize mixing indicators
  // mixind is also used by cholesky_sv
  arma::umat mixind(T, (M+factors));
  //initialize d_sqrt
  arma::mat d_sqrt = arma::exp(logvar.cols(0,M-1)/2);
  
  
  //------------------------------------STORAGE-------------------------------//
  
  const int nsave = std::floor(draws/thin);
  
  Rcpp::NumericVector logvar_draws_rcpp((tvp_keep == "all" ? logvar.n_rows : 1) * logvar.n_cols * nsave);
  logvar_draws_rcpp.attr("dim") = Rcpp::IntegerVector::create(tvp_keep == "all" ? logvar.n_rows : 1, logvar.n_cols, nsave);
  arma::cube logvar_draws(logvar_draws_rcpp.begin(), tvp_keep == "all" ? logvar.n_rows : 1, logvar.n_cols, nsave, false);
  Rcpp::NumericVector  sv_para_draws_rcpp(sv_para_in.nrow()*sv_para_in.ncol()*nsave);
  sv_para_draws_rcpp.attr("dim") = Rcpp::IntegerVector::create(sv_para_in.nrow(), sv_para_in.ncol(), nsave);
  arma::cube sv_para_draws(sv_para_draws_rcpp.begin(), sv_para_in.nrow(), sv_para_in.ncol(), nsave, false);
  Rcpp::NumericVector PHI_draws_rcpp(PHI.n_rows*PHI.n_cols*nsave);
  PHI_draws_rcpp.attr("dim") = Rcpp::IntegerVector::create(PHI.n_rows, PHI.n_cols,nsave);
  arma::cube PHI_draws(PHI_draws_rcpp.begin(),PHI.n_rows, PHI.n_cols,nsave, false);
  
  Rcpp::NumericMatrix U_draws_rcpp(sigma_type=="cholesky" ? U_upper_indices.n_elem : 0, sigma_type=="cholesky" ? nsave : 0);
  arma::mat U_draws(U_draws_rcpp.begin(), U_draws_rcpp.nrow(), U_draws_rcpp.ncol(), false);
  
  Rcpp::NumericVector facload_draws_rcpp((sigma_type=="factor" ? facload.nrow() : 0) *
    (sigma_type=="factor" ? facload.ncol() : 0) * (sigma_type=="factor" ? nsave : 0));
  facload_draws_rcpp.attr("dim") = Rcpp::IntegerVector::create(sigma_type=="factor" ? facload.nrow() : 0, sigma_type=="factor" ? facload.ncol() : 0, sigma_type=="factor" ? nsave : 0);
  arma::cube facload_draws(facload_draws_rcpp.begin(), sigma_type=="factor" ? facload.nrow() : 0, sigma_type=="factor" ? facload.ncol() : 0, sigma_type=="factor" ? nsave : 0, false);
  Rcpp::NumericVector fac_draws_rcpp((sigma_type=="factor" ? fac.nrow() : 0) * (sigma_type=="factor" ? (tvp_keep == "all" ? fac.ncol() : 1) : 0) * (sigma_type=="factor" ? nsave : 0));
  fac_draws_rcpp.attr("dim") = Rcpp::IntegerVector::create(sigma_type=="factor" ? fac.nrow() : 0, sigma_type=="factor" ? (tvp_keep == "all" ? fac.ncol() : 1) : 0, sigma_type=="factor" ? nsave : 0);
  arma::cube fac_draws(fac_draws_rcpp.begin(), sigma_type=="factor" ? fac.nrow() : 0, sigma_type=="factor" ? (tvp_keep == "all" ? fac.ncol() : 1) : 0, sigma_type=="factor" ? nsave : 0, false);
  
  int tvp_keep_start = 0;
  if(tvp_keep == "last") tvp_keep_start += (logvar.n_rows - 1);
  
  const bool isLambda = priorPHI == "DL" || priorPHI == "GT" || priorPHI == "HS" || priorPHI == "SSVS";
  Rcpp::NumericVector lambda_draws_rcpp(lambda_in.nrow()*lambda_in.ncol()*nsave*isLambda);
  lambda_draws_rcpp.attr("dim") = Rcpp::Dimension(lambda_in.nrow()*isLambda, lambda_in.ncol()*isLambda, nsave*isLambda);
  arma::cube lambda_draws(lambda_draws_rcpp.begin(), lambda_in.nrow()*isLambda, lambda_in.ncol()*isLambda, nsave*isLambda, false);
  
  const bool isPsi = priorPHI == "DL" || (priorPHI == "GT" && GT_priorkernel == "exponential") || priorPHI == "HS";
  Rcpp::NumericVector psi_draws_rcpp(psi_in.nrow()*psi_in.ncol()*nsave*isPsi);
  psi_draws_rcpp.attr("dim") = Rcpp::Dimension(psi_in.nrow()*isPsi, psi_in.ncol()*isPsi, nsave*isPsi);
  arma::cube psi_draws(psi_draws_rcpp.begin(), psi_in.nrow()*isPsi, psi_in.ncol()*isPsi, nsave*isPsi, false);
  
  const bool isXi = priorPHI == "GT" || (priorPHI == "DL" && DL_plus) || priorPHI == "HS";
  Rcpp::NumericMatrix xi_draws_rcpp(xi_in.length()*isXi, nsave*isXi);
  arma::mat xi_draws(xi_draws_rcpp.begin(), xi_draws_rcpp.nrow(), xi_draws_rcpp.ncol(), false);
  
  const bool isA = (priorPHI == "DL" && a_hyper) || (priorPHI == "GT" && a_hyper)  || (priorPHI == "SSVS" && a_hyper);
  Rcpp::NumericMatrix a_draws_rcpp(a_in.length()*isA, nsave*isA);
  arma::mat a_draws(a_draws_rcpp.begin(), a_draws_rcpp.nrow(), a_draws_rcpp.ncol(), false);
  
  const bool isVarpi = priorPHI == "HS";
  Rcpp::NumericMatrix varpi_draws_rcpp(varpi_in.length()*isVarpi, nsave*isVarpi);
  arma::mat varpi_draws(varpi_draws_rcpp.begin(), varpi_draws_rcpp.nrow(), varpi_draws_rcpp.ncol(), false);
  
  const bool isHMP = priorPHI == "HMP";
  Rcpp::NumericVector lambda_1_draws_rcpp(nsave*isHMP);
  Rcpp::NumericVector lambda_2_draws_rcpp(nsave*isHMP);
  
  const bool isLambda_U = (priorU == "DL" || priorU == "GT" || priorU == "HS" || priorU == "SSVS") && sigma_type == "cholesky";
  Rcpp::NumericMatrix lambda_U_draws_rcpp(lambda_U_in.length()*isLambda_U, nsave*isLambda_U);
  arma::mat lambda_U_draws(lambda_U_draws_rcpp.begin(), lambda_U_draws_rcpp.nrow(), lambda_U_draws_rcpp.ncol(), false);
  
  const bool isPsi_U = (priorU == "DL" || (priorU == "GT" && GT_priorkernel_U == "exponential") || priorU == "HS") && sigma_type == "cholesky";
  Rcpp::NumericMatrix psi_U_draws_rcpp(psi_U_in.length()*isPsi_U, nsave*isPsi_U);
  arma::mat psi_U_draws(psi_U_draws_rcpp.begin(), psi_U_draws_rcpp.nrow(), psi_U_draws_rcpp.ncol(), false);
  
  const bool isXi_U = priorU == "GT" || (priorU == "DL" && DL_plus_U) || priorU == "HS";
  Rcpp::NumericVector xi_U_draws_rcpp(nsave*isXi_U);
  
  const bool isA_U = (priorU == "DL" && a_hyper_U) || (priorU == "GT" && a_hyper_U)  || (priorU == "SSVS" && a_hyper_U);
  Rcpp::NumericVector a_U_draws_rcpp(nsave*isA_U);
  
  const bool isVarpi_U = priorU == "HS";
  Rcpp::NumericVector varpi_U_draws_rcpp(nsave*isVarpi_U);
  
  const bool isHMP_U = priorU == "HMP";
  Rcpp::NumericVector lambda_3_draws_rcpp(nsave*isHMP_U);
  
  const bool isFsvNg = sigma_type == "factor" && ngprior;
  Rcpp::NumericMatrix factor_lambda2_rcpp(armalambda2.n_elem*isFsvNg, nsave*isFsvNg);
  arma::mat factor_lambda2(factor_lambda2_rcpp.begin(), factor_lambda2_rcpp.nrow(), factor_lambda2_rcpp.ncol(), false);
  Rcpp::NumericVector factor_tau2_rcpp(tau2.nrow()*tau2.ncol()*nsave*isFsvNg);
  factor_tau2_rcpp.attr("dim") = Rcpp::Dimension(tau2.nrow()*isFsvNg, tau2.ncol()*isFsvNg, nsave*isFsvNg);
  arma::cube factor_tau2(factor_tau2_rcpp.begin(), tau2.nrow()*isFsvNg, tau2.ncol()*isFsvNg, nsave*isFsvNg, false);

  //indicator vector needed for DL and R2D2 on U
  arma::uvec ind_U(n_U);
  for(arma::uword i=0; i<n_U; ++i){
    ind_U(i)=i;
  }
  
  
  //-----------------------------------SAMPLER--------------------------------//
  const int tot = draws + burnin;
  // number of digits needed for printing progress
  const int printwidth = std::floor(std::log10(tot + .1))+ 1;
  
  for(int rep = 0; rep < tot; rep++){

    // Check for user interrupts
    if (rep % 16 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    if(progressbar){
      Rprintf("\r###  %*i / %*i ### (%3.0f%%) ###",
              printwidth, rep, printwidth, tot, 100.*rep/tot);
    }

    single_bvar(PHI, PHI_priorvariance, a, c, psi, lambda, xi, 
                varpi, lambda_1, lambda_2,
                U, u, d_sqrt, U_priorvariance, a_U, b_U, c_U, psi_U, lambda_U, xi_U, 
                lambda_3, varpi_U, armafacload, 
                armafac, armatau2, armalambda2, logvar, logvar0, sv_para, mixind,
                rep, burnin, Y, X, K, T, M, n_ol, n_cl, i_mat,
                i_ol, i_cl, groups, PHI_priormean, priorPHI, intercept, priorIntercept, a_hyper, 
                DL_plus, b, a_vec, a_weight, norm_consts, c_vec, 
                c_rel_a, GT_priorkernel, GT_vs, tau_0, tau_1, SSVS_s_a, SSVS_s_b, 
                s_r_1, s_r_2, V_i_prep, huge, sigma_type, U_upper_indices,
                ind_U, priorU, a_hyper_U, DL_plus_U,
                a_vec_U, a_weight_U, norm_consts_U, c_vec_U, c_rel_a_U, GT_priorkernel_U, 
                GT_vs_U, tau_0_U, tau_1_U, SSVS_s_a_U, SSVS_s_b_U, 
                s_r_3, aShrink, cShrink, dShrink, factor_facloadtol, factor_armarestr, 
                armafacloadtunrestrictedelements, nonzerospercol, nonzerosperrow,
                sv_priorh0, ngprior, columnwise, factor_homoskedastic, factor_offset,
                heteroscedastic, factor_interweaving, expert_sv, expert_fac,
                prior_specs, B011inv, priorHomoscedastic, sv_offset, abs_floor, stable,
                stable_lim);

    //-------Store draws after burnin
    if(rep >= burnin && ( (rep+1-burnin) % thin == 0 )){
      const int nthsave = (rep+1-burnin)/thin - 1;
      PHI_draws.slice(nthsave) = PHI;
      logvar_draws.slice(nthsave) = logvar.rows(tvp_keep_start ,logvar.n_rows-1);
      sv_para_draws.slice((nthsave)) = sv_para;
      
      if(isLambda) {
        lambda_draws.slice(nthsave) = lambda;
      }
      if(isPsi) {
        psi_draws.slice(nthsave) = psi;
      }
      
      if(isXi) {
        xi_draws.col(nthsave) = xi;
      }
      if(isA) {
        a_draws.col(nthsave) = a;
      }
      if(isVarpi) {
        varpi_draws.col(nthsave) = varpi;
      }
      if(isHMP) {
        lambda_1_draws_rcpp(nthsave) = lambda_1;
        lambda_2_draws_rcpp(nthsave) = lambda_2;
      }

      if(sigma_type == "cholesky"){
        U_draws.col(nthsave) = u;
        if(isLambda_U) {
          lambda_U_draws.col(nthsave) = lambda_U;
        }
        if(isPsi_U) {
          psi_U_draws.col(nthsave) = psi_U;
        }
        if(isXi_U) {
          xi_U_draws_rcpp(nthsave) = xi_U;
        }
        if(isA_U) {
          a_U_draws_rcpp(nthsave) = a_U;
        }
        if(isVarpi_U) {
          varpi_U_draws_rcpp(nthsave) = varpi_U;
        }
        if(isHMP_U) {
          lambda_3_draws_rcpp(nthsave) = lambda_3;
        }
      }else if(sigma_type == "factor" && factors > 0){
        facload_draws.slice(nthsave) = armafacload;
        fac_draws.slice(nthsave) = armafac.cols(tvp_keep_start ,armafac.n_cols-1);
        if(isFsvNg) {
          factor_lambda2.col(nthsave) = armalambda2;
          factor_tau2.slice(nthsave) = armatau2;
        }
      }
    }
  }
  
  if(progressbar){
    Rprintf("\r###  %*i / %*i ### (%3.0f%%) ###\n",
            printwidth, tot, printwidth, tot, 100.);
  }

  List out = List::create(
    Named("PHI") = PHI_draws_rcpp,
    Named("U") = U_draws_rcpp,
    Named("logvar") = logvar_draws_rcpp,
    Named("phi_para") = Rcpp::List::create(
      Rcpp::Named("lambda") = lambda_draws_rcpp,
      Rcpp::Named("psi") = psi_draws_rcpp,
      Rcpp::Named("a") = a_draws_rcpp,
      Rcpp::Named("xi") = xi_draws_rcpp,
      Rcpp::Named("varpi") = varpi_draws_rcpp,
      Rcpp::Named("lambda_1") = lambda_1_draws_rcpp,
      Rcpp::Named("lambda_2") = lambda_2_draws_rcpp
    ),
    Named("sigma_para") =  Rcpp::List::create(
      Rcpp::Named("U") = Rcpp::List::create(
        Rcpp::Named("lambda") = lambda_U_draws_rcpp,
        Rcpp::Named("psi") = psi_U_draws_rcpp,
        Rcpp::Named("a") = a_U_draws_rcpp,
        Rcpp::Named("xi") = xi_U_draws_rcpp,
        Rcpp::Named("varpi") = varpi_U_draws_rcpp,
        Rcpp::Named("lambda_3") = lambda_3_draws_rcpp
      ),
      Rcpp::Named("facload") = Rcpp::List::create(
        Rcpp::Named("tau2") = factor_tau2_rcpp,
        Rcpp::Named("lambda2") = factor_lambda2_rcpp
      ),
      Rcpp::Named("sv_para") = sv_para_draws_rcpp
    ),
    Named("facload") = facload_draws_rcpp,
    Named("fac") = fac_draws_rcpp,
    Named("lastdraw") = List::create(
      Named("phi") = List::create(
        Named("PHI") = PHI_in,
        Named("lambda") = lambda_in,
        Named("psi") = psi_in,
        Named("xi") = xi_in,
        Named("a") = a_in,
        Named("varpi") = varpi_in,
        Named("lambda_1") = lambda_1,
        Named("lambda_2") = lambda_2 
      ),
      Named("sigma") = List::create(
        Named("U") = U_in,
        Named("lambda_U") = lambda_U_in,
        Named("psi_U") = psi_U_in,
        Named("xi_U") = xi_U,
        Named("a_U") = a_U,
        Named("varpi_U") = varpi_U,
        Named("lambda_3") = lambda_3,
        Named("facload") = facload,
        Named("fac") = fac,
        Named("tau2") = tau2,
        Named("sv_logvar") = logvar_in,
        Named("sv_para")= sv_para_in,
        Named("sv_logvar0") = logvar0_in
      )
    )
  );

  return out;
}
