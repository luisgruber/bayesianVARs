#include <RcppArmadillo.h>
#include <factorstochvol.h>
#include <stochvol.h>
#include "bvar_single_update.h"

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
              const int& intercept,
              const arma::vec& priorIntercept,
              arma::mat& PHI0,
              const List priorPHI_in,
              const List priorSigma_in,
              const List Rstartvals_in,
              const arma::imat& i_mat,
              const bool& progressbar,
              const bool& huge,
              const double& abs_floor){
  
  // Rcpp objects reuse memory
  // arma make copy by default unless deactivated
  
  //-------------------------Preliminaries--------------------------------------//
  
  const double n = K*M; // number of VAR coefficients without intercept
  const arma::uvec i_ocl= arma::find(i_mat != 0); // indicator for all coefficients except intercept
  const arma::uvec i_i= arma::find(i_mat == 0); // indicator for intercepts
  
  const std::string sigma_type = priorSigma_in["type"];
  
  //////////////////
  // indicators for coefficients without intercept
  const arma::ivec i_vec_small=i_mat(i_ocl);
  // indicators for cross-/own-lags (relevant for HM prior)
  const arma::uvec i_ol = arma::find(i_vec_small > 0);
  const arma::uvec i_cl = arma::find(i_vec_small < 0);
  const int n_ol = i_ol.size(); // nr. of ownlags
  const int n_cl = i_cl.size(); // nr. of crosslags
  /////////////////
  
  //--------------------Initialization -----------------------//
  
  List startvals_in(clone(Rstartvals_in)); // clone makes copy
  //---- PHI
  
  Rcpp::List phi_startvals = startvals_in["phi"];
  NumericMatrix PHI_in = phi_startvals["PHI"];
  arma::mat PHI(PHI_in.begin(), PHI_in.nrow(), PHI_in.ncol(), false);
  arma::mat PHI_diff = PHI - PHI0; // will hold PHI - PHI0
  std::string priorPHI = priorPHI_in["prior"];
  
  // V_i holds prior variances (without intercepts)
  arma::vec V_i(n);
  
  arma::vec V_i_long(n+M*intercept); //V_i plus intercept prior variances
  V_i_long(i_i) = priorIntercept;
  
  if(priorPHI == "normal"){
    //arma::vec V_i_in = priorPHI_in["V_i"];
    NumericVector V_i_in = priorPHI_in["V_i"];
    arma::vec armaV_i_in(V_i_in.begin(), V_i_in.size(), false);
    //in case of 'normal' V_i is fixed at user specified values
    V_i = armaV_i_in;
  }
  
  //  sub-groups for semi-global-local GL priors and SSVS prior
  const int n_groups = priorPHI_in["n_groups"]; // number of distinct groups
  IntegerVector groups_in = priorPHI_in["groups"];
  arma::ivec groups(groups_in.begin(), groups_in.size(),false);
  
  // parameters shared by GT (i.e. NG & R2D2) and DL priors
  NumericVector lambda_in = phi_startvals["lambda"];
  arma::vec lambda(lambda_in.begin(), lambda_in.length(), false);//(n, fill::value(1/static_cast<double>(n)));
  NumericVector psi_in = phi_startvals["psi"];
  arma::vec psi(psi_in.begin(), psi_in.length(), false);//(n, fill::ones); // must be 1 for NG prior (GT with normal kernel) ? isn't the code clever enough
  NumericVector xi_in = phi_startvals["xi"];
  arma::vec xi(xi_in.begin(), xi_in.length(), false);//(n_groups, fill::value(0.5));
  NumericVector a_in = phi_startvals["a"];
  arma::vec a(a_in.begin(), a_in.length(), false);//priorPHI_in["a"]; // arma makes copy...is that necessary?
  const arma::vec b = priorPHI_in["b"]; // arma makes copy...is that necessary?
  arma::vec c = priorPHI_in["c"]; // arma makes copy...is that necessary?
  
  NumericVector a_vec_in = priorPHI_in["a_vec"]; //const
  arma::vec a_vec(a_vec_in.begin(), a_vec_in.size(), false);
  NumericVector a_weight_in = priorPHI_in["a_weight"]; //const
  arma::vec a_weight(a_weight_in.begin(), a_weight_in.size(), false);
  NumericVector norm_consts_in = priorPHI_in["norm_consts"];
  arma::vec norm_consts(norm_consts_in.begin(), norm_consts_in.size(), false); //const
  NumericVector c_vec_in = priorPHI_in["c_vec"];
  arma::vec c_vec(c_vec_in.begin(), c_vec_in.size(), false); // const
  const bool c_rel_a = priorPHI_in["c_rel_a"];
  
  //---- DL prior specific parameters
  const bool DL_hyper = priorPHI_in["DL_hyper"];
  const bool DL_plus = priorPHI_in["DL_plus"];
  //arma::vec zeta(n_groups); zeta.fill(10);
  if(priorPHI == "DL" ){
    V_i = psi % lambda % lambda;// * zeta(0) * zeta(0);
  }
  
  //----GT on PHI (GT refers to gamma type priors a la normal gamma or r2d2)
  const bool GT_hyper = priorPHI_in["GT_hyper"];
  const double GT_vs = priorPHI_in["GT_vs"];
  const std::string GT_priorkernel = priorPHI_in["GT_priorkernel"];
  if(priorPHI == "GT"){
    V_i = lambda*GT_vs;
    if(GT_priorkernel == "exponential"){
      V_i %= psi;// /2;
    }
    // else if(GT_priorkernel == "normal"){
    // V_i = lambda*GT_vs;
    // }
  }
  
  if(c_rel_a && (GT_hyper || DL_hyper)) {
    for(int g=0; g<n_groups; ++g){
      arma::vec tmp = c_vec.elem(arma::find(a_vec == a(g)));
      c(g) = tmp(0);
    }
  }
  
  //---- Horseshoe on PHI
  // initialization of global, local and auxiliary scaling parameters
  // NumericVector theta_hs_in = phi_startvals["theta_hs"];
  // arma::vec theta_hs(theta_hs_in.begin(), theta_hs_in.length(), false); //theta_hs.fill(1/static_cast<double>(n));
  // NumericVector zeta_hs_in = phi_startvals["zeta_hs"];
  // arma::vec zeta_hs(zeta_hs_in.begin(), zeta_hs_in.length(), false);//(n_groups); zeta_hs.fill(10);
  // NumericVector nu_in = phi_startvals["nu"];
  // arma::vec nu(nu_in.begin(), nu_in.length(), false);//(n); nu.fill(1);
  NumericVector varpi_in = phi_startvals["varpi"];
  arma::vec varpi(varpi_in.begin(), varpi_in.length(), false);//(n_groups); varpi.fill(1);
  if(priorPHI == "HS"){
    // V_i = lambda*xi(0);
    arma::ivec::iterator g;
    int j=0;
    
    for(g=groups.begin(); g!=groups.end(); ++g){
      arma::uvec ind = arma::find(i_vec_small == *g);
      V_i(ind) = lambda(ind) * xi(j);
      j++;
    }
  }
  
  //---- SSVS on PHI
  arma::vec tau_0 = priorPHI_in["SSVS_tau0"];
  arma::vec tau_1 = priorPHI_in["SSVS_tau1"];
  double SSVS_s_a = priorPHI_in["SSVS_s_a"];
  double SSVS_s_b = priorPHI_in["SSVS_s_b"];
  bool SSVS_hyper = priorPHI_in["SSVS_hyper"];
  // NumericVector p_i_in = phi_startvals["p_i"];
  // arma::vec p_i(p_i_in.begin(), p_i_in.length(), false);//= priorPHI_in["SSVS_p"];
  
  if(priorPHI == "SSVS"){
    // if(SSVS_hyper){
    //   V_i = tau_1 % tau_1;
    // }else{
    //   V_i = tau_0 % tau_0; // !!! tau_0 % tau_0
    // }
    V_i = arma::square(tau_1%lambda + (lambda-1)%tau_0);
  }
  // Rcpp::NumericVector gammas_in = phi_startvals["gammas"];
  // arma::vec gammas(gammas_in.begin(), gammas_in.length(), false);//(n, fill::zeros);
  
  //---- HMP on PHI
  double lambda_1= phi_startvals["lambda_1"];//0.04;
  double lambda_2= phi_startvals["lambda_2"];//0.0016;
  arma::vec V_i_prep = priorPHI_in["V_i_prep"];
  arma::vec s_r_1 = priorPHI_in["lambda_1"];
  arma::vec s_r_2 = priorPHI_in["lambda_2"];
  if(priorPHI == "HMP"){
    //V_i_long(i_i) = priorIntercept % V_i_prep(i_i);
    // arma::vec V_i_prep_small = V_i_prep(i_ocl);
    V_i(i_ol) = lambda_1*V_i_prep(i_ol);//_small(i_ol);
    V_i(i_cl) = lambda_2*V_i_prep(i_cl);//_small(i_cl);
  }
  
  //Fill V_i_long with remaining prior variances
  V_i_long(i_ocl) = V_i;
  arma::mat V_prior(V_i_long.begin(), K + intercept, M, false); // pointer, reuses memory, more efficient than reshape!!!
  //arma::mat V_prior = arma::reshape(V_i_long, K + intercept, M);
  
  //-------------- U
  List sigma_startvals = startvals_in["sigma"];
  NumericMatrix U_in = sigma_startvals["U"];
  arma::mat U(U_in.begin(), U_in.nrow(), U_in.ncol(), false);
  
  std::string priorU = priorSigma_in["cholesky_U_prior"];
  
  arma::uvec U_upper_indices(1);
  if(sigma_type == "cholesky"){
    U_upper_indices = trimatu_ind( size(U),  1);
  }
  arma::vec u = U(U_upper_indices);
  const double n_U = u.size();
  arma::vec V_i_U(n_U);
  
  if(priorU == "normal"){
    arma::vec V_i_U_in = priorSigma_in["cholesky_V_i"];
    V_i_U = V_i_U_in;
  }
  
  /// Initialize GL parameters
  NumericVector lambda_U_in = sigma_startvals["lambda_U"];
  arma::vec lambda_U(lambda_U_in.begin(), lambda_U_in.length(), false);//(n_U, fill::value(1/static_cast<double>(n_U)));
  NumericVector psi_U_in = sigma_startvals["psi_U"];
  arma::vec psi_U(psi_U_in.begin(), psi_U_in.length(), false);//(n_U, fill::ones);
  double xi_U = sigma_startvals["xi_U"];//0.5;
  double a_U = sigma_startvals["a_U"];//priorSigma_in["cholesky_a"];
  const double b_U = priorSigma_in["cholesky_b"];
  double c_U = priorSigma_in["cholesky_c"];
  const arma::vec c_vec_U = priorSigma_in["cholesky_c_vec"];
  const bool c_rel_a_U = priorSigma_in["cholesky_c_rel_a"];
  ///
  
  //---- DL prior on U
  bool DL_hyper_U = priorSigma_in["cholesky_DL_hyper"];
  
  const bool DL_plus_U = priorSigma_in["cholesky_DL_plus"];
  const arma::vec a_vec_U = priorSigma_in["cholesky_a_vec"];
  const arma::vec a_weight_U = priorSigma_in["cholesky_a_weight"];
  const arma::vec norm_consts_U = priorSigma_in["cholesky_norm_consts"];
  if(priorU == "DL"){
    V_i_U= psi_U % lambda_U % lambda_U;
  }
  
  //----GT on U (GT refers to gamma type priors a la normal gamma or r2d2)
  const double GT_vs_U = priorSigma_in["cholesky_GT_vs"];
  const bool GT_hyper_U = priorSigma_in["cholesky_GT_hyper"];
  const std::string GT_priorkernel_U = priorSigma_in["cholesky_GT_priorkernel"];
  if(priorU == "GT"){
    V_i_U = lambda_U*GT_vs_U;
    if(GT_priorkernel_U == "exponential"){
      V_i_U %= psi_U;
    }
    // else if(GT_priorkernel_U == "normal"){
    //   V_i_U = lambda_U;
    // }
  }
  if(c_rel_a_U && (GT_hyper_U || DL_hyper_U)) {
    arma::vec tmp = c_vec_U.elem(arma::find(a_vec_U == a_U));
    c_U = tmp(0);
  }
  
  //---- Horseshoe on U
  // initialization of global, local and auxiliary scaling parameters
  // NumericVector theta_hs_U_in = sigma_startvals["theta_hs_U"];
  // arma::vec theta_hs_U(theta_hs_U_in.begin(), theta_hs_U_in.length(), false);//(n_U); theta_hs_U.fill(1/static_cast<double>(n_U));
  // double zeta_hs_U = sigma_startvals["zeta_hs_U"];//= 10;
  // NumericVector nu_U_in = sigma_startvals["nu_U"];
  // arma::vec nu_U(nu_U_in.begin(), nu_U_in.length(), false);//(n_U, fill::ones);
  double varpi_U = sigma_startvals["varpi_U"];//= 1;
  if(priorU == "HS"){
    V_i_U = lambda_U*xi_U;
  }
  
  //---- SSVS on U
  bool SSVS_hyper_U = priorSigma_in["cholesky_SSVS_hyper"];
  // NumericVector p_i_U_in = sigma_startvals["p_i_U"];
  // arma::vec p_i_U(p_i_U_in.begin(), p_i_U_in.length(), false);//= priorSigma_in["cholesky_SSVS_p"];
  arma::vec tau_0_U = priorSigma_in["cholesky_SSVS_tau0"];
  arma::vec tau_1_U = priorSigma_in["cholesky_SSVS_tau1"];
  double SSVS_s_a_U = priorSigma_in["cholesky_SSVS_s_a"];
  double SSVS_s_b_U = priorSigma_in["cholesky_SSVS_s_b"];
  if(priorU == "SSVS"){
    // V_i_U = tau_0_U % tau_0_U;
    V_i_U = arma::square(tau_1_U%lambda_U + (lambda_U-1)%tau_0_U);
  }
  // NumericVector gammas_U_in = sigma_startvals["gammas_U"];
  // arma::vec gammas_U(gammas_U_in.begin(), gammas_U_in.length(), false);//(n_U, fill::zeros);
  
  //---- HMP on U
  double lambda_3 = sigma_startvals["lambda_3"];//0.001;
  NumericVector s_r_3 = priorSigma_in["cholesky_lambda_3"];
  if(priorU == "HMP"){
    arma::vec V_i_U_tmp(n_U); V_i_U_tmp.fill(1.0);
    V_i_U= lambda_3*V_i_U_tmp;
  }
  //-----------------------factorstochvol----------------------//
  const int factors = priorSigma_in["factor_factors"];
  if(sigma_type == "cholesky" && factors != 0){
    Rcpp::stop("cholesky with factors,...,");
  }
  
  // The objects preceded by //!!!/// are needed as inputs for factorstochvol::update_fsv()
  // initialization of factor loadings, factors, and tau2(shrinkage hyperparameter for factorloadings under ng-prior)
  // List factor_startval = startvals_in["factor_startval"];
  NumericMatrix facload = sigma_startvals["facload"];
  //!!!///
  arma::mat armafacload(facload.begin(), facload.nrow(), facload.ncol(), false);
  NumericMatrix fac = sigma_startvals["fac"];
  //!!!///
  arma::mat armafac(fac.begin(), fac.nrow(), fac.ncol(), false);
  
  NumericMatrix tau2 = sigma_startvals["tau2"];
  //!!!///
  arma::mat armatau2(tau2.begin(), tau2.nrow(), tau2.ncol(), false);
  
  //----------  prior settings
  // shrinkage hyperparamters for factor loadings
  const List factor_shrinkagepriors = priorSigma_in["factor_shrinkagepriors"];
  //!!!///
  const NumericVector aShrink = factor_shrinkagepriors["a"]; // a is hyperparameter of local prior
  //!!!///
  const NumericVector cShrink = factor_shrinkagepriors["c"]; // c and d are hyperparameters of global prior
  //!!!///
  const NumericVector dShrink = factor_shrinkagepriors["d"];
  
  //!!!///
  const bool ngprior = priorSigma_in["factor_ngprior"];
  //!!!///
  const bool columnwise = priorSigma_in["factor_columnwise"];
  //!!!//
  arma::vec armalambda2((ngprior && columnwise) ? factors : ngprior ? M : 0 ); // conditional statement: if ngprior && columnwise, then factors, if ngprior (not columnwise), then M, else 0
  
  // underbound of absolute value of factorloadings for better preventing numerical issues
  //!!!///
  const double factor_facloadtol = priorSigma_in["factor_facloadtol"];
  
  //!!!///
  /* Could there be any situation where I need an offset? */
  const double factor_offset = 0; // unify with cholesky_sv ???
  
  // restriction on factor loadings
  IntegerMatrix factor_restriction = priorSigma_in["factor_restrinv"];
  //!!!///
  arma::imat factor_armarestr(factor_restriction.begin(), factor_restriction.nrow(), factor_restriction.ncol(), false);
  //!!!///
  arma::uvec armafacloadtunrestrictedelements = arma::find(factor_armarestr.t() != 0);
  //!!!///
  arma::irowvec nonzerospercol = arma::sum(factor_armarestr, 0);
  //!!!///
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
  //!!!///
  NumericVector heteroscedastic = priorSigma_in["sv_heteroscedastic"];
  
  // in case of constant idiosyncratic variance, homoscedastic collects shape and scale parameters of inverse-gamma prior of idi variances
  //!!!///
  NumericMatrix factor_homoskedastic = priorSigma_in["factor_priorhomoskedastic"]; // unify with cholesky_sv ???
  //!!!///
  const int factor_interweaving = priorSigma_in["factor_interweaving"]; // factor
  
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
  /// NumericVector sv_offset = sv_spec["sv_offset"];
  NumericVector sv_offset = priorSigma_in["cholesky_sv_offset"];
  
  //!!!///
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
  //!!!//
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
  //!!!//
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
  //!!!// // mixind is also used by cholesky_sv
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
  Rcpp::NumericVector V_prior_draws_rcpp(V_prior.n_rows * V_prior.n_cols * nsave);
  V_prior_draws_rcpp.attr("dim") = Rcpp::IntegerVector::create(V_prior.n_rows, V_prior.n_cols, nsave);
  arma::cube V_prior_draws(V_prior_draws_rcpp.begin(), V_prior.n_rows, V_prior.n_cols, nsave, false);
  
  //arma::cube U_draws(sigma_type=="cholesky" ? nsave : 0, sigma_type=="cholesky" ? U.n_rows : 0, sigma_type=="cholesky" ? U.n_cols : 0);
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
  
  // maybe add storage for tau2 and lambda2, shrinkage hyperparameters of ngpior for factorloadings
  
  int phi_hyperparameter_size(0);
  if(priorPHI == "DL" ){
    phi_hyperparameter_size += 2*n; // a + xi + n(lambda + psi)
    if(DL_plus){
      phi_hyperparameter_size += n_groups;
    }
    if(DL_hyper){
      phi_hyperparameter_size += n_groups;
    }
  }else if(priorPHI == "GT"){
    phi_hyperparameter_size += n_groups + n; // a+xi + n(lambda + psi)
    if(GT_priorkernel == "exponential"){
      phi_hyperparameter_size += n;
    }
    if(GT_hyper) {
      phi_hyperparameter_size += n_groups;
    }
  }else if(priorPHI == "HS"){
    phi_hyperparameter_size += 2*n_groups + 2*n; // zeta + varpi + n(theta + nu)
  }else if(priorPHI == "SSVS"){
    phi_hyperparameter_size +=  n; // n_groups a + n gammas
    if(SSVS_hyper){
      phi_hyperparameter_size += n_groups;
    }
  }else if(priorPHI == "HMP"){
    phi_hyperparameter_size += 2; // lambda_1 + lambda_2
  }
  Rcpp::NumericMatrix phi_hyperparameter_draws_rcpp(phi_hyperparameter_size,
                                                    priorPHI != "normal" ? nsave : 0);
  arma::mat phi_hyperparameter_draws(phi_hyperparameter_draws_rcpp.begin(),
                                     phi_hyperparameter_draws_rcpp.nrow(),
                                     phi_hyperparameter_draws_rcpp.ncol(), false);
  // store info, see store after burnin
  const int offset_a = (DL_hyper || GT_hyper) ? 1 : 0;
  const int offset_xi = (DL_plus || priorPHI == "GT") ? 1 : 0;

  int u_hyperparameter_size(0);
  if(sigma_type == "cholesky"){
    if(priorU == "DL" ){
      u_hyperparameter_size += 2*n_U;
      if(DL_plus_U) {
        u_hyperparameter_size += 1;
      }
      if(DL_hyper_U) {
        u_hyperparameter_size += 1;
      }
    }else if(priorU == "GT"){
      u_hyperparameter_size += 1 + n_U; // a+xi + n(lambda + psi)
      if(GT_priorkernel_U == "exponential") {
        u_hyperparameter_size += n_U;
      }
      if(GT_hyper_U) {
        u_hyperparameter_size += 1;
      }
    }else if(priorU == "SSVS"){
      u_hyperparameter_size += n_U;
      if(SSVS_hyper_U) {
        u_hyperparameter_size += 1;
      }
    }else if(priorU == "HMP"){
      u_hyperparameter_size += 1;
    }else if(priorU == "HS"){
      u_hyperparameter_size += 2+2*n_U;
    }
  }
  Rcpp::NumericMatrix u_hyperparameter_draws_rcpp(sigma_type=="cholesky" ? u_hyperparameter_size : 0, (sigma_type=="cholesky" && priorU != "normal") ? nsave : 0);
  arma::mat u_hyperparameter_draws(u_hyperparameter_draws_rcpp.begin(),
                                   u_hyperparameter_draws_rcpp.nrow(),
                                   u_hyperparameter_draws_rcpp.ncol(), false);
  // store info, see store after burnin
  const int offset_a_U = (DL_hyper_U || GT_hyper_U) ? 1 : 0;
  const int offset_xi_U = (DL_plus_U || priorU == "GT") ? 1 : 0;

  //indicator vector needed for DL and R2D2 on U
  arma::uvec ind_U(n_U);
  for(int i=0; i<n_U; ++i){
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

    single_bvar(PHI, PHI_diff, V_prior, V_i, V_i_long, a, c, psi, lambda, xi, 
                varpi, lambda_1, lambda_2,
                U, u, d_sqrt, V_i_U, a_U, b_U, c_U, psi_U, lambda_U, xi_U, 
                lambda_3, varpi_U, armafacload, 
                armafac, armatau2, armalambda2, logvar, logvar0, sv_para, mixind,
                rep, burnin, Y, X, K, T, M, n_ol, n_cl, i_vec_small, i_ocl,
                i_ol, i_cl, groups, PHI0, priorPHI, DL_hyper, 
                DL_plus, GT_hyper, b, a_vec, a_weight, norm_consts, c_vec, 
                c_rel_a, GT_priorkernel, GT_vs, tau_0, tau_1, SSVS_s_a, SSVS_s_b, 
                SSVS_hyper, s_r_1, s_r_2, V_i_prep, huge, sigma_type, U_upper_indices,
                ind_U, priorU, DL_hyper_U, DL_plus_U, GT_hyper_U,
                a_vec_U, a_weight_U, norm_consts_U, c_vec_U, c_rel_a_U, GT_priorkernel_U, 
                GT_vs_U, tau_0_U, tau_1_U, SSVS_s_a_U, SSVS_s_b_U, SSVS_hyper_U, 
                s_r_3, aShrink, cShrink, dShrink, factor_facloadtol, factor_armarestr, 
                armafacloadtunrestrictedelements, nonzerospercol, nonzerosperrow,
                sv_priorh0, ngprior, columnwise, factor_homoskedastic, factor_offset,
                heteroscedastic, factor_interweaving, expert_sv, expert_fac,
                prior_specs, B011inv, priorHomoscedastic, sv_offset, abs_floor);

    //-------Store draws after burnin
    if(rep >= burnin && ( (rep+1-burnin) % thin == 0 )){
      const int nthsave = (rep+1-burnin)/thin - 1;
      PHI_draws.slice(nthsave) = PHI;
      logvar_draws.slice(nthsave) = logvar.rows(tvp_keep_start ,logvar.n_rows-1);
      sv_para_draws.slice((nthsave)) = sv_para;
      V_prior_draws.slice(nthsave) = V_prior;

      if(priorPHI == "DL" || priorPHI== "GT"){

        if(DL_hyper || GT_hyper) {
          phi_hyperparameter_draws(span(0,(n_groups-1)), nthsave) = a;
        }
        if(DL_plus || priorPHI == "GT") {
          phi_hyperparameter_draws(span(offset_a*n_groups,((1+offset_a)*(n_groups)-1)), nthsave) = xi;
          // offset_a += 1;
        }
        phi_hyperparameter_draws(span((offset_xi+offset_a)*n_groups,((offset_xi+offset_a)*n_groups+n-1)), nthsave) = lambda;
        if(priorPHI == "DL" || GT_priorkernel == "exponential"){
          phi_hyperparameter_draws(span((phi_hyperparameter_size-n),(phi_hyperparameter_size-1)), nthsave) = psi;
        }

      }else if(priorPHI == "HS"){

        phi_hyperparameter_draws(span(0,(n_groups-1)), nthsave) = xi;
        phi_hyperparameter_draws(span(n_groups,(2*n_groups-1)), nthsave) = varpi;
        phi_hyperparameter_draws(span(2*n_groups,(2*n_groups+n-1)), nthsave) = lambda;
        phi_hyperparameter_draws(span((2*n_groups+n),(phi_hyperparameter_size-1)), nthsave) = psi;

      }else if(priorPHI == "SSVS"){

        if(SSVS_hyper) {
          phi_hyperparameter_draws(span(0,(n_groups-1)), nthsave) = a;

        }
        phi_hyperparameter_draws(span(phi_hyperparameter_size-n, phi_hyperparameter_size-1), nthsave) = lambda;

      }else if(priorPHI == "HMP"){
        phi_hyperparameter_draws(0, nthsave) = lambda_1;
        phi_hyperparameter_draws(1, nthsave) = lambda_2;
      }

      if(sigma_type == "cholesky"){

        U_draws.col(nthsave) = u;

        if(priorU == "DL" || priorU == "GT"){

          if(DL_hyper_U || GT_hyper_U) {
            u_hyperparameter_draws(0, nthsave) = a_U;
          }
          if(DL_plus_U || priorU == "GT"){
            u_hyperparameter_draws(offset_a_U, nthsave) = xi_U;
          }
          u_hyperparameter_draws(span((offset_a_U + offset_xi_U),(offset_a_U +  offset_xi_U + n_U -1)), nthsave) = lambda_U;
          if(priorU == "DL" || GT_priorkernel_U == "exponential"){
            u_hyperparameter_draws(span(u_hyperparameter_size-n_U,u_hyperparameter_size-1), nthsave) = psi_U;
          }

        }else if(priorU == "SSVS"){

          if(SSVS_hyper_U) {
            u_hyperparameter_draws(0, nthsave) = a_U;
          }
          u_hyperparameter_draws(span(u_hyperparameter_size-n_U, u_hyperparameter_size-1), nthsave) = lambda_U;

        }else if(priorU == "HMP"){

          u_hyperparameter_draws(0, nthsave) = lambda_3;

        }else if(priorU == "HS"){

          u_hyperparameter_draws(0, nthsave) = xi_U ;
          u_hyperparameter_draws(1, nthsave) = varpi_U;
          u_hyperparameter_draws(span(2, n_U+1), nthsave) = lambda_U;
          u_hyperparameter_draws(span((n_U+2),(u_hyperparameter_size-1)), nthsave)= psi_U;

        }
      }else if(sigma_type == "factor" && factors > 0){
        facload_draws.slice(nthsave) = armafacload;
        fac_draws.slice(nthsave) = armafac.cols(tvp_keep_start ,armafac.n_cols-1);
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
    Named("sv_para") = sv_para_draws_rcpp,
    Named("phi_hyperparameter") = phi_hyperparameter_draws_rcpp,
    Named("u_hyperparameter") = u_hyperparameter_draws_rcpp,
    Named("V_prior") = V_prior_draws_rcpp,
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
