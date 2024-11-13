#include <RcppArmadillo.h>
#include <factorstochvol.h>
#include <progress.hpp>
#include <Rcpp/Benchmark/Timer.h>
#include "sample_coefficients.h"
// <stochvol.h> is already included in <factorstochvol.h>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
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
              const arma::vec priorIntercept,
              arma::mat& PHI0, // prior mean
              const List priorPHI_in,
              const List priorSigma_in, // priorSigma_in
              const List Rstartvals_in,
              const arma::imat& i_mat,
              const arma::ivec& i_vec,
              const bool& progressbar,
              const double& PHI_tol,
              const double& L_tol,
              const bool& huge){

  // Rcpp objects reuse memory
  // arma make copy by default unless deactivated

//-------------------------Preliminaries--------------------------------------//

  const double n = K*M; // number of VAR coefficients without intercept
  const arma::uvec i_ocl= arma::find(i_vec != 0); // indicator for all coefficients except intercept
  const arma::uvec i_i= arma::find(i_vec == 0); // indicator for intercepts

  const std::string sigma_type = priorSigma_in["type"];

  //////////////////
  // indicators for coefficients without intercept
  const arma::ivec i_vec_small=i_vec(i_ocl);
  // indicators for cross-/own-lags (relevant for HM prior)
  const arma::uvec i_ol = arma::find(i_vec_small > 0);
  const arma::uvec i_cl = arma::find(i_vec_small < 0);
  const int n_ol = i_ol.size(); // nr. of ownlags
  const int n_cl = i_cl.size(); // nr. of crosslags
  /////////////////

//--------------------Initialization -----------------------//

  List startvals_in(clone(Rstartvals_in)); // clone makes copy
//---- PHI

  NumericMatrix PHI_in = startvals_in["PHI"];
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
  arma::vec lambda(n, fill::value(1/static_cast<double>(n)));
  arma::vec psi(n, fill::ones); // must be 1 for NG prior (GT with normal kernel)
  arma::vec xi(n_groups, fill::value(0.5));
  arma::vec a = priorPHI_in["a"]; // arma makes copy...is that necessary?
  arma::vec b = priorPHI_in["b"]; // arma makes copy...is that necessary?
  arma::vec c = priorPHI_in["c"]; // arma makes copy...is that necessary?
  const double GL_tol = priorPHI_in["GL_tol"];

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
    if(GT_priorkernel == "exponential"){
      V_i = psi%lambda/2;
    }else if(GT_priorkernel == "normal"){
      V_i = lambda;
    }
  }

  //---- Horseshoe on PHI
  // initialization of global, local and auxiliary scaling parameters
  arma::vec theta_hs(n); theta_hs.fill(1/static_cast<double>(n));
  arma::vec zeta_hs(n_groups); zeta_hs.fill(10);
  arma::vec nu(n); nu.fill(1);
  arma::vec varpi(n_groups); varpi.fill(1);
  if(priorPHI == "HS"){
    V_i = theta_hs*zeta_hs(0);
  }

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
  V_i_long(i_ocl) = V_i;
  arma::mat V_prior(V_i_long.begin(), K + intercept, M, false); // pointer, reuses memory, more efficient than reshape!!!
  //arma::mat V_prior = arma::reshape(V_i_long, K + intercept, M);

//-------------- U

  NumericMatrix U_in = startvals_in["U"];
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
  arma::vec lambda_U(n_U, fill::value(1/static_cast<double>(n_U)));
  arma::vec psi_U(n_U, fill::ones);
  double xi_U = 0.5;
  double a_U = priorSigma_in["cholesky_a"];
  double b_U = priorSigma_in["cholesky_b"];
  double c_U = priorSigma_in["cholesky_c"];
  const double GL_tol_U = priorSigma_in["cholesky_GL_tol"];
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
    V_i_U= psi_U % lambda_U;
  }

  //----GT on U (GT refers to gamma type priors a la normal gamma or r2d2)
  const double GT_vs_U = priorSigma_in["cholesky_GT_vs"];
  const double GT_hyper_U = priorSigma_in["cholesky_GT_hyper"];
  const std::string GT_priorkernel_U = priorSigma_in["cholesky_GT_priorkernel"];
  if(priorU == "GT"){
    if(GT_priorkernel_U == "exponential"){
      V_i_U = psi_U%lambda_U/2;
    }else if(GT_priorkernel_U == "normal"){
      V_i_U = lambda_U;
    }
  }

  //---- Horseshoe on U
  // initialization of global, local and auxiliary scaling parameters
  arma::vec theta_hs_U(n_U); theta_hs_U.fill(1/static_cast<double>(n_U));
  double zeta_hs_U = 10;
  arma::vec nu_U(n_U, fill::ones);
  double varpi_U = 1;
  if(priorU == "HS"){
    V_i_U = theta_hs_U*zeta_hs_U;
  }

  //---- SSVS on U
  bool SSVS_hyper_U = priorSigma_in["cholesky_SSVS_hyper"];
  arma::vec p_i_U = priorSigma_in["cholesky_SSVS_p"];
  arma::vec tau_0_U = priorSigma_in["cholesky_SSVS_tau0"];
  arma::vec tau_1_U = priorSigma_in["cholesky_SSVS_tau1"];
  double SSVS_s_a_U = priorSigma_in["cholesky_SSVS_s_a"];
  double SSVS_s_b_U = priorSigma_in["cholesky_SSVS_s_b"];
  if(priorU == "SSVS"){
    V_i_U = tau_0_U % tau_0_U;
  }
  arma::vec gammas_U(n_U, fill::zeros);

  //---- HMP on U
  double lambda_3 = 0.001;
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
  List factor_startval = startvals_in["factor_startval"];
  NumericMatrix facload = factor_startval["facload"];
  //!!!///
  arma::mat armafacload(facload.begin(), facload.nrow(), facload.ncol(), false);
  NumericMatrix fac = factor_startval["fac"];
  //!!!///
  arma::mat armafac(fac.begin(), fac.nrow(), fac.ncol(), false);

  NumericMatrix tau2 = factor_startval["tau2"];
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

  NumericMatrix logvar_in = startvals_in["sv_logvar"];
  arma::mat logvar(logvar_in.begin(), logvar_in.nrow(), logvar_in.ncol(), false);
  NumericMatrix sv_para_in = startvals_in["sv_para"];
  arma::mat sv_para(sv_para_in.begin(), sv_para_in.nrow(), sv_para_in.ncol(), false);
  NumericVector logvar0_in = startvals_in["sv_logvar0"];
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

  arma::cube logvar_draws(tvp_keep == "all" ? logvar.n_rows : 1, logvar.n_cols, nsave);
  arma::cube sv_para_draws(sv_para_in.nrow(), sv_para_in.ncol(), nsave);
  arma::cube PHI_draws(PHI.n_rows, PHI.n_cols,nsave);
  arma::cube V_prior_draws(V_prior.n_rows, V_prior.n_cols, nsave);

  //arma::cube U_draws(sigma_type=="cholesky" ? nsave : 0, sigma_type=="cholesky" ? U.n_rows : 0, sigma_type=="cholesky" ? U.n_cols : 0);
  arma::mat U_draws(sigma_type=="cholesky" ? U_upper_indices.n_elem : 0, sigma_type=="cholesky" ? nsave : 0);

  arma::cube facload_draws(sigma_type=="factor" ? facload.nrow() : 0, sigma_type=="factor" ? facload.ncol() : 0, sigma_type=="factor" ? nsave : 0);
  arma::cube fac_draws(sigma_type=="factor" ? fac.nrow() : 0, sigma_type=="factor" ? (tvp_keep == "all" ? fac.ncol() : 1) : 0, sigma_type=="factor" ? nsave : 0);

  int tvp_keep_start = 0;
  if(tvp_keep == "last") tvp_keep_start += (logvar.n_rows - 1);

  // maybe add storage for tau2 and lambda2, shrinkage hyperparameters of ngpior for factorloadings

  int phi_hyperparameter_size(0);
  if(priorPHI == "DL" ){
    phi_hyperparameter_size += 2*n_groups + 2*n; // a + xi + n(lambda + psi)
  }else if(priorPHI == "GT"){
    phi_hyperparameter_size += 2*n_groups + 2*n; // a+xi + n(lambda + psi)
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
  arma::mat phi_hyperparameter_draws(phi_hyperparameter_size, priorPHI != "normal" ? nsave : 0);

  int u_hyperparameter_size(0);
  if(sigma_type == "cholesky"){
    if(priorU == "DL" ){
      u_hyperparameter_size += 2 + 2*n_U;
    }else if(priorU == "GT"){
      u_hyperparameter_size += 2 + 2*n_U; // a+xi + n(lambda + psi)
    }else if(priorU == "R2D2"){
      u_hyperparameter_size += 4 + 2*n_U; // b + api +xi + zeta + n(theta + psi)
    }else if(priorU == "SSVS"){
      u_hyperparameter_size += 2*n_U;
    }else if(priorU == "HMP"){
      u_hyperparameter_size += 1;
    }else if(priorU == "HS"){
      u_hyperparameter_size += 2+2*n_U;
    }
  }
  arma::mat u_hyperparameter_draws(sigma_type=="cholesky" ? u_hyperparameter_size : 0, (sigma_type=="cholesky" && priorU != "normal") ? nsave : 0);

  //indicator vector needed for DL and R2D2 on U
  arma::uvec ind_U(n_U);
  for(int i=0; i<n_U; ++i){
    ind_U(i)=i;
  }
  //-----------------------------------SAMPLER--------------------------------//
  const int tot = draws + burnin;
  // Initialize progressbar
  //Progress p(tot, progressbar);
  Timer timer;
  timer.step("start");

  for(int rep = 0; rep < tot; rep++){

    // Check for user interrupts
    if (rep % 256 == 0) {
      Rcpp::checkUserInterrupt();
    }

    //----1) Draw PHI (reduced form VAR coefficients)
    if(sigma_type == "cholesky"){
      //arma::mat PHI_old = PHI; // zombie???
      try{
        sample_PHI(PHI, PHI0, Y, X, U, d_sqrt, V_prior, M);
      } catch(...){
        ::Rf_error("Couldn't sample PHI in rep %i.", rep);
      }
    }else if(sigma_type == "factor"){
      sample_PHI_factor(PHI, PHI0, Y, X, logvar.cols(0,M-1), V_prior,
                        armafacload, armafac, huge);
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

            sample_V_i_DL(V_i, PHI_diff(i_ocl), a(j),b(j),c(j), a_vec, a_weight, psi,
                          lambda, xi(j), ind, DL_hyper, norm_consts, GL_tol, DL_plus,
                          c_vec, c_rel_a);

        }else if(priorPHI == "GT"){ // R2D2 is GT with exponential kernel, NG is GT with normal kernel
          sample_V_i_GT(V_i, PHI_diff(i_ocl), psi, lambda, xi(j), a(j), b(j), c(j),
                        ind, GL_tol, GT_priorkernel, GT_vs, norm_consts,
                        a_vec, a_weight, c_vec, GT_hyper, c_rel_a);
        }else if(priorPHI == "SSVS"){

          if(rep > 0.1*burnin || SSVS_hyper){

            sample_V_i_SSVS_beta(V_i, gammas, p_i, PHI_diff(i_ocl), tau_0,
                                 tau_1, SSVS_s_a, SSVS_s_b, SSVS_hyper, ind);

          }

        }else if(priorPHI == "HS"){

          sample_V_i_HS(V_i, PHI_diff(i_ocl), theta_hs, zeta_hs(j), nu, varpi(j), ind);

        }
        j++;
      }

    }else if(priorPHI == "HMP"){

      if(rep > 0.1*burnin){
        sample_V_i_HMP(lambda_1, lambda_2, V_i, s_r_1(0), s_r_1(1), s_r_2(0),
                       s_r_2(1), PHI_diff(i_ocl), V_i_prep(i_ocl), n_ol, n_cl,
                       i_ol, i_cl);
      }
    }

    V_i_long(i_ocl) = V_i;

    arma::mat resid = Y - X*PHI;

    //----3) Draw Sigma_t

    if(sigma_type == "factor"){
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

    }else if(sigma_type == "cholesky"){

      try{
        sample_U(U, resid, V_i_U, d_sqrt);
      }
      catch(...){
        ::Rf_error("Couldn't sample U in rep %i.", rep);
      }

      if(priorU == "DL" || priorU == "R2D2" || priorU == "GT"){

        for (int jj = 0; jj<M; jj++){
          for (int ii = 0; ii<jj; ii++) {
            if(U(ii,jj) == 0) {
              if(R::rbinom( 1, 0.5 )==0){
                U(ii,jj) = L_tol;
              }else{
                U(ii,jj) = -L_tol;
              }
            }else
              if(U(ii,jj) < L_tol && U(ii,jj) > 0){
                U(ii,jj) = L_tol;
              }else if (U(ii,jj) > -L_tol && U(ii,jj) < 0){
                U(ii,jj) = -L_tol;
              }
          }
        }
      }

      //----3b) Draw hyperparameters of hierarchical priors
      u = U(U_upper_indices);
      if(priorU == "DL" ){


        try{
            sample_V_i_DL(V_i_U, u, a_U, b_U, c_U, a_vec_U, a_weight_U,
                        psi_U, lambda_U, xi_U, ind_U, DL_hyper_U, norm_consts_U,
                        GL_tol_U, DL_plus_U, c_vec_U, c_rel_a_U);
        } catch (...) {
          ::Rf_error("Couldn't sample V_i_U (DL prior)  in run %i", rep);

        }

      }else if (priorU == "GT"){

        sample_V_i_GT(V_i_U, u, psi_U, lambda_U, xi_U, a_U, b_U, c_U, ind_U,
                      GL_tol_U, GT_priorkernel_U, GT_vs_U, norm_consts_U,
                      a_vec_U, a_weight_U, c_vec_U, GT_hyper_U, c_rel_a_U);

      }else if(priorU == "SSVS"){

        sample_V_i_SSVS_beta(V_i_U, gammas_U, p_i_U, u, tau_0_U, tau_1_U, SSVS_s_a_U,
                             SSVS_s_b_U, SSVS_hyper_U, ind_U);

      }else if(priorU == "HMP"){

        sample_V_i_U_HMP(lambda_3, V_i_U, s_r_3(0), s_r_3(1), u);
      }else if(priorU == "HS"){

        sample_V_i_HS(V_i_U, u, theta_hs_U, zeta_hs_U, nu_U, varpi_U, ind_U);

      }

      //----3c) Draw elements of D_t
      //        in case of SV use package stochvol
      arma::mat str_resid = resid*U; // structural (orthogonalized) residuals

      for(int j =0; j<M; j++){
        if(heteroscedastic(j) == false){
          double s_p = priorHomoscedastic(j,1) + 0.5*accu(square(str_resid.col(j)));
          double d_i = 1. / R::rgamma(priorHomoscedastic(j,0)+.5*T, 1./s_p);
          d_sqrt.col(j).fill(sqrt(d_i));
          logvar.col(j).fill(log(d_i));
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
          d_sqrt.col(j) = exp(h_j/2);
        }
      }
    }

    //-------Store draws after burnin
    if(rep >= burnin && ( (rep+1-burnin) % thin == 0 )){

      PHI_draws.slice((rep+1-burnin)/thin - 1) = PHI;
      logvar_draws.slice((rep+1-burnin)/thin - 1) = logvar.rows(tvp_keep_start ,logvar.n_rows-1);
      sv_para_draws.slice(((rep+1-burnin)/thin - 1)) = sv_para;
      V_prior_draws.slice((rep+1-burnin)/thin - 1) = V_prior;

      if(priorPHI == "DL" || priorPHI== "GT"){

        phi_hyperparameter_draws(span(0,(n_groups-1)), (rep+1-burnin)/thin - 1) = a;
        phi_hyperparameter_draws(span(n_groups,(2*n_groups-1)), (rep+1-burnin)/thin - 1) = xi;
        phi_hyperparameter_draws(span(2*n_groups,(2*n_groups+n-1)), (rep+1-burnin)/thin - 1) = lambda;
        phi_hyperparameter_draws(span((2*n_groups+n),(phi_hyperparameter_size-1)), (rep+1-burnin)/thin - 1) = psi;

        // if(priorPHI == "GT" ){
        //
        //   phi_hyperparameter_draws((rep+1-burnin)/thin - 1, span(0,(n_groups-1))) = a;
        //   phi_hyperparameter_draws((rep+1-burnin)/thin - 1, span(n_groups,(2*n_groups-1))) = xi;
        //   phi_hyperparameter_draws((rep+1-burnin)/thin - 1, span(2*n_groups,(2*n_groups+n-1))) = lambda;
        //   phi_hyperparameter_draws((rep+1-burnin)/thin - 1, span((2*n_groups+n),(phi_hyperparameter_size-1))) = psi;
        //
        // }

      }else if(priorPHI == "HS"){

        phi_hyperparameter_draws(span(0,(n_groups-1)), (rep+1-burnin)/thin - 1) = zeta_hs;
        phi_hyperparameter_draws(span(n_groups,(2*n_groups-1)), (rep+1-burnin)/thin - 1) = varpi;
        phi_hyperparameter_draws(span(2*n_groups,(2*n_groups+n-1)), (rep+1-burnin)/thin - 1) = theta_hs;
        phi_hyperparameter_draws(span((2*n_groups+n),(phi_hyperparameter_size-1)), (rep+1-burnin)/thin - 1) = nu;

      }else if(priorPHI == "SSVS"){

        phi_hyperparameter_draws(span(0, (n-1)), (rep+1-burnin)/thin - 1) = gammas;
        phi_hyperparameter_draws(span(n, (phi_hyperparameter_size-1)), (rep+1-burnin)/thin - 1) = p_i;

      }else if(priorPHI == "HMP"){
        phi_hyperparameter_draws(0, (rep+1-burnin)/thin - 1) = lambda_1;
        phi_hyperparameter_draws(1, (rep+1-burnin)/thin - 1) = lambda_2;
      }

      if(sigma_type == "cholesky"){

        U_draws.col((rep+1-burnin)/thin - 1) = u;

        if(priorU == "DL" || priorU == "GT"){

          u_hyperparameter_draws(0, (rep+1-burnin)/thin - 1) = a_U;
          u_hyperparameter_draws(1, (rep+1-burnin)/thin - 1) = xi_U;
          u_hyperparameter_draws(span(2,n_U+1), (rep+1-burnin)/thin - 1) = lambda_U;
          u_hyperparameter_draws(span(n_U+2,u_hyperparameter_size-1), (rep+1-burnin)/thin - 1) = psi_U;

          // u_hyperparameter_draws((rep+1-burnin)/thin - 1, 0) = a_U;
          // u_hyperparameter_draws((rep+1-burnin)/thin - 1, span(1,n_U)) = psi_U.as_row();
          // u_hyperparameter_draws((rep+1-burnin)/thin - 1, span(n_U+1,2*n_U)) = lambda_U.as_row();
          // u_hyperparameter_draws((rep+1-burnin)/thin - 1, u_hyperparameter_size-1) = xi_U;
          // else if(priorU == "GT"){
          //
          //   u_hyperparameter_draws((rep+1-burnin)/thin - 1, 0) = a_U;
          //   u_hyperparameter_draws((rep+1-burnin)/thin - 1, 1) = xi_U;
          //   u_hyperparameter_draws((rep+1-burnin)/thin - 1, span(2,n_U+1)) = psi_U;
          //   u_hyperparameter_draws((rep+1-burnin)/thin - 1, span(n_U+2,u_hyperparameter_size-1)) = lambda_U;
          //
          // }

        }else if(priorU == "SSVS"){

          u_hyperparameter_draws(span(0, (n_U-1)), (rep+1-burnin)/thin - 1) = gammas_U;
          u_hyperparameter_draws(span(n_U, (u_hyperparameter_size-1)), (rep+1-burnin)/thin - 1) = p_i_U;

        }else if(priorU == "HMP"){

          u_hyperparameter_draws(0, (rep+1-burnin)/thin - 1) = lambda_3;

        }else if(priorU == "HS"){

          u_hyperparameter_draws(0, (rep+1-burnin)/thin - 1) = zeta_hs_U ;
          u_hyperparameter_draws(1, (rep+1-burnin)/thin - 1) = varpi_U;
          u_hyperparameter_draws(span(2, n_U+1), (rep+1-burnin)/thin - 1) = theta_hs_U;
          u_hyperparameter_draws(span((n_U+2),(u_hyperparameter_size-1)), (rep+1-burnin)/thin - 1)= nu_U;

        }
      }else if(sigma_type == "factor"){
        facload_draws.slice((rep+1-burnin)/thin - 1) = armafacload;
        fac_draws.slice((rep+1-burnin)/thin - 1) = armafac.cols(tvp_keep_start ,armafac.n_cols-1);
      }
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
    Named("U") = U_draws,
    Named("logvar") = logvar_draws,
    Named("sv_para") = sv_para_draws,
    Named("phi_hyperparameter") = phi_hyperparameter_draws,
    Named("u_hyperparameter") = u_hyperparameter_draws,
    Named("bench") = time,
    Named("V_prior") = V_prior_draws,
    Named("facload") = facload_draws,
    Named("fac") = fac_draws
  );

  return out;
}
