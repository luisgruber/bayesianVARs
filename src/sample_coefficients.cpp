#include <RcppArmadilloExtensions/sample.h>
#include "univariate_regression.h"
#include "utilities_cpp.h"

using namespace Rcpp;
using namespace arma;

/*see https://stackoverflow.com/a/23487774
 * arma::vec rnormvec(p);
 * rnormvec = Rcpp::rnorm(p); // !bad!, creation (and memory allocation) of temporary NumericVector and arma::vec
 * rnormvec.imbue(R::norm_rand); // !good!, .imbue() is similar to .fill() but handles functors or lambda functions. R::norm_rand is the underlying function of wrappers like Rcpp::rnorm and R::rnorm.
 * std::generate(vec.begin(),vec.end(), ::norm_rand); // !good!, STL style
 * for(int j = 0; j < p; ++j){ // !okayish!. in terms of memory allocation as good as arma::vec.imbue(::norm_rand)/std::generate(vec.begin(),vec.end(), ::norm_rand), but slightly slower
 *    rnormvec(j) = R::rnorm(0,1);
 * }
*/

double do_rgig(double lambda, double chi, double psi) {

    SEXP (*fun)(int, double, double, double) = NULL;
    if (!fun) fun = (SEXP(*)(int, double, double, double)) R_GetCCallable("GIGrvg", "do_rgig");

    return as<double>(fun(1, lambda, chi, psi));
}

//Kastner and Huber 2019 JoF equation per equation
void sample_PHI_factor(arma::mat& PHI, const arma::mat& PHI_prior,
                       const arma::mat& Y, const arma::mat& X,
                       const arma::mat& logvaridi, const arma::mat& V_prior,
                       const arma::mat& facload, const arma::mat& fac,
                       const bool& huge, const double floor){

  const int M = Y.n_cols;
  const int r = fac.n_rows;

  const arma::mat Y_hat = (r > 0) ? Y - arma::trans(facload*fac) : Y; // conditioning on factors leads to independent equations
  
  const arma::mat normalizer = exp(-logvaridi/2.);

  for(int j = 0; j < M; j++){

    arma::mat X_norm = X.each_col() % normalizer.col(j);
    arma::vec y_norm = Y_hat.col(j) % normalizer.col(j);
    arma::vec PHI_j = PHI.unsafe_col(j);
    arma::vec PHI_prior_j = PHI_prior.unsafe_col(j);
    arma::vec V_prior_j = V_prior.unsafe_col(j);

      if(huge) {
        try {
          PHI_j = sample_coef_blm_known1var_large_p_small_n(y_norm, X_norm,
                                                            PHI_prior_j,
                                                            V_prior_j);  
        } catch(const std::exception& e) {
          Rcpp::stop("sample_coef_blm_known1var_large_p_small_n() failed in %ith eqation: %s", j+1, e.what());
        } catch (...) {
          Rprintf("\nsample_coef_blm_known1var_large_p_small_n() failed in %ith eqation. Rethrowing exception:", j+1);
          throw;
        }
      } else {
        try {
          PHI_j = sample_coef_blm_known1var(X_norm.t()*y_norm, X_norm.t()*X_norm,
                                            PHI_prior_j,
                                            V_prior_j);
        } catch (...) {
          try {
            PHI_j = sample_coef_blm_known1var_qr(y_norm, X_norm,
                                                 PHI_prior_j,
                                                 V_prior_j);  
          } catch(const std::exception& e) {
            Rcpp::stop("sample_coef_blm_known1var_qr() failed in %ith eqation: %s", j+1, e.what());
          } catch (...) {
            Rprintf("\nsample_coef_blm_known1var_qr() failed in %ith eqation. Rethrowing exception:", j+1);
            throw;
          }
        }
      }
    if(floor>0) {
      for (arma::uword i=0; i<PHI_j.n_elem; ++i) {
        abs_floor2(PHI_j(i), PHI_prior_j(i), floor);
      }
    }
  }

}

// sample_PHI samples VAR coefficients using the corrected triangular
// algorithm as in Carrier, Chan, Clark & Marcellino (2021)
void sample_PHI(arma::mat& PHI, const arma::mat& PHI_prior, const arma::mat& Y,
                const arma::mat& X, const arma::mat& U, const arma::mat& d_sqrt,
                const arma::mat& V_prior, const int M, const double floor) {

  for(int j = 0; j < M; j++){

    arma::mat PHI_0 = PHI;
    PHI_0.col(j).zeros();
    arma::mat normalizer_mat = d_sqrt.cols(j, M-1);
    arma::vec normalizer(normalizer_mat.begin(), normalizer_mat.size(), false);
    arma::mat Y_new_mat = (Y - X * PHI_0) * U.cols(j, M-1);
    arma::vec Y_new(Y_new_mat.begin(), Y_new_mat.size(), false);
    Y_new /= normalizer;
    arma::mat X_new_tmp = arma::kron(U( span(j,j), span(j, M-1)).t(), X);
    arma::mat X_new = X_new_tmp.each_col() / normalizer;

    arma::vec PHI_j = PHI.unsafe_col(j);
    arma::vec PHI_prior_j = PHI_prior.unsafe_col(j);
    arma::vec V_prior_j = V_prior.unsafe_col(j);
    // generate posterior draw
    try {
      PHI_j = sample_coef_blm_known1var(X_new.t()*Y_new, X_new.t()*X_new,
                                        PHI_prior_j,
                                        V_prior_j);
      
    } catch (...) {
      try {
        PHI_j = sample_coef_blm_known1var_qr(Y_new, X_new,
                                             PHI_prior_j,
                                             V_prior_j);  
      } catch(const std::exception& e) {
        Rcpp::stop("sample_coef_blm_known1var_qr() failed in %ith eqation:\n%s", j+1, e.what());
      } catch (...) {
        Rprintf("\nsample_coef_blm_known1var_qr() failed in %ith eqation. Rethrowing exception:", j+1);
        throw;
      }
    }  
      
    if(floor>0) {
      for (arma::uword i=0; i<PHI_j.n_elem; ++i) {
        abs_floor2(PHI_j(i), PHI_prior_j(i), floor);
      }
    }
  }
}

// [[Rcpp::export]]
arma::mat sample_PHI_cholesky(const arma::mat& PHI, const arma::mat& PHI_prior,
                              const arma::mat& Y, const arma::mat& X,
                              const arma::mat& U, const arma::mat& d_sqrt,
                              const arma::mat& V_prior,
                              const double floor){
  const int M = Y.n_cols;
  arma::mat PHI_ret = PHI;
  sample_PHI(PHI_ret, PHI_prior, Y, X, U, d_sqrt, V_prior, M, floor);
  return(PHI_ret);
}


void sample_U(arma::mat& U, const arma::mat& Ytilde, const arma::vec& V_i, const arma::mat& d_sqrt,
              const double floor) {

  int M = Ytilde.n_cols;

  int ind = 0;
  // arma::vec prior_mean(1); // proxy needed for univariate_regression_update
  for(int i = 1; i < M; i++){

    // create matrix with covariables (residuals of preceeding equations)
    arma::mat Z = Ytilde.cols(0, i-1);
    Z = -1.0 * (Z.each_col() / d_sqrt.col(i));
    arma::vec prior_mean(i, arma::fill::zeros); // proxy needed for univariate_regression_update
    // create 'dependent' variable
    arma::vec c = Ytilde.col(i) / d_sqrt.col(i);
    // update free off-diagonal elements in column i
    arma::vec U_j = U.unsafe_col(i);
    try{
      arma::vec tmp;
      try {
        tmp = sample_coef_blm_known1var(Z.t()*c, Z.t()*Z, prior_mean, V_i.subvec(ind, ind+i-1));
        // V_i probably not passed by reference in that case, because subvec is a subview
      } catch (...) {
        tmp = sample_coef_blm_known1var_qr(c, Z, prior_mean, V_i.subvec(ind, ind+i-1));
      }
      if(floor>0) {
        for(arma::uword k=0; k<tmp.n_elem; ++k){
          abs_floor(tmp(k), floor);
        }  
      }
      U_j.subvec(0,i-1) = tmp;
    } catch (const std::exception& e) {
      Rcpp::stop("sample_coef_blm_known1var() failed in equation %i: %s", i+1, e.what());
    } catch (...) {
      Rprintf("\nsample_coef_blm_known1var() failed in equation %i. Rethrowing exception:", i+1);
      throw;
    }
    
    ind = ind+i;

  }

}

void sample_a_c(double& a, double& c, const arma::vec& a_vec, 
                const arma::vec& a_weight, const arma::vec& norm_consts,
                const arma::vec& c_vec, const bool c_rel_a, const double xi,
                const double b, const arma::vec& lambda, const int n){
  
  arma::vec logprobs = log(a_weight) + n*(a_vec*log(xi) - norm_consts) +
    (a_vec-1.0)*arma::accu(log(lambda)) ;
  if(!c_rel_a){
    logprobs += (-1.0)*2.0*xi*c/a_vec- b*log(a_vec);
  }
  
  arma::vec w = exp(logprobs - logprobs.max());
  arma::vec weights = w/sum(w);
  const arma::uword gridlength = a_vec.n_elem;
  const arma::uvec idx = arma::linspace<arma::uvec>(0, gridlength-1, gridlength);
  //arma::vec which = Rcpp::RcppArmadillo::sample(a_vec, 1, true, weights); // if(c_rel_a) this would require calling arma::find later
  arma::uvec which = Rcpp::RcppArmadillo::sample(idx, 1, true, weights);
  a = a_vec(which(0));
  if(c_rel_a){
    c = c_vec(which(0));
  }
}

void sample_V_i_GT(arma::vec& V_i, const arma::vec& coefs, arma::vec& psi,
                   arma::vec& lambda, double& xi, double& a, const double b,
                   double& c, const arma::uvec& ind,
                   const std::string priorkernel,
                   const double vs, const arma::vec& norm_consts,
                   const arma::vec& a_vec, const arma::vec& a_weight,
                   const arma::vec& c_vec, const bool hyper, const bool c_rel_a){

  const int n = ind.size();
  const arma::vec coefs_squared = arma::square(coefs);
  for(arma::uvec::const_iterator it = ind.begin(); it != ind.end(); ++it){
    // For R2d2exp vs must be 1/2, for NG 1
    double lambdascale = vs;
    // sample auxiliary scaling parameters
    if(priorkernel == "exponential"){
      psi(*it) = 1./do_rgig(-0.5, 1.0, coefs_squared(*it)/(lambda(*it)*vs));
      lambdascale *= psi(*it);
    }
    // sample local scales
    lambda(*it) = do_rgig(a - 0.5, coefs_squared(*it)/lambdascale, 2.0*xi); 
  }
  double zeta = arma::accu(lambda(ind));
  //sample global scale
  xi = R::rgamma(b+n*a, 1./(2.0*c/a + zeta));

  V_i(ind) = lambda(ind)*vs;
  if(priorkernel == "exponential"){
    V_i(ind) %= psi(ind);
  }

  if(hyper){
    sample_a_c(a, c, a_vec, a_weight, norm_consts, c_vec, c_rel_a, xi, b, lambda(ind),
               n);
  }
}

void sample_V_i_HS(arma::vec& V_i, const arma::vec& coefs, arma::vec& lambda,
                   double& xi, arma::vec& psi, double& varpi, const arma::uvec& ind){

  int n = ind.size();
  for(arma::uvec::const_iterator it = ind.begin(); it != ind.end(); ++it){
    // R::rgamma uses scale!!!!
    lambda(*it) = 1./(R::rgamma(1.0, 1./( 1.0/psi(*it) + (coefs(*it)*coefs(*it))/(2.0*xi) )));
    psi(*it) = 1./(R::rgamma(1.0, 1./(1.0 + 1./lambda(*it)) ));
  }
  xi = 1./(R::rgamma((n+1.0)/2., 1./(1./varpi + 0.5*arma::accu(square(coefs(ind))/lambda(ind))) ));
  varpi = 1./(R::rgamma(1.0, 1./(1.0 + 1.0/xi) ));
  V_i(ind) = lambda(ind)*xi;

}

void sample_V_i_DL(arma::vec& V_i, const arma::vec& coefs, double& a ,
                   const double b, double& c,
                   const arma::vec& a_vec, const arma::vec& a_weight,
                   arma::vec& psi, arma::vec& lambda, double& xi, const arma::uvec& ind,
                   const bool hyper,const arma::vec& norm_consts,
                   const bool DL_plus,
                   const arma::vec& c_vec, const bool c_rel_a){ //, bool hyper

  const int n = ind.size();
  arma::vec coefs_abs = arma::abs(coefs);

  for(arma::uvec::const_iterator it = ind.begin(); it != ind.end(); ++it){
    //important: joint update of p(lambda,psi|...): first lambda, than psi!!! (is wrong in Bhattacharya et al 2015!)
    lambda(*it) = do_rgig(a-1., 2.0*coefs_abs(*it), 2.0*xi);
    psi(*it) =  do_rgig(0.5, (coefs(*it)*coefs(*it)) /
      (lambda(*it)*lambda(*it)), 1.0 );
  }
  
  if(DL_plus){
    const double zeta = arma::accu(lambda(ind));
    xi = R::rgamma(b+n*a, 1./(2.0*c/a + zeta));
  }

  V_i(ind) = psi(ind) % square(lambda(ind));

  if(hyper){
    sample_a_c(a, c, a_vec, a_weight, norm_consts, c_vec, c_rel_a, xi, b, lambda(ind),
               n);
  }

}

void sample_V_i_SSVS_beta(arma::vec& V_i, arma::vec& lambda, double& a,
                     const arma::vec& coeffs, const arma::vec& tau_0,
                     const arma::vec& tau_1, const double s_a, const double s_b,
                     const bool hyper, const arma::uvec& ind){

  int n = ind.size();

  for(arma::uvec::const_iterator it = ind.begin(); it != ind.end(); ++it){
    // taking squares
    const double phi2_i = coeffs(*it)*coeffs(*it);
    const double tau2_1_i = tau_1(*it)*tau_1(*it);
    const double tau2_0_i = tau_0(*it)*tau_0(*it);
    // Compute conditional posterior inclusion probability
    // taking logs is much stabler!!!
    const double u_i1 = -std::log(tau_1(*it)) - 0.5 * phi2_i/tau2_1_i + std::log(a);
    const double u_i2 = -std::log(tau_0(*it)) - 0.5 * phi2_i/tau2_0_i + std::log(1.0 - a);
    const double logdif = u_i2 - u_i1;
    const double gst = 1.0/(1.0 + exp(logdif)); // == exp(u_i1)/(exp(u_i2) + exp(u_i1))
    // Draw lambda
    lambda(*it) = R::rbinom(1,gst);

    // Compute prior variances
    if(lambda(*it)==1){
      V_i(*it) = tau2_1_i;
    }else{
      V_i(*it) = tau2_0_i;
    }
  }

  if(hyper){
    int incl = accu(lambda(ind));
    a = R::rbeta(s_a + incl, s_b + n - incl);
  }
}

void sample_V_i_HMP(double& lambda_1, double& lambda_2, arma::vec& V_i, const double& s1,
                    const double& r1, const double& s2, const double& r2,
                    const arma::vec& PHI_diff, const arma::vec& V_i_prep,
                    const int& n_ol, const int& n_cl, const arma::uvec& i_ol,
                    const arma::uvec& i_cl){

  lambda_1 = do_rgig(s1 - 0.5*n_ol, sum(square(PHI_diff(i_ol))/V_i_prep(i_ol)),
                     2*r1 );
  lambda_2 = do_rgig(s2 - 0.5*n_cl, sum(square(PHI_diff(i_cl))/V_i_prep(i_cl)), 2.0*r2 );

  V_i(i_ol) = lambda_1 * V_i_prep(i_ol);
  V_i(i_cl) = lambda_2 * V_i_prep(i_cl);

}

void sample_V_i_U_HMP(double& lambda_3, arma::vec& V_i_U, const double& s1,
                    const double& r1, const arma::vec& u){

  int n = u.size();
  lambda_3 = do_rgig(s1 - 0.5*n, sum(square(u)), 2.0*r1 );
  V_i_U.fill(lambda_3);

}
