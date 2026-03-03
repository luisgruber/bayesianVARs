#include <RcppArmadillo.h>

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

// univariate_regression_update updates the coefficientes of an univariate
// linear regression with known variance equaling one, i.e. both the dependent
// and the independent variables have to be scaled with the regression variance.
void univariate_regression_update(arma::vec& coef_draw,
                    const arma::vec& prior_mean,
                    const arma::vec& prior_variance,
                    const arma::vec& y_norm,
                    const arma::mat& X_norm,
                    const bool& typeU,
                    const bool& huge){

  const int K = prior_variance.size(); // coef_draw.size() will not work if typeU==true
  const int T = y_norm.size();

  arma::vec coef_draw_tmp(K);

  if(huge){

    arma::vec u(K); u.imbue(R::norm_rand); // 
    u %= arma::sqrt(prior_variance);
    if(!typeU){
      u += prior_mean;
    }
    arma::vec delta(T); delta.imbue(R::norm_rand);

    arma::vec v = X_norm*u + delta;
    const arma::mat X_norm_pv = X_norm.each_row() % prior_variance.as_row();
    arma::vec w = arma::solve(X_norm_pv*X_norm.t() + arma::eye(T,T), y_norm - v);
    coef_draw_tmp = u + X_norm_pv.t()*w;

  }else{
    arma::mat Omega_post = X_norm.t()*X_norm;
    Omega_post.diag() += 1./prior_variance;

    arma::mat Omega_post_chol =  chol(Omega_post, "upper");
    arma::mat Omega_post_chol_inv = arma::inv(trimatu(Omega_post_chol));

    arma::vec coef_post = X_norm.t()*y_norm;
    if(!typeU) { // U has by construction prior_mean of zero
      coef_post += prior_mean / prior_variance;
    }
    coef_post = Omega_post_chol_inv * Omega_post_chol_inv.t() * coef_post;

    arma::colvec rnormvec(K);
    rnormvec.imbue(R::norm_rand);

    coef_draw_tmp = coef_post + Omega_post_chol_inv*rnormvec;
  }

  if(typeU){
    //pad draw with zeros (and one 1 for the diagonal element)
    coef_draw.zeros();
    coef_draw.subvec(0,K-1) = coef_draw_tmp;
    coef_draw(K) = 1.;
  }else{
    coef_draw = coef_draw_tmp;
  }
}

//Kastner and Huber 2019 JoF equation per equation
void sample_PHI_factor(arma::mat& PHI, const arma::mat& PHI_prior,
                       const arma::mat& Y, const arma::mat& X,
                       const arma::mat& logvaridi, const arma::mat& V_prior,
                       const arma::mat& facload, const arma::mat& fac,
                       const bool& huge){

  const int M = Y.n_cols;
  const int r = fac.n_rows;

  const arma::mat Y_hat = (r > 0) ? Y - arma::trans(facload*fac) : Y; // conditioning on factors leads to independent equations
  
  const arma::mat normalizer = exp(-logvaridi/2.);

  for(int j = 0; j < M; j++){

    arma::mat X_norm = X.each_col() % normalizer.col(j);
    arma::vec y_norm = Y_hat.col(j) % normalizer.col(j);
    arma::vec PHI_j = PHI.unsafe_col(j);
    try{
      univariate_regression_update(PHI_j, PHI_prior.col(j), V_prior.col(j),
                                   y_norm, X_norm, false, huge);
    } catch(const std::exception& e) {
      Rcpp::stop("univariate_regression_update() failed in %ith eqation: %s", j+1, e.what());
    } catch (...) {
      Rprintf("\nunivariate_regression_update() failed in %ith eqation. Rethrowing exception:", j+1);
      throw;
    }

  }

}

// sample_PHI samples VAR coefficients using the corrected triangular
// algorithm as in Carrier, Chan, Clark & Marcellino (2021)
void sample_PHI(arma::mat& PHI, const arma::mat PHI_prior, const arma::mat Y,
                const arma::mat X, const arma::mat U, const arma::mat d_sqrt,
                const arma::mat V_prior, const int M) {

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
    // generate posterior draw
    try{
      univariate_regression_update(PHI_j, PHI_prior.col(j), V_prior.col(j),
                                   Y_new, X_new, false, false);
    } catch(const std::exception& e) {
      Rcpp::stop("univariate_regression_update() failed in %ith eqation:\n%s", j+1, e.what());
    } catch (...) {
      Rprintf("\nunivariate_regression_update() failed in %ith eqation. Rethrowing exception:", j+1);
      throw;
    }
  }
}

// [[Rcpp::export]]
arma::mat sample_PHI_cholesky(const arma::mat PHI, const arma::mat& PHI_prior,
                              const arma::mat& Y, const arma::mat& X,
                              const arma::mat& U, const arma::mat& d_sqrt,
                              const arma::mat& V_prior){
  const int M = Y.n_cols;
  arma::mat PHI_ret = PHI;
  sample_PHI(PHI_ret, PHI_prior, Y, X, U, d_sqrt, V_prior, M);
  return(PHI_ret);
}


void sample_U(arma::mat& U, const arma::mat& Ytilde, const arma::vec& V_i, const arma::mat& d_sqrt) {

  int M = Ytilde.n_cols;

  int ind = 0;
  arma::vec prior_mean(1); // proxy needed for univariate_regression_update
  for(int i = 1; i < M; i++){

    // create matrix with covariables (residuals of preceeding equations)
    arma::mat Z = Ytilde.cols(0, i-1);
    Z = -1 * (Z.each_col() / d_sqrt.col(i));
    // create 'dependent' variable
    arma::vec c = Ytilde.col(i) / d_sqrt.col(i);
    // update free off-diagonal elements in column i
    arma::vec U_j = U.unsafe_col(i);
    try{
      univariate_regression_update(U_j, prior_mean, V_i.subvec(ind, ind+i-1), c, Z, true, false);  
    } catch (const std::exception& e) {
      Rcpp::stop("univariate_regression_update() failed in equation %i: %s", i+1, e.what());
    } catch (...) {
      Rprintf("\nunivariate_regression_update() failed in equation %i. Rethrowing exception:", i+1);
      throw;
    }
    
    ind = ind+i;

  }

}

void sample_V_i_GT(arma::vec& V_i, const arma::vec coefs, arma::vec& psi,
                   arma::vec& lambda, double& xi, double& a, const double b,
                   double& c, arma::uvec ind, const double tol,
                   const std::string priorkernel,
                   const double vs, const arma::vec norm_consts,
                   const arma::vec a_vec, const arma::vec a_weight,
                   const arma::vec c_vec, const bool hyper, const bool c_rel_a){

  const int n = ind.size();
  const int gridlength = a_vec.size();
  arma::vec logprobs(gridlength);
  const arma::vec coefs_squared = arma::square(coefs);
  arma::uvec::iterator it;
  for(it = ind.begin(); it != ind.end(); ++it){
    // sample auxiliary scaling parameters
    if(priorkernel == "exponential"){
      psi(*it) = 1./do_rgig(-0.5, 1, coefs_squared(*it)/(lambda(*it)*vs)); //For R2d2exp vs must be 1/2
    }
    // sample local scales
    lambda(*it) = do_rgig(a - 0.5, coefs_squared(*it)/(psi(*it)*vs), 2*xi); // For R2d2exp vs must be 1/2, for NG psi must be 1
  }
  double zeta = arma::accu(lambda(ind));
  if(tol>0){
    arma::vec theta = lambda(ind)/zeta;
    arma::uvec theta_ind = find( theta < tol);
    theta(theta_ind) = arma::vec(theta_ind.size(), fill::value(tol));
    lambda(ind) = theta * zeta;
  }
  //sample global scale
  xi = R::rgamma(b+n*a, 1./(2*c/a + zeta));

  if(priorkernel == "exponential"){
    V_i(ind) = psi(ind)%lambda(ind)*vs;
  }else if(priorkernel == "normal"){
    V_i(ind) = lambda(ind)*vs;
  }

  if(hyper){

    logprobs = log(a_weight) +  n*(a_vec*log(xi) - norm_consts) +
      (a_vec-1)*arma::accu(log(lambda(ind))) ;
    if(!c_rel_a){//if c is a fixed proportion or multiple of a, the following terms cancel:
      logprobs += (-1)*2*xi*c/a_vec- b*log(a_vec);
    }
    ///
    arma::vec w = exp(logprobs - logprobs.max());
    arma::vec weights = w/sum(w);
    arma::ivec iv(gridlength);
    R::rmultinom(1, weights.begin(), gridlength, iv.begin());
    arma::uvec i = arma::find(iv == 1,1); // reports only the first value that meets the condition (by construction there is only one 1)
    a = a_vec(i(0));
    if(c_rel_a){
      c = c_vec(i(0));
    }
  }
}

void sample_V_i_HS(arma::vec& V_i, const arma::vec coefs, arma::vec& theta,
                   double& zeta, arma::vec& nu, double& varpi ,arma::uvec ind){

  int n = ind.size();
  arma::uvec::iterator it;
  for(it = ind.begin(); it != ind.end(); ++it){
    // R::rgamma uses scale!!!!
    theta(*it) = 1./(R::rgamma(1, 1./( 1/nu(*it) + (coefs(*it)*coefs(*it))/(2*zeta) )));
    nu(*it) = 1./(R::rgamma(1, 1./(1 + 1./theta(*it)) ));
  }
  zeta = 1./(R::rgamma((n+1)/2., 1./(1./varpi + 0.5*arma::accu(square(coefs(ind))/theta(ind))) ));
  varpi = 1./(R::rgamma(1, 1./(1 + 1/zeta) ));
  V_i(ind) = theta(ind)*zeta;

}

void sample_V_i_DL(arma::vec& V_i, const arma::vec coefs, double& a ,
                   const double b, double& c,
                   const arma::vec a_vec, const arma::vec a_weight,
                   arma::vec& psi, arma::vec& lambda, double& xi, arma::uvec ind,
                   const bool hyper,const arma::vec norm_consts,
                   const double tol, const bool DL_plus,
                   const arma::vec c_vec, const bool c_rel_a){ //, bool hyper

  const int n = ind.size();
  arma::vec coefs_abs = arma::abs(coefs);
  int gridlength = a_vec.size();
  arma::vec logprobs(gridlength);

  arma::uvec::iterator it;
  for(it = ind.begin(); it != ind.end(); ++it){
    //important: joint update of p(lambda,psi|...): first lambda, than psi!!! (is wrong in Bhattacharya et al 2015!)
    lambda(*it) = do_rgig(a-1., 2*coefs_abs(*it), 2*xi);
    psi(*it) =  do_rgig(0.5, (coefs(*it)*coefs(*it)) /
      (lambda(*it)*lambda(*it)), 1 );
  }
  double zeta = arma::accu(lambda(ind));
  if(DL_plus){
    xi = R::rgamma(b+n*a, 1./(2*c/a + zeta));
  }

  if(tol>0){
    arma::vec theta = lambda(ind) / zeta;
    uvec th_f = find(theta<tol);
    theta(th_f) = vec(th_f.size(), fill::value(tol));
    lambda(ind) = theta * zeta;
    }

  V_i(ind) = psi(ind) % square(lambda(ind));

  if(hyper){

    logprobs = log(a_weight) + n*(a_vec*log(xi) - norm_consts) +
      (a_vec-1)*arma::accu(log(lambda(ind))) ;
    if(!c_rel_a){
      logprobs += (-1)*2*xi*c/a_vec- b*log(a_vec);
    }

    arma::vec w = exp(logprobs - logprobs.max());
    arma::vec weights = w/sum(w);
    arma::ivec iv(gridlength);
    R::rmultinom(1, weights.begin(), gridlength, iv.begin());
    arma::uvec i = arma::find(iv == 1,1); // reports only the first value that meets the condition (by construction there is only one 1)
    a = a_vec(i(0));
    if(c_rel_a){
      c = c_vec(i(0));
    }
  }

}

void sample_V_i_SSVS_beta(arma::vec& V_i, arma::vec& gammas, arma::vec& p_i,
                     const arma::vec coeffs, const arma::vec tau_0,
                     const arma::vec tau_1, const double s_a, const double s_b,
                     const bool hyper, arma::uvec ind){

  int n = ind.size();
  // Compute conditional posterior inclusion probability
  // taking logs is much stabler!!!
  arma::vec u_i1 = -log(tau_1) - 0.5 * square((coeffs)/tau_1) + log(p_i);
  arma::vec u_i2 = -log(tau_0) - 0.5 * square((coeffs)/tau_0) + log(1 - p_i);
  arma::vec logdif = u_i2 - u_i1;
  arma::vec gst = 1/(1 + exp(logdif)); // == exp(u_i1)/(exp(u_i2) + exp(u_i1))

  arma::uvec::iterator it;
  //double j = 0;
  for(it = ind.begin(); it != ind.end(); ++it){

    // Draw gammas
    gammas(*it) = R::rbinom(1,gst(*it));

    // Compute prior variances
    if(gammas(*it)==1){
      V_i(*it) = tau_1(*it)*tau_1(*it);
    }else{
      V_i(*it) = tau_0(*it)*tau_0(*it);
    }

  }

  if(hyper){
    int incl = accu(gammas(ind));
    double p = R::rbeta(s_a + incl, s_b + n - incl);
    //Rcout << "The value of n : " << n << "\n";
    p_i(ind).fill(p);
  }
}

void sample_V_i_HMP(double& lambda_1, double& lambda_2, arma::vec& V_i, const double& s1,
                    const double& r1, const double& s2, const double& r2,
                    const arma::vec& PHI_diff, const arma::vec& V_i_prep,
                    const int& n_ol, const int& n_cl, const arma::uvec& i_ol,
                    const arma::uvec& i_cl){

  lambda_1 = do_rgig(s1 - n_ol/2, sum(square(PHI_diff(i_ol))/V_i_prep(i_ol)),
                     2*r1 );
  lambda_2 = do_rgig(s2 - n_cl/2, sum(square(PHI_diff(i_cl))/V_i_prep(i_cl)), 2*r2 );

  V_i(i_ol) = lambda_1 * V_i_prep(i_ol);
  V_i(i_cl) = lambda_2 * V_i_prep(i_cl);

}

void sample_V_i_U_HMP(double& lambda_3, arma::vec& V_i_U, const double& s1,
                    const double& r1, const arma::vec& u){

  int n = u.size();
  lambda_3 = do_rgig(s1 - n/2, sum(square(u)), 2*r1 );
  V_i_U.fill(lambda_3);

}
