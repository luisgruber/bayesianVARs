//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

arma::vec mvrnorm1(arma::vec& mu, arma::mat& Sigma, double tol = 1e-06){

  int p = mu.size();
  arma::rowvec rnormvec(p);
  for(int j = 0; j < p; ++j){
    rnormvec(j) = R::rnorm(0,1);
  }

  arma::mat U(p,p);
  bool chol_succes = arma::chol(U, Sigma, "upper" );
  arma::rowvec X(p);
  if(chol_succes){

    X = mu.t() + rnormvec * U;

  }else{
    // imitate R's eigen(x, symmetric = TRUE)
    arma::mat Sigma_sym = symmatl(Sigma);
    arma::vec eigval; // arma stores in ascending order!!!
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, Sigma_sym);
    // check positive definiteness
    for(int j=0; j<p; ++j ){
      if(eigval(j) <= -tol*abs(eigval((p-1.)))){
        Rcpp::stop("'Sigma' is not positive definite");
      }
    }

    arma::colvec pmax(p);
    arma::mat pmat(p,2, arma::fill::zeros);
    pmat.col(0) = eigval;
    pmax = max(pmat, 1);

    arma::mat rmat = diagmat(sqrt(pmax));

    X = mu.t() +  rnormvec * eigvec * rmat;
  }

  return X.t();
}

// Import rgig from R package GIGrvg
double do_rgig1(double lambda, double chi, double psi) {

  if ( !(R_FINITE(lambda) && R_FINITE(chi) && R_FINITE(psi)) ||
       (chi <  0. || psi < 0)      ||
       (chi == 0. && lambda <= 0.) ||
       (psi == 0. && lambda >= 0.) ) {
    Rcpp::stop("invalid parameters for GIG distribution: lambda=%g, chi=%g, psi=%g",
               lambda, chi, psi);
  }

  double res;
  // circumvent GIGrvg in these cases
  if ((chi < (10 * DBL_EPSILON))) {
    /* special cases which are basically Gamma and Inverse Gamma distribution */
    if (lambda > 0.0) {
      res = R::rgamma(lambda, 2.0/psi);
    }
  else {
    res = 1.0/R::rgamma(-lambda, 2.0/chi); // fixed
  }
} else if ((psi < (10 * DBL_EPSILON)) ) {
  /* special cases which are basically Gamma and Inverse Gamma distribution */
  if (lambda > 0.0) {
    res = R::rgamma(lambda, 2.0/psi);  // fixed
  } else {
    res = 1.0/R::rgamma(-lambda, 2.0/chi); // fixed
  }
} else {
    SEXP (*fun)(int, double, double, double) = NULL;
    if (!fun) fun = (SEXP(*)(int, double, double, double)) R_GetCCallable("GIGrvg", "do_rgig");

    res = as<double>(fun(1, lambda, chi, psi));
  }
  return res;
}


// draw_PHI samples VAR coefficients using the corrected triangular
// algorithm as in Carrier, Chan, Clark & Marcellino (2021)
// Rcpp::export --> for use within a sampler written in R

// [[Rcpp::export]]
arma::mat draw_PHI(arma::mat& PHI, arma::mat& PHI_prior, arma::mat& Y, arma::mat& X,
                   arma::mat& L, arma::mat& d, arma::vec& V_i, int& K, int& M) {

  arma::mat V_prior = reshape(V_i, K, M);

  for(int i = 0; i < M; i++){

    arma::mat V_p_inv = diagmat(1/V_prior.col(i));

    arma::mat PHI_0 = PHI;
    PHI_0.col(i).zeros();
    arma::vec normalizer = vectorise( pow(d.cols(i, M-1), 0.5)  );
    arma::vec Y_new = vectorise( (Y - X * PHI_0) * L.cols(i, M-1) ) / normalizer;
    arma::mat X_new_tmp = arma::kron(L( span(i,i), span(i, M-1)).t(), X);
    arma::mat X_new = X_new_tmp.each_col() / normalizer;

    // SL
//    for(int j = 0; j<K; j++){
//      //gammas
//      arma::vec Y_new_star = Y_new - X_new * PHI.col(i) + X_new.col(j)*PHI(j,i);
//      double V_post = 1./(1/V_prior(j,i) + trans(X_new.col(j))*X_new.col(j));
//      double phi_j_post = V_post*(1/V_prior(j,i)*PHI_prior(j,i) + trans(X_new.col(j))*Y_new_star);
//      PHI(j,i) = R::rnorm(phi_j_post,V_post);
//    }
    // Compute V_post
    //arma::mat V_post = ( V_p_inv + X_new.t()*X_new ).i(); //unstable
    // compute inverse via the QR decomposition (similar to chol2inv in R)
    arma::mat V_post_tmp = chol(V_p_inv + X_new.t()*X_new, "upper");
    mat Q;
    mat R;
    qr(Q, R, V_post_tmp);
    mat S = inv(trimatu(R));
    mat V_post= S * S.t();

    arma::colvec phi_post = V_post * (V_p_inv*PHI_prior.col(i) + X_new.t()*Y_new);

    PHI.col(i) = mvrnorm1(phi_post, V_post, 1e-06);

  }

  return(PHI);
}


// draw_PHI samples VAR coefficients using the corrected triangular
// algorithm as in Carrier, Chan, Clark & Marcellino (2021)
// for use within smapler written in Rcpp
void sample_PHI(arma::mat& PHI, const arma::mat& PHI_prior, const arma::mat& Y,
                const arma::mat& X, const arma::mat& L, const arma::mat& d_sqrt,
                const arma::mat& V_prior, const int& K, const int& M, bool subs) {

  for(int i = 0; i < M; i++){

    arma::mat V_prior_inv = diagmat(1/V_prior.col(i));

    arma::mat PHI_0 = PHI;
    PHI_0.col(i).zeros();
    arma::vec normalizer = vectorise( d_sqrt.cols(i, M-1)  );
    arma::vec Y_new = vectorise( (Y - X * PHI_0) * L.cols(i, M-1) ) / normalizer;
    arma::mat X_new_tmp = arma::kron(L( span(i,i), span(i, M-1)).t(), X);
    arma::mat X_new = X_new_tmp.each_col() / normalizer;

    // 'subs' avoids the direct computation of the inverse needed for V_post
    // maybe in very large dimensions a little bit more efficient
    // not thouroughly tested, hence not used in the sampler
    arma::mat V_post_inv_chol = chol(V_prior_inv + X_new.t()*X_new, "upper");
    if(subs){
      arma::vec phi_post_prep = V_prior_inv*PHI_prior.col(i) + X_new.t()*Y_new;
      arma::colvec phi_post_tmp = arma::solve(trimatl(V_post_inv_chol.t()), phi_post_prep);
      arma::colvec phi_post = arma::solve(trimatu(V_post_inv_chol), phi_post_tmp);

      arma::colvec rnormvec(K);
      for(int j = 0; j < K; ++j){
        rnormvec(j) = R::rnorm(0,1);
      }
      PHI.col(i) = phi_post + solve(trimatu(V_post_inv_chol), rnormvec);
    }else{

      //arma::mat V_post = ( V_prior_inv + X_new.t()*X_new ).i(); //unstable
      // compute inverse via the QR decomposition (similar to chol2inv in R)
      //arma::mat V_post_tmp = chol(V_prior_inv + X_new.t()*X_new, "upper");
      mat Q;
      mat R;
      qr(Q, R, V_post_inv_chol);
      mat S = inv(trimatu(R));
      mat V_post= S * S.t();

      arma::colvec phi_post = V_post * (V_prior_inv*PHI_prior.col(i) + X_new.t()*Y_new);

      PHI.col(i) = mvrnorm1(phi_post, V_post, 1e-06);
    }
  }
}


// draw_L samples the free off-diagonal elements in L, where L is an
// upper unitriangular matrix resulting from the L.inv*D*L.inv' factorization
// of the variance-covariance matrix as in Cogley & Sargent (2005)
// Rcpp::export --> for use within a sampler written in R

// [[Rcpp::export]]
arma::mat draw_L(arma::mat Ytilde, arma::vec& V_i, arma::mat& d) {

  int M = Ytilde.n_cols;
  mat L(M,M,fill::eye);

  int ind = 0;

  for(int i = 1; i < M; i++){

    vec normalizer = vectorise( pow(d.col(i), 0.5));
    mat Z = -1 * (Ytilde.cols(0, i-1).each_col() / normalizer);
    vec c = vectorise(Ytilde.col(i)) / normalizer;

    mat V_p_inv = diagmat(1/V_i.subvec(ind, ind+i-1));

    // Compute V_post
    //mat V_post = (V_p_inv + Z.t()*Z).i(); // unstable
    // compute inverse via the QR decomposition (similar to chol2inv in R)
    arma::mat V_post_tmp = chol(V_p_inv + Z.t()*Z, "upper"); //V_p_inv
    mat Q;
    mat R;
    qr(Q, R, V_post_tmp);
    mat S = inv(trimatu(R));
    mat V_post= S * S.t();

    colvec l_post = V_post * (Z.t()*c);

    L(span(0,i-1), span(i,i)) = mvrnorm1(l_post, V_post);

    ind = ind+i;

  }

  return(L);
}

void sample_L(arma::mat& L, arma::mat& Ytilde, const arma::vec& V_i, const arma::mat& d_sqrt) {

  int M = Ytilde.n_cols;
  //mat L(M,M,fill::eye);

  int ind = 0;

  for(int i = 1; i < M; i++){

    vec normalizer = vectorise( d_sqrt.col(i) );
    mat Z = -1 * (Ytilde.cols(0, i-1).each_col() / normalizer);
    vec c = vectorise(Ytilde.col(i)) / normalizer;

    mat V_p_inv = diagmat(1/V_i.subvec(ind, ind+i-1));

    // Compute V_post
    // mat V_post = (V_p_inv + Z.t()*Z).i(); // unstable
    // compute inverse via the QR decomposition (similar to chol2inv in R)
    arma::mat V_post_tmp = chol(V_p_inv + Z.t()*Z, "upper"); //V_p_inv
    mat Q;
    mat R;
    qr(Q, R, V_post_tmp);
    mat S = inv(trimatu(R));
    mat V_post= S * S.t();

    colvec l_post = V_post * (Z.t()*c);

    L(span(0,i-1), span(i,i)) = mvrnorm1(l_post, V_post);

    ind = ind+i;

  }

}

arma::colvec ddir_prep(const arma::colvec& x, const arma::vec& prep1, const arma::rowvec& prep2){

  //arma::rowvec logd = sum(prep1.each_col() % log(x), 0) + prep2;
  arma::rowvec logd(prep2.size());
  for(int j=0; j<prep2.size(); j++){
    logd(j) = sum(prep1(j) * log(x));
  }
  logd += prep2;

  return(logd.t());
}

void sample_V_i_DL(arma::vec& V_i, const arma::vec coefs, double& a ,
                   const arma::vec a_vec, const arma::vec prep1,
                   const arma::vec prep2, double& zeta, arma::vec& psi,
                   arma::vec& theta, arma::uvec ind, const bool hyper){ //, bool hyper

  double n = coefs.size();
  arma::vec coefs_abs = arma::abs(coefs);
  double zeta2 = zeta*zeta;
  arma::vec theta_prep(n);
  arma::uvec::iterator it;
  double j = 0;
    for(it = ind.begin(); it != ind.end(); ++it){ //int j = 0; j < n; j++
      psi(*it) = 1./do_rgig1(-0.5, 1, (coefs_abs(j) * coefs_abs(j)) /
        ( zeta2 * theta(*it) * theta(*it)));
      theta_prep(j) = do_rgig1(a-1., 2*coefs_abs(j), 1);
      j += 1;
    }

  double tmp4samplingzeta = arma::accu(coefs_abs / theta(ind));
  zeta = do_rgig1(n*(a-1.), 2*tmp4samplingzeta,1);
  theta(ind) = theta_prep / arma::accu(theta_prep);
  V_i(ind) = psi(ind) % theta(ind)%theta(ind) * zeta*zeta;

  if(hyper){
    int gridlength = a_vec.size();
    arma::vec logprobs = ddir_prep(theta(ind), prep1, prep2);
    for(int j=0; j<gridlength; j++){
      logprobs(j) += R::dgamma(zeta, n*a_vec(j), 1./0.5, true); // R::dgamma uses scale
    }

    arma::vec w = exp(logprobs - logprobs.max());
    arma::vec weights = w/sum(w);
    int k = weights.size();
    arma::ivec iv(k);
    R::rmultinom(1, weights.begin(), k, iv.begin());
    arma::uvec i = arma::find(iv == 1,1); // reports only the first value that meets the condition (by construction there is only one 1)
    a = a_vec(i(0));
  }

}

double logddirichlet(arma::vec x, double a){
  // symmetric dirichlet distribution
  int n = x.size();
  double d = lgamma(n*a) - n*lgamma(a) +
    sum((a - 1) * log(x));
  return(d);
}

void sample_V_i_R2D2(arma::vec& V_i, const arma::vec coefs, double& api,
                    const arma::vec api_vec, double& zeta, arma::vec& psi,
                    arma::vec& theta, double& xi, double& b,
                    const arma::vec b_vec, arma::uvec ind, const bool hyper){ //, bool hyper

  double n = coefs.size();
  arma::vec theta_prep(n);

  arma::uvec::iterator it;
  double j = 0;
  for(it = ind.begin(); it != ind.end(); ++it){ //int j = 0; j < n; j++
    psi(*it) = 1./do_rgig1(-0.5, 1, (coefs(j) * coefs(j)) /
      ( zeta * theta(*it)/2));
    theta_prep(j) = do_rgig1(api - .5, 2*coefs(j)*coefs(j)/psi(*it), 2*xi);
    j += 1;
  }

  double tmp4samplingzeta = arma::accu(square(coefs) / (theta(ind)%psi(ind)));
  zeta = do_rgig1(n*api - n/2, 2*tmp4samplingzeta, 2*xi);
  xi = R::rgamma(n*api + b, 1/(1+ zeta)); // uses scale
  theta(ind) = theta_prep / arma::accu(theta_prep);
  V_i(ind) = psi(ind) % theta(ind) * zeta / 2;

  if(hyper){
    int gridlength = b_vec.size();
    arma::vec logprobs(gridlength);
    for(int i=0; i<gridlength; ++i){

      logprobs(i) = R::dgamma(xi, b_vec(i), 1, true) +
        logddirichlet(theta(ind), api_vec(i));

    }
    arma::vec w_tmp = exp(logprobs - logprobs.max());
    arma::vec w = w_tmp/sum(w_tmp);

    int k = w.size();
    arma::ivec iv(k);
    R::rmultinom(1, w.begin(), k, iv.begin());
    arma::uvec ii = arma::find(iv == 1,1); // reports only the first value that meets the condition (by construction there is only one 1)
    b = b_vec(ii(0));
    api = api_vec(ii(0));
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
  double j = 0;
  for(it = ind.begin(); it != ind.end(); ++it){

    // Draw gammas
    gammas(*it) = R::rbinom(1,gst(*it));

    // Compute prior variances
    if(gammas(*it)==1){
      V_i(*it) = tau_1(*it)*tau_1(*it);
    }else{
      V_i(*it) = tau_0(*it)*tau_0(*it);
    }

//    if(type=="local"){
//      // Draw prior inclusion probabilities
//      p_i(*it) = R::rbeta(s_a + gammas(*it), s_b + 1 - gammas(*it));
//    }

  }

  if(hyper){
    int incl = accu(gammas(ind));
    double p = R::rbeta(s_a + incl, s_b + n - incl);
    //Rcout << "The value of n : " << n << "\n";
    p_i(ind).fill(p);
  }


}

void sample_V_i_SSVS(arma::vec& V_i, arma::vec& gammas, arma::vec& p_i,
                     const arma::vec& coeffs, const arma::vec& tau_0,
                     const arma::vec& tau_1, const double& s_a, const double& s_b,
                     const std::string type){

  int n = coeffs.size();
  // Compute conditional posterior inclusion probability
  // taking logs is much stabler!!!
  arma::vec u_i1 = -log(tau_1) - 0.5 * square((coeffs)/tau_1) + log(p_i);
  arma::vec u_i2 = -log(tau_0) - 0.5 * square((coeffs)/tau_0) + log(1 - p_i);
  arma::vec logdif = u_i2 - u_i1;
  arma::vec gst = 1/(1 + exp(logdif)); // == exp(u_i1)/(exp(u_i2) + exp(u_i1))

  for(int j=0; j<n; ++j){

    // Draw gammas
    gammas(j) = R::rbinom(1,gst(j));

    // Compute prior variances
    if(gammas(j)==1){
      V_i(j) = tau_1(j)*tau_1(j);
    }else{
      V_i(j) = tau_0(j)*tau_0(j);
    }

    if(type=="local"){
      // Draw prior inclusion probabilites
      p_i(j) = R::rbeta(s_a + gammas(j), s_b + 1 - gammas(j));
    }

  }

  if(type=="global"){
    int incl = accu(gammas);
    double p = R::rbeta(s_a + incl, s_b + n - incl);

    p_i.fill(p);
  }


}

void sample_V_i_HMP(double& lambda_1, double& lambda_2, arma::vec& V_i, const double& s1,
                    const double& r1, const double& s2, const double& r2,
                    const arma::vec& PHI_diff, const arma::vec& V_i_prep,
                    const int& n_ol, const int& n_cl, const arma::uvec& i_ol,
                    const arma::uvec& i_cl){

  lambda_1 = do_rgig1(s1 - n_ol/2, sum(square(PHI_diff(i_ol))/V_i_prep(i_ol)),
                     2*r1 );
  lambda_2 = do_rgig1(s2 - n_cl/2, sum(square(PHI_diff(i_cl))/V_i_prep(i_cl)), 2*r2 );

  V_i(i_ol) = lambda_1 * V_i_prep(i_ol);
  V_i(i_cl) = lambda_2 * V_i_prep(i_cl);

}

void sample_V_i_L_HMP(double& lambda_3, arma::vec& V_i_L, const double& s1,
                    const double& r1, const arma::vec& l){

  int n = l.size();
  lambda_3 = do_rgig1(s1 - n/2, sum(square(l)), 2*r1 );

  //for(int j=0; j<n; ++j){
  //  V_i_L(j) = lambda_3;
  //}
  V_i_L.fill(lambda_3);

}



void sample_PHI_SL(arma::mat& PHI, const arma::mat& PHI_prior, const arma::mat& Y,
                   const arma::mat& X, const arma::mat& L, const arma::mat& d_sqrt,
                   arma::mat& Gamma, const int& K, const int& M,
                   const double& nu_a, const double& nu_b) {

  for(int i = 0; i < M; i++){

    arma::mat PHI_0 = PHI;
    PHI_0.col(i).zeros();
    arma::vec normalizer = vectorise( d_sqrt.cols(i, M-1)  );
    arma::vec Y_new = vectorise( (Y - X * PHI_0) * L.cols(i, M-1) ) / normalizer;
    arma::mat X_new_tmp = arma::kron(L( span(i,i), span(i, M-1)).t(), X);
    arma::mat X_new = X_new_tmp.each_col() / normalizer;

    for(int j = 0; j<K; j++){
      //gammas
      double p_i = R::rbeta(0.5 + Gamma(j,i), 0.5 + 1 - Gamma(j,i));
      double v_i = 1./R::rgamma((nu_a + Gamma(j,i))/2, 1./((nu_b + PHI(j,i)*PHI(j,i))/2));

      arma::vec Y_new_star = Y_new - X_new * PHI.col(i) + X_new.col(j)*PHI(j,i);

      //vec pre1 = trans(X_new.col(j))*X_new.col(j); // pre1 is 1x1 dimensional
      double pre1 = accu(X_new.col(j)%X_new.col(j)) ;
      double V_post = 1./(1./v_i + pre1); // conversion to double via pre1(0)
      //vec pre2 = trans(X_new.col(j))*Y_new_star;
      double pre2 = accu(X_new.col(j)%Y_new_star) ;
      double phi_j_post = V_post*(PHI_prior(j,i)/v_i + pre2);

      double u_i1 = -0.5*log(v_i) - 0.5*(PHI_prior(j,i)*PHI_prior(j,i)/v_i) -
        (-0.5*log(V_post) - 0.5*(phi_j_post*phi_j_post/V_post)) + log(p_i);
      double u_i2 =  log(1 - p_i);

      double numericalnormalizer = u_i2;
      if(u_i1>u_i2){
        numericalnormalizer = u_i1;
      }

      double pp1 = exp(u_i1 - numericalnormalizer);
      double pp2 = exp(u_i2 - numericalnormalizer);
      double gst = pp1/(pp1 + pp2);
      //double logdiff = u_i2 - u_i1;
      //double gst = (1./(1+exp(logdiff)));

      Gamma(j,i) = R::rbinom(1,gst);

      //if(!(Gamma(j,i)==0 || Gamma(j,i)==1)){
      //  ::Rf_error("Gamma_ji: %f, u_i1: %f, u_i2: %f, v_i: %f, V_post: %f, phi_j_post: %f, pre1: %f", Gamma(j,i), u_i1, u_i2, v_i, V_post, phi_j_post, pre1);
      //}

      if(Gamma(j,i)==1){
        PHI(j,i) = R::rnorm(phi_j_post, sqrt(V_post));
      }else {
        PHI(j,i) = 0;
      }

    }

  }
}

void sample_L_SL(arma::mat& L, arma::mat& Ytilde, const arma::mat& d_sqrt,
              vec& omega, const double& nu_a, const double& nu_b) {

  int M = Ytilde.n_cols;
  //mat L(M,M,fill::eye);

  int ind = 0;

  for(int i = 1; i < M; i++){


    vec normalizer = vectorise( d_sqrt.col(i) );
    mat Z = -1 * (Ytilde.cols(0, i-1).each_col() / normalizer);
    vec c = vectorise(Ytilde.col(i)) / normalizer;

    for(int j=0; j<Z.n_cols; j++){

      double p_i = R::rbeta(0.5 + omega(ind), 0.5 + 1 - omega(ind));
      double v_i = 1./R::rgamma((nu_a + omega(ind))/2, 1./((nu_b + L(j,i)*L(j,i))/2));

      vec c_star = c - Z*L(span(0,i-1), span(i,i)) + Z.col(j)*L(j,i);

      vec pre1 = (trans(Z.col(j))*Z.col(j)); // pre1 is 1x1 dimensional
      double V_post = 1./(1./v_i + pre1(0)); // conversion to double via pre1(0)
      vec pre2 = (trans(Z.col(j))*c_star);
      double l_j_post = V_post*(pre2(0));

      double u_i1 = -0.5*log(v_i) - (-0.5*log(V_post) - 0.5*(l_j_post*l_j_post/V_post)) + log(p_i);
      double u_i2 =  log(1 - p_i);

      double logdiff = u_i2 - u_i1;
      double gst = (1./(1+exp(logdiff)));

      omega(ind) = R::rbinom(1,gst);

      if(omega(ind)==1){
        L(j,i) = R::rnorm(l_j_post, sqrt(V_post));
      }else {
        L(j,i) = 0;
      }
      ind = ind+1;
    }

  }

}

void sample_prior_mean(vec& m_i ,const vec& coefs, const vec& v_i, const double& mu0,
                       const double& b0){
  double n = m_i.n_elem;
  for(int i=0; i<n; i++){
    double b = 1./(1./b0 + 1./v_i(i));
    double mu = b * (mu0/b0 + coefs(i)/v_i(i));
    m_i(i) = R::rnorm(mu, sqrt(b));
    }
}

