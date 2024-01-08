#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// currently not in use
/*
arma::vec mvrnorm1(arma::vec& mu, arma::mat& Sigma, double tol = 1e-06){

  int p = mu.size();
  arma::rowvec rnormvec(p);
  rnormvec.imbue(R::norm_rand);
  // for(int j = 0; j < p; ++j){ loop is okayish. in terms of memory allocation as good as arma::vec.imbue(::norm_rand)/std::generate(vec.begin(),ve.end(), ::norm_rand), but slightly slower
  //   rnormvec(j) = R::rnorm(0,1);
  // }

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

    arma::vec u(K); u.imbue(R::norm_rand); // .imbue() is similar to .fill() but handles functors or lambda functions. R::norm_rand is the underlying function of wrappers like Rcpp::rnorm and R::rnorm.
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
    //rnormvec = Rcpp::rnorm(K); !!!inefficient: Rcpp::rnorm creates NumericVector, however rnormvec is arma::colvec. Hence there will be a copy
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

  const arma::mat Y_hat = Y - arma::trans(facload*fac); // conditioning on factors leads to independent equations
  const arma::mat normalizer = exp(-logvaridi/2.);

  for(int j = 0; j < M; j++){

    arma::mat X_norm = X.each_col() % normalizer.col(j);
    arma::vec y_norm = Y_hat.col(j) % normalizer.col(j);
    arma::vec PHI_j = PHI.unsafe_col(j);
    univariate_regression_update(PHI_j, PHI_prior.col(j), V_prior.col(j),
                                 y_norm, X_norm, false, huge);

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
    univariate_regression_update(PHI_j, PHI_prior.col(j), V_prior.col(j),
                                 Y_new, X_new, false, false);

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
    univariate_regression_update(U_j, prior_mean, V_i.subvec(ind, ind+i-1), c, Z, true, false);

    ind = ind+i;

  }

}


// void sample_U(arma::mat& U, arma::mat& Ytilde, const arma::vec& V_i, const arma::mat& d_sqrt) {
//
//   int M = Ytilde.n_cols;
//   //mat U(M,M,fill::eye);
//
//   int ind = 0;
//
//   for(int i = 1; i < M; i++){
//
//     vec normalizer = vectorise( d_sqrt.col(i) );
//     // create matrix with covariables (residuals of preceeding equations)
//     mat Z = -1 * (Ytilde.cols(0, i-1).each_col() / normalizer);
//     // create 'dependent' variable
//     vec c = vectorise(Ytilde.col(i)) / normalizer;
//     // prepare posterior precision
//     arma::mat Omegal_post = Z.t()*Z;
//     // add prior precision
//     Omegal_post.diag() += 1./V_i.subvec(ind, ind+i-1);
//     // prepare posterior mean
//     arma::vec l_post_prep = Z.t()*c;
//     // generate one posterior draw
//     arma::vec U_j = U.unsafe_col(i);
//     draw_post(U_j, Omegal_post, l_post_prep, true);
// ///    arma::vec l_tmp(size(l_post_prep));
// ///    draw_post(l_tmp, Omegal_post, l_post_prep, false);
// ///   U(span(0,i-1), span(i,i)) = l_tmp ;//mvrnorm1(l_post, V_post);
//
// /// old...
// //    mat V_p_inv = diagmat(1/V_i.subvec(ind, ind+i-1));
// //    // cholesky factor of posterior precision
// //    arma::mat V_post_tmp = chol(V_p_inv + Z.t()*Z, "upper"); //V_p_inv
// //    // Compute V_post
// //    // mat V_post = (V_p_inv + Z.t()*Z).i(); // unstable
// //    // compute inverse via the QR decomposition (similar to chol2inv in R)
// //    mat Q;
// //    mat R;
// //    qr(Q, R, V_post_tmp);
// //    mat S = inv(trimatu(R));
// //    mat V_post= S * S.t();
// //    colvec l_post = V_post * (Z.t()*c);
// //    U(span(0,i-1), span(i,i)) = mvrnorm1(l_post, V_post);
//
//     ind = ind+i;
//
//   }
//
// }

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

// void sample_V_i_NG(arma::vec& V_i, const arma::vec coefs, arma::vec& theta_tilde,
//                    double& zeta, double& a , const arma::vec a_vec,
//                    const double varrho0, const double varrho1, arma::uvec ind,
//                    const bool hyper, const double tol){
//
//   const int n = ind.size();
//   int gridlength = a_vec.size();
//   arma::vec logprobs(gridlength); logprobs.fill(0);
//
//   arma::uvec::iterator it;
//   for(it = ind.begin(); it != ind.end(); ++it){
//     theta_tilde(*it) = do_rgig(a-0.5, coefs(*it)*coefs(*it), a/zeta);
//     if(!theta_tilde.is_finite()){
//                     ::Rf_error("Non-finite theta_tilde: PHI(*): %e, a(j): %f, z(j): %e, *: %i, giglambda: %f, gigchi:%e, gigpsi:%e", coefs(*it), a, zeta, *it,a-0.5,coefs(*it)*coefs(*it),a/zeta);
//                   }
//     if(hyper){
//       for(int i=0; i<gridlength; ++i){
//           logprobs(i) += R::dgamma(theta_tilde(*it),a_vec(i), 2*zeta/a_vec(i), true); // scale!!!
//       }
//       }
//   }
//   if(tol>0){
//     double th_sum = accu(theta_tilde);
//     arma::vec theta = theta_tilde/th_sum;
//     arma::uvec th_uvec= find(theta<tol);
//     theta(th_uvec) = vec(th_uvec.size(), fill::value(tol));
//     theta_tilde = theta * th_sum;
//   }
//
//   zeta = 1./R::rgamma(varrho0 + a*n, 1./(varrho1 + 0.5*a*arma::accu(theta_tilde(ind))));
//
//   if(hyper){
//     arma::vec w_tmp = exp(logprobs - logprobs.max());
//     arma::vec w = w_tmp/sum(w_tmp);
//
//     if(!w.is_finite()){
//       ::Rf_error("sample_V_i_NG (NG_a = 'hyperprior'): non-finite weights");
//     }else if(sum(w)==0){
//       ::Rf_error("sample_V_i_NG (NG_a = 'hyperprior'): zero weights");
//     }
//
//     int k = w.size();
//     arma::ivec iv(k);
//     R::rmultinom(1, w.begin(), k, iv.begin());
//     arma::uvec i = arma::find(iv == 1,1); // reports only the first value that meets the condition (by construction there is only one 1)
//     a = a_vec(i(0));
//   }
//
//
//   V_i(ind) = theta_tilde(ind);
//
// }

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

// arma::colvec ddir_prep(const arma::colvec& x, const arma::vec& prep1, const arma::rowvec& prep2){
//
//   //arma::rowvec logd = sum(prep1.each_col() % log(x), 0) + prep2;
//   arma::rowvec logd(prep2.size());
//   for(unsigned int j=0; j<prep2.size(); j++){
//     logd(j) = sum(prep1(j) * log(x));
//   }
//   logd += prep2;
//
//   return(logd.t());
// }
//
// void sample_V_i_DL_deprecated(arma::vec& V_i, const arma::vec coefs, double& a ,
//                    const arma::vec a_vec, const arma::vec prep1,
//                    const arma::vec prep2, double& zeta, arma::vec& psi,
//                    arma::vec& theta, arma::uvec ind, const bool hyper,
//                    const int method, const double tol){ //, bool hyper
//
//   const int n = ind.size();
//   arma::vec lambda(n);
//   arma::vec coefs_abs = arma::abs(coefs);
//   double zeta2 = zeta*zeta;
//   arma::vec lambda0(theta.size());
//   arma::uvec::iterator it;
//     for(it = ind.begin(); it != ind.end(); ++it){ //int j = 0; j < n; j++
//       psi(*it) = 1./do_rgig(-0.5, 1, (coefs_abs(*it) * coefs_abs(*it)) /
//         ( zeta2 * theta(*it) * theta(*it)));
//       lambda0(*it) = do_rgig(a-1., 2*coefs_abs(*it), 1);
//     }
//
//   theta(ind) = lambda0(ind) / arma::accu(lambda0(ind));
//     if(tol>0){
//       uvec th_f = find(theta<tol);
//       theta(th_f) = vec(th_f.size(), fill::value(tol));
//     }
//   if(method == 1.){
//     double tmp4samplingzeta = arma::accu(coefs_abs(ind) / theta(ind));
//     zeta = do_rgig(n*(a-1.), 2*tmp4samplingzeta,1);
//     lambda = theta(ind) * zeta;
//   }else if(method == 2.){
//     zeta = arma::accu(lambda0(ind));
//     lambda = lambda0(ind);
//     if(tol>0){
//       lambda = theta(ind) * zeta;
//     }
//   }
//
//   V_i(ind) = psi(ind) % square(lambda);//theta(ind)%theta(ind) * zeta*zeta
//
//   if(hyper){
//     int gridlength = a_vec.size();
//     arma::vec logprobs = ddir_prep(theta(ind), prep1, prep2);
//     for(int j=0; j<gridlength; j++){
//       logprobs(j) += R::dgamma(zeta, n*a_vec(j), 1./0.5, true); // R::dgamma uses scale
//     }
//
//     arma::vec w = exp(logprobs - logprobs.max());
//     arma::vec weights = w/sum(w);
//     int k = weights.size();
//     arma::ivec iv(k);
//     R::rmultinom(1, weights.begin(), k, iv.begin());
//     arma::uvec i = arma::find(iv == 1,1); // reports only the first value that meets the condition (by construction there is only one 1)
//     a = a_vec(i(0));
//   }
//
// }
//
// double logddirichlet(arma::vec x, double a){
//   // symmetric dirichlet distribution
//   int n = x.size();
//   double d = lgamma(n*a) - n*lgamma(a) +
//     sum((a - 1) * log(x));
//   return(d);
// }
//
// void sample_V_i_R2D2(arma::vec& V_i, const arma::vec coefs, double& api,
//                     const arma::vec api_vec, double& zeta, arma::vec& psi,
//                     arma::vec& theta, double& xi, double& b,
//                     const arma::vec b_vec, arma::uvec ind, const bool hyper,
//                     const int method, const std::string kernel){ //, bool hyper
//
//   double n = coefs.size();
//   arma::vec theta_prep(n);
//
//   arma::uvec::iterator it;
//   double j = 0;
//   for(it = ind.begin(); it != ind.end(); ++it){ //int j = 0; j < n; j++
//     if(kernel=="laplace"){
//       psi(*it) = 1./do_rgig(-0.5, 1, (coefs(j) * coefs(j)) /
//         ( zeta * theta(*it)/2));
//       theta_prep(j) = do_rgig(api - .5, 2*coefs(j)*coefs(j)/psi(*it), 2*xi);
//     }else if(kernel == "normal"){
//       theta_prep(j) = do_rgig(api - .5, coefs(j)*coefs(j)/psi(*it), 2*xi);
//     }
//
//     j += 1;
//   }
//
//   if(method==1.){
//     double tmp4samplingzeta = arma::accu(square(coefs) / (theta(ind)%psi(ind)));
//     zeta = do_rgig(n*api - n/2, 2*tmp4samplingzeta, 2*xi);
//   }else if(method==2.){
//     zeta = arma::accu(theta_prep);
//   }
//
//   xi = R::rgamma(n*api + b, 1/(1+ zeta)); // uses scale
//   theta(ind) = theta_prep / arma::accu(theta_prep);
//   if(kernel=="laplace"){
//     V_i(ind) = psi(ind) % theta(ind) * zeta / 2;
//   }else if(kernel == "normal"){
//     V_i(ind) = theta(ind) * zeta;
//   }
//
//   if(hyper){
//     int gridlength = b_vec.size();
//     arma::vec logprobs(gridlength);
//     for(int i=0; i<gridlength; ++i){
//
//       logprobs(i) = R::dgamma(xi, b_vec(i), 1, true) +
//         logddirichlet(theta(ind), api_vec(i));
//
//     }
//     arma::vec w_tmp = exp(logprobs - logprobs.max());
//     arma::vec w = w_tmp/sum(w_tmp);
//
//     int k = w.size();
//     arma::ivec iv(k);
//     R::rmultinom(1, w.begin(), k, iv.begin());
//     arma::uvec ii = arma::find(iv == 1,1); // reports only the first value that meets the condition (by construction there is only one 1)
//     b = b_vec(ii(0));
//     api = api_vec(ii(0));
//   }
//
// }

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

// void sample_V_i_SSVS(arma::vec& V_i, arma::vec& gammas, arma::vec& p_i,
//                      const arma::vec& coeffs, const arma::vec& tau_0,
//                      const arma::vec& tau_1, const double& s_a, const double& s_b,
//                      const std::string type){
//
//   int n = coeffs.size();
//   // Compute conditional posterior inclusion probability
//   // taking logs is much stabler!!!
//   arma::vec u_i1 = -log(tau_1) - 0.5 * square((coeffs)/tau_1) + log(p_i);
//   arma::vec u_i2 = -log(tau_0) - 0.5 * square((coeffs)/tau_0) + log(1 - p_i);
//   arma::vec logdif = u_i2 - u_i1;
//   arma::vec gst = 1/(1 + exp(logdif)); // == exp(u_i1)/(exp(u_i2) + exp(u_i1))
//
//   for(int j=0; j<n; ++j){
//
//     // Draw gammas
//     gammas(j) = R::rbinom(1,gst(j));
//
//     // Compute prior variances
//     if(gammas(j)==1){
//       V_i(j) = tau_1(j)*tau_1(j);
//     }else{
//       V_i(j) = tau_0(j)*tau_0(j);
//     }
//
//     if(type=="local"){
//       // Draw prior inclusion probabilites
//       p_i(j) = R::rbeta(s_a + gammas(j), s_b + 1 - gammas(j));
//     }
//
//   }
//
//   if(type=="global"){
//     int incl = accu(gammas);
//     double p = R::rbeta(s_a + incl, s_b + n - incl);
//
//     p_i.fill(p);
//   }
//
//
// }

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

  //for(int j=0; j<n; ++j){
  //  V_i_U(j) = lambda_3;
  //}
  V_i_U.fill(lambda_3);

}



// void sample_PHI_SL(arma::mat& PHI, const arma::mat& PHI_prior, const arma::mat& Y,
//                    const arma::mat& X, const arma::mat& U, const arma::mat& d_sqrt,
//                    arma::mat& Gamma, const int& K, const int& M,
//                    const double& nu_a, const double& nu_b) {
//
//   for(int i = 0; i < M; i++){
//
//     arma::mat PHI_0 = PHI;
//     PHI_0.col(i).zeros();
//     arma::vec normalizer = vectorise( d_sqrt.cols(i, M-1)  );
//     arma::vec Y_new = vectorise( (Y - X * PHI_0) * U.cols(i, M-1) ) / normalizer;
//     arma::mat X_new_tmp = arma::kron(U( span(i,i), span(i, M-1)).t(), X);
//     arma::mat X_new = X_new_tmp.each_col() / normalizer;
//
//     for(int j = 0; j<K; j++){
//       //gammas
//       double p_i = R::rbeta(0.5 + Gamma(j,i), 0.5 + 1 - Gamma(j,i));
//       double v_i = 1./R::rgamma((nu_a + Gamma(j,i))/2, 1./((nu_b + PHI(j,i)*PHI(j,i))/2));
//
//       arma::vec Y_new_star = Y_new - X_new * PHI.col(i) + X_new.col(j)*PHI(j,i);
//
//       //vec pre1 = trans(X_new.col(j))*X_new.col(j); // pre1 is 1x1 dimensional
//       double pre1 = accu(X_new.col(j)%X_new.col(j)) ;
//       double V_post = 1./(1./v_i + pre1); // conversion to double via pre1(0)
//       //vec pre2 = trans(X_new.col(j))*Y_new_star;
//       double pre2 = accu(X_new.col(j)%Y_new_star) ;
//       double phi_j_post = V_post*(PHI_prior(j,i)/v_i + pre2);
//
//       double u_i1 = -0.5*log(v_i) - 0.5*(PHI_prior(j,i)*PHI_prior(j,i)/v_i) -
//         (-0.5*log(V_post) - 0.5*(phi_j_post*phi_j_post/V_post)) + log(p_i);
//       double u_i2 =  log(1 - p_i);
//
//       double numericalnormalizer = u_i2;
//       if(u_i1>u_i2){
//         numericalnormalizer = u_i1;
//       }
//
//       double pp1 = exp(u_i1 - numericalnormalizer);
//       double pp2 = exp(u_i2 - numericalnormalizer);
//       double gst = pp1/(pp1 + pp2);
//       //double logdiff = u_i2 - u_i1;
//       //double gst = (1./(1+exp(logdiff)));
//
//       Gamma(j,i) = R::rbinom(1,gst);
//
//       //if(!(Gamma(j,i)==0 || Gamma(j,i)==1)){
//       //  ::Rf_error("Gamma_ji: %f, u_i1: %f, u_i2: %f, v_i: %f, V_post: %f, phi_j_post: %f, pre1: %f", Gamma(j,i), u_i1, u_i2, v_i, V_post, phi_j_post, pre1);
//       //}
//
//       if(Gamma(j,i)==1){
//         PHI(j,i) = R::rnorm(phi_j_post, sqrt(V_post));
//       }else {
//         PHI(j,i) = 0;
//       }
//
//     }
//
//   }
// }
//
// void sample_U_SL(arma::mat& U, arma::mat& Ytilde, const arma::mat& d_sqrt,
//               vec& omega, const double& nu_a, const double& nu_b) {
//
//   int M = Ytilde.n_cols;
//   //mat U(M,M,fill::eye);
//
//   int ind = 0;
//
//   for(int i = 1; i < M; i++){
//
//
//     vec normalizer = vectorise( d_sqrt.col(i) );
//     mat Z = -1 * (Ytilde.cols(0, i-1).each_col() / normalizer);
//     vec c = vectorise(Ytilde.col(i)) / normalizer;
//
//     for(unsigned int j=0; j<Z.n_cols; j++){
//
//       double p_i = R::rbeta(0.5 + omega(ind), 0.5 + 1 - omega(ind));
//       double v_i = 1./R::rgamma((nu_a + omega(ind))/2, 1./((nu_b + U(j,i)*U(j,i))/2));
//
//       vec c_star = c - Z*U(span(0,i-1), span(i,i)) + Z.col(j)*U(j,i);
//
//       vec pre1 = (trans(Z.col(j))*Z.col(j)); // pre1 is 1x1 dimensional
//       double V_post = 1./(1./v_i + pre1(0)); // conversion to double via pre1(0)
//       vec pre2 = (trans(Z.col(j))*c_star);
//       double l_j_post = V_post*(pre2(0));
//
//       double u_i1 = -0.5*log(v_i) - (-0.5*log(V_post) - 0.5*(l_j_post*l_j_post/V_post)) + log(p_i);
//       double u_i2 =  log(1 - p_i);
//
//       double logdiff = u_i2 - u_i1;
//       double gst = (1./(1+exp(logdiff)));
//
//       omega(ind) = R::rbinom(1,gst);
//
//       if(omega(ind)==1){
//         U(j,i) = R::rnorm(l_j_post, sqrt(V_post));
//       }else {
//         U(j,i) = 0;
//       }
//       ind = ind+1;
//     }
//
//   }
//
// }

