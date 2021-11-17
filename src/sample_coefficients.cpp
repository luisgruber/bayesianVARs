// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Import rgig from R package GIGrvg
double do_rgig1(double lambda, double chi, double psi) {

  double res;
  // circumvent GIGrvg in these cases
  if ((chi < (11 * DOUBLE_EPS)) & (lambda != 0)) {
    /* special cases which are basically Gamma and Inverse Gamma distribution */
    if (lambda > 0.0) {
      res = R::rgamma(lambda, 2.0/psi);
    }
    else {
      res = 1.0/R::rgamma(-lambda, 2.0/chi); // fixed
    }
  } else if ((psi < (11 * DOUBLE_EPS)) & (lambda != 0)) {
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
arma::mat draw_PHI(arma::mat PHI, arma::mat PHI_prior, arma::mat Y, arma::mat X,
                   arma::mat L, arma::mat d, arma::vec V_i, int K, int M) {

  // Import MASS::mvrnorm function
  Environment MASS = Environment::namespace_env("MASS");
  Function Mrmvnorm = MASS["mvrnorm"];
  //Environment base = Environment("package:base");
  //Function Rchol = base["chol"];


  arma::mat V_prior = reshape(V_i, K, M);

  for(int i = 0; i < M; i++){

    arma::mat V_p_inv = diagmat(1/V_prior.col(i));

    arma::mat PHI_0 = PHI;
    PHI_0.col(i).zeros();
    arma::vec normalizer = vectorise( pow(d.cols(i, M-1), 0.5)  );
    arma::vec Y_new = vectorise( (Y - X * PHI_0) * L.cols(i, M-1) ) / normalizer;
    arma::mat X_new_tmp = arma::kron(L( span(i,i), span(i, M-1)).t(), X);
    arma::mat X_new = X_new_tmp.each_col() / normalizer;

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

    mat V_post_chol(K,K);
    bool chol_success = chol(V_post_chol, V_post,"lower");

    // Fall back on MASS::mvnorm  if armadillo fails
    if(chol_success == true){

      arma::colvec randnvec(K);

      for(int j = 0; j < K; j++) {

        randnvec(j) = R::rnorm(0, 1);

      }

      PHI.col(i) = phi_post + V_post_chol * randnvec;

    }else if(chol_success == false){

      NumericVector tmp = Mrmvnorm(1, Named("mu")=phi_post, Named("Sigma")=V_post); //, Named("tol")=100
      PHI.col(i) = as<arma::colvec>(tmp);//arma::mvnrnd(phi_post, V_post)


    }

  }

  return(PHI);
}


// draw_PHI samples VAR coefficients using the corrected triangular
// algorithm as in Carrier, Chan, Clark & Marcellino (2021)
// for use within smapler written in Rcpp
void sample_PHI(arma::mat& PHI, const arma::mat& PHI_prior, const arma::mat& Y,
                const arma::mat& X, const arma::mat& L, const arma::mat& d_sqrt,
                const arma::mat& V_prior, const int& K, const int& M) {

  // Import MASS::mvrnorm function
  Environment MASS = Environment::namespace_env("MASS");
  Function Mrmvnorm = MASS["mvrnorm"];
  //Environment base = Environment("package:base");
  //Function Rchol = base["chol"];

  for(int i = 0; i < M; i++){

    arma::mat V_p_inv = diagmat(1/V_prior.col(i));

    arma::mat PHI_0 = PHI;
    PHI_0.col(i).zeros();
    arma::vec normalizer = vectorise( d_sqrt.cols(i, M-1)  );
    arma::vec Y_new = vectorise( (Y - X * PHI_0) * L.cols(i, M-1) ) / normalizer;
    arma::mat X_new_tmp = arma::kron(L( span(i,i), span(i, M-1)).t(), X);
    arma::mat X_new = X_new_tmp.each_col() / normalizer;

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

    mat V_post_chol(K,K);
    bool chol_success = chol(V_post_chol, V_post,"lower");

    // Fall back on MASS::mvnorm  if armadillo fails
    if(chol_success == true){

      arma::colvec randnvec(K);

      for(int j = 0; j < K; j++) {

        randnvec(j) = R::rnorm(0, 1);

      }

      PHI.col(i) = phi_post + V_post_chol * randnvec;

    }else if(chol_success == false){

      NumericVector tmp = Mrmvnorm(1, Named("mu")=phi_post, Named("Sigma")=V_post); //, Named("tol")=100
      PHI.col(i) = as<arma::colvec>(tmp);//arma::mvnrnd(phi_post, V_post)


    }

  }
}


// draw_L samples the free off-diagonal elements in L, where L is an
// upper unitriangular matrix resulting from the L.inv*D*L.inv' factorization
// of the variance-covariance matrix as in Cogley & Sargent (2005)
// Rcpp::export --> for use within a sampler written in R

// [[Rcpp::export]]
arma::mat draw_L(arma::mat Ytilde, arma::vec V_i, arma::mat d) {

  // Import MASS::mvrnorm function
  Environment MASS = Environment::namespace_env("MASS");
  Function Mrmvnorm = MASS["mvrnorm"];

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

    mat V_post_chol(i,i);
    bool chol_success = chol(V_post_chol, V_post,"lower");

    // Fall back on MASS::mvnorm  if armadillo fails
    if(chol_success == false){

      NumericVector tmp = Mrmvnorm(1, l_post, V_post);
      L(span(0,i-1), span(i,i)) = as<arma::colvec>(tmp);//arma::mvnrnd(l_post, V_post)

    }else if(chol_success == true){

      arma::colvec randnvec(i);

      for(int j = 0; j < i; j++) {

        randnvec(j) = R::rnorm(0, 1);

      }

      L(span(0,i-1), span(i,i)) = l_post + V_post_chol * randnvec;

    }

    ind = ind+i;

  }

  return(L);
}

void sample_L(arma::mat& L, arma::mat& Ytilde, arma::vec& V_i, arma::mat& d_sqrt) {

  // Import MASS::mvrnorm function
  Environment MASS = Environment::namespace_env("MASS");
  Function Mrmvnorm = MASS["mvrnorm"];

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

    mat V_post_chol(i,i);
    bool chol_success = chol(V_post_chol, V_post,"lower");

    // Fall back on MASS::mvnorm  if armadillo fails
    if(chol_success == false){

      NumericVector tmp = Mrmvnorm(1, l_post, V_post);
      L(span(0,i-1), span(i,i)) = as<arma::colvec>(tmp);//arma::mvnrnd(l_post, V_post)

    }else if(chol_success == true){

      arma::colvec randnvec(i);

      for(int j = 0; j < i; j++) {

        randnvec(j) = R::rnorm(0, 1);

      }

      L(span(0,i-1), span(i,i)) = l_post + V_post_chol * randnvec;

    }

    ind = ind+i;

  }

}

void sample_V_i_DL(arma::vec& V_i, const arma::vec& coefs, const double& a ,
                   double& zeta, arma::vec& psi, arma::vec&theta){

  double n = coefs.size();
  double tmp4samplingzeta = 0;
  arma::vec theta_prep(n);

    for(int j = 0; j < n; j++){
      psi(j) = 1./do_rgig1(-0.5, 1, (coefs(j) * coefs(j)) /
        ( zeta * zeta * theta(j) * theta(j)));
      tmp4samplingzeta += fabs(coefs(j))/theta(j) ;
      theta_prep(j) = do_rgig1(a-1., 2*fabs(coefs(j)), 1);
    }

  zeta = do_rgig1(n*(a-1.), 2*tmp4samplingzeta,1);
  theta = theta_prep / arma::accu(theta_prep);
  V_i = psi % theta%theta * zeta*zeta  ;
}

arma::colvec ddir_prep(arma::colvec x, arma::mat prep1, arma::rowvec prep2){

  arma::rowvec logd = sum(prep1.each_col() % log(x), 0) + prep2;

  return(logd.t());
}

void sample_DL_hyper(double& a, const arma::vec& theta, const arma::mat& prep1,
                     const arma::rowvec& prep2, const double& zeta,
                     arma::vec& a_vec){
  const int n = theta.size();
  arma::vec logprobs = ddir_prep(theta, prep1, prep2);
  for(int j=0; j<1000; j++){
    logprobs(j) += R::dgamma(zeta, n*a_vec(j), 1/0.5, true); // R::dgamma uses scale
  }

  arma::vec w = exp(logprobs - logprobs.max());
  arma::vec weights = w/sum(w);
  arma::ivec iv(1000);
  R::rmultinom(1, weights.begin(), 1000, iv.begin());
  arma::uvec i = arma::find(iv == 1,1); // reports only the first value that meets the condition (by construction there is only one 1)
  a = a_vec(i(0));
}
