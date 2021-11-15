// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

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
                const arma::mat& X, const arma::mat& L, const arma::mat& d,
                const arma::vec& V_i, const int& K, const int& M) {

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

void sample_L(arma::mat& L, arma::mat& Ytilde, arma::vec& V_i, arma::mat& d) {

  // Import MASS::mvrnorm function
  Environment MASS = Environment::namespace_env("MASS");
  Function Mrmvnorm = MASS["mvrnorm"];

  int M = Ytilde.n_cols;
  //mat L(M,M,fill::eye);

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

}
