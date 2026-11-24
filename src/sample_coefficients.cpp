#include <RcppArmadillo.h>
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

//Kastner and Huber 2019 JoF equation per equation
void sample_PHI_factor(arma::mat& PHI, const arma::mat& PHI_prior,
                       const arma::mat& Y, const arma::mat& X,
                       const arma::mat& logvaridi, const arma::mat& V_prior,
                       const arma::mat& facload, const arma::mat& fac,
                       const bool& huge, const double floor,
                       const bool stable, const double stable_lim){

  const int M = Y.n_cols;
  const int r = fac.n_rows;

  const arma::mat Y_hat = (r > 0) ? Y - arma::trans(facload*fac) : Y; // conditioning on factors leads to independent equations
  
  const arma::mat normalizer = exp(-logvaridi/2.);
  
  const arma::mat PHI_old = PHI;
  for(int j = 0; j < M; j++){

    arma::mat X_norm = X.each_col() % normalizer.col(j);
    arma::vec y_norm = Y_hat.col(j) % normalizer.col(j);
    arma::vec PHI_j = PHI.unsafe_col(j);
    arma::vec PHI_prior_j = PHI_prior.unsafe_col(j);
    arma::vec V_prior_j = V_prior.unsafe_col(j);
    
    try {
      PHI_j = sample_coef_blm_known1var(y_norm, X_norm, PHI_prior_j, V_prior_j,
                                        huge);
    } catch(const std::exception& e) {
      Rcpp::stop("sample_coef_blm_known1var failed in %ith eqation: %s", j+1, e.what());
    } catch (...) {
      Rprintf("\nsample_coef_blm_known1var failed in %ith eqation. Rethrowing exception:", j+1);
      throw;
    }
    if(floor>0) {
      for (arma::uword i=0; i<PHI_j.n_elem; ++i) {
        abs_floor2(PHI_j(i), PHI_prior_j(i), floor);
      }
    }
  }
  if(stable) {
    stable_PHI(PHI, PHI_old, M, stable_lim);
  }
}

// sample_PHI samples VAR coefficients using the corrected triangular
// algorithm as in Carrier, Chan, Clark & Marcellino (2021)
void sample_PHI_chol(arma::mat& PHI, const arma::mat& PHI_prior, const arma::mat& Y,
                     const arma::mat& X, const arma::mat& U, const arma::mat& logvar,
                     const arma::mat& V_prior, const int M, const double floor,
                     const bool stable, const double stable_lim) {
  const arma::mat PHI_old = PHI;
  for(int j = 0; j < M; j++){
    
    arma::mat PHI_0 = PHI;
    PHI_0.col(j).zeros();
    arma::mat normalizer_mat = arma::exp(0.5*logvar.cols(j, M-1));//d_sqrt.cols(j, M-1);
    arma::vec normalizer(normalizer_mat.begin(), normalizer_mat.size(), false);
    arma::mat Y_new_mat = (Y - X * PHI_0) * U.cols(j, M-1);
    arma::vec Y_new(Y_new_mat.begin(), Y_new_mat.size(), false);
    Y_new /= normalizer;
    arma::mat X_new_tmp = arma::kron(U( span(j,j), span(j, M-1)).t(), X);
    arma::mat X_new = X_new_tmp.each_col() / normalizer;
    
    arma::vec PHI_j = PHI.unsafe_col(j);
    arma::vec PHI_prior_j = PHI_prior.unsafe_col(j);
    arma::vec V_prior_j = V_prior.unsafe_col(j);
    
    try {
      PHI_j = sample_coef_blm_known1var(Y_new, X_new, PHI_prior_j, V_prior_j,
                                        false);
    } catch(const std::exception& e) {
      Rcpp::stop("sample_coef_blm_known1var failed in %ith eqation: %s", j+1, e.what());
    } catch (...) {
      Rprintf("\nsample_coef_blm_known1var failed in %ith eqation. Rethrowing exception:", j+1);
      throw;
    }
    
    if(floor>0) {
      for (arma::uword i=0; i<PHI_j.n_elem; ++i) {
        abs_floor2(PHI_j(i), PHI_prior_j(i), floor);
      }
    }
  }
  if(stable) {
    stable_PHI(PHI, PHI_old, M, stable_lim);
  }
}

void sample_PHI(arma::mat& PHI, const arma::mat& PHI_prior,
                const arma::mat& Y, const arma::mat& X,
                const arma::mat& V_prior, const double floor, const bool stable, 
                const double stable_lim, const std::string sigma_type,
                const arma::mat& logvar, const arma::mat& facload, 
                const arma::mat& fac, const bool& huge, 
                const arma::mat& U, const int M) {
  if(sigma_type == "factor") {
    const arma::mat logvaridi = logvar.cols(0,M-1);
    sample_PHI_factor(PHI, PHI_prior, Y, X, logvaridi, V_prior,
                      facload, fac, huge, floor, stable, stable_lim);
  } else if (sigma_type == "cholesky") {
    sample_PHI_chol(PHI, PHI_prior, Y, X, U, logvar, V_prior, M, floor,
                    stable, stable_lim);
  } else {
    Rcpp::stop("sample_PHI() unaware of 'sigma_type' " + sigma_type);
  }
  
}

// [[Rcpp::export]]
arma::mat sample_PHI_cholesky(const arma::mat& PHI, const arma::mat& PHI_prior,
                              const arma::mat& Y, const arma::mat& X,
                              const arma::mat& U, const arma::mat& logvar,
                              const arma::mat& V_prior,
                              const double floor, const bool stable,
                              const double stable_lim){
  const int M = Y.n_cols;
  arma::mat PHI_ret = PHI;
  sample_PHI_chol(PHI_ret, PHI_prior, Y, X, U, logvar, V_prior, M, floor, stable, stable_lim);
  return(PHI_ret);
}


void sample_U(arma::mat& U, const arma::mat& Ytilde, const arma::vec& V_i, const arma::mat& logvar,
              const double floor) {
  int M = Ytilde.n_cols;
  int ind = 0;
  for(int i = 1; i < M; i++){
    arma::vec normalizer = arma::exp(0.5*logvar.col(i));
    // create matrix with covariables (residuals of preceeding equations)
    arma::mat Z = Ytilde.cols(0, i-1);
    Z = -1.0 * (Z.each_col() /normalizer);
    arma::vec prior_mean(i, arma::fill::zeros);
    // create 'dependent' variable
    arma::vec c = Ytilde.col(i) / normalizer;
    // update free off-diagonal elements in column i
    arma::vec U_j = U.unsafe_col(i);
    try {
      arma::vec tmp = sample_coef_blm_known1var(c, Z, prior_mean, V_i.subvec(ind, ind+i-1), false);
      // V_i.subvec not passed by reference in that case, because subvec is a subview
      // @ huge=false: could true be an option if M>>T?
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
