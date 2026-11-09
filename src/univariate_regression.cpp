#include <RcppArmadillo.h>

// linear regression, gaussian error term with known variance of 1
// multivariate normal prior on coefs
// prior covariance of coefficients is assumed to be diagonal
// quite robust algorithm, relying on forward and backward substitution of 
//  cholesky factors instead of naively inverting the full conditional posterior
//  covariance matrix
arma::vec sample_coef_blm_known1var_trad(const arma::vec& y, const arma::mat& X,
                                         const arma::vec& coef_prior_mean,
                                         const arma::vec& coef_prior_variance){
  arma::vec coef_draw;
  const int p = X.n_cols;
  Rcpp::NumericVector rcpp_normrandvec = Rcpp::rnorm(p, 0.0, 1.0);
  arma::vec normrandvec(rcpp_normrandvec.begin(), rcpp_normrandvec.length(), false);
  arma::vec Xy_plus = X.t()*y + coef_prior_mean/coef_prior_variance;
  try {
    // see, Appendix A.1.1.2 in https://doi.org/10.1016/j.jeconom.2018.11.006
    // this alternative computation of the full conditional posterior precision
    // matrix is numerically more robust if some elements of coef_prior_variance
    // are close to zero (shrinkage prior)
    arma::vec scales = arma::sqrt(coef_prior_variance);
    arma::mat XX_scaled_plus1 = X.t()*X;
    XX_scaled_plus1.each_col() %= scales;
    XX_scaled_plus1.each_row() %= scales.as_row();
    XX_scaled_plus1.diag() += 1.0;
    arma::mat XX_scaled_plus1_chol = arma::chol(XX_scaled_plus1);
    arma::vec rhs1 = Xy_plus%scales;
    arma::vec coef_post_tmp = arma::solve(arma::trimatl(XX_scaled_plus1_chol.t()),
                                          rhs1, arma::solve_opts::no_approx);
    arma::mat rhs2 = arma::join_rows(coef_post_tmp, normrandvec);
    arma::mat coef_post_mat = arma::solve(arma::trimatu(XX_scaled_plus1_chol),
                                          rhs2, arma::solve_opts::no_approx);
    coef_draw = coef_post_mat.col(0)%scales + coef_post_mat.col(1)%scales;
  } catch (...) {
    // fall back to traditional approach
    arma::mat coef_post_precision = X.t()*X;
    coef_post_precision.diag() += 1.0/coef_prior_variance;
    arma::mat coef_post_precision_chol = arma::chol(coef_post_precision);
    arma::vec coef_post_tmp = arma::solve(arma::trimatl(coef_post_precision_chol.t()),
                                          Xy_plus, arma::solve_opts::no_approx);
    arma::vec rhs = coef_post_tmp + normrandvec;
    coef_draw = arma::solve(arma::trimatu(coef_post_precision_chol),
                            rhs, arma::solve_opts::no_approx);
  }
  return coef_draw;
}

// linear regression, gaussian error term with known variance of 1
// multivariate normal prior on coefs
// prior covariance of coefficients is assumed to be diagonal
// this algo exploits qr decomposition of augmented system
arma::vec sample_coef_blm_known1var_qr(const arma::vec& y, const arma::mat& X,
                                    const arma::vec& coef_prior_mean,
                                    const arma::vec& coef_prior_variance){
  const int p = X.n_cols;
  
  arma::vec coef_prior_variance_sqrt = arma::sqrt(coef_prior_variance);
  arma::mat augmentedX = arma::join_cols(X, arma::diagmat(1.0/coef_prior_variance_sqrt));
  arma::vec augmentedy = arma::join_cols(y, coef_prior_mean/coef_prior_variance_sqrt);
  
  arma::mat Q;
  arma::mat R;
  arma::uvec P;
  const bool s0 = arma::qr_econ(Q, R, augmentedX);
  if(!s0) {
    if(!arma::qr(Q, R, P, augmentedX, "vector")) {
      throw std::runtime_error("Pivoted QR failed.");
    }
  }
  Rcpp::NumericVector rcpp_rhs = Rcpp::rnorm(p, 0.0, 1.0);
  arma::vec rhs(rcpp_rhs.begin(), rcpp_rhs.length(), false);
  rhs += Q.t()*augmentedy;
  arma::vec coef_draw = arma::solve(arma::trimatu(R), rhs);
  if(!s0) {
    coef_draw(P) = coef_draw;
  }
    
  return coef_draw;
}

// linear regression, gaussian error term with known variance of 1
// multivariate normal prior on coefs
// prior covariance of coefficients is assumed to be diagonal
// variant of the algorithm of doi: 10.1093/biomet/asw042 that can handle a
// non-zero prior mean
arma::vec sample_coef_blm_known1var_large_p_small_n(const arma::vec& y, const arma::mat& X,
                                                    const arma::vec& coef_prior_mean,
                                                    const arma::vec& coef_prior_variance){
  arma::uword p = X.n_cols;
  arma::uword n = X.n_rows;
  arma::vec u(p); u.imbue(R::norm_rand);
  u %= arma::sqrt(coef_prior_variance);
  u += coef_prior_mean; 
  arma::vec delta(n); delta.imbue(R::norm_rand);
  arma::mat X_scaled_trans = X.t(); 
  X_scaled_trans.each_col() %= coef_prior_variance;
  arma::mat tosolve = X*X_scaled_trans;
  tosolve.diag() += 1.0;
  arma::vec rhs = y- X*u + delta;
  arma::vec w = arma::solve(tosolve, rhs, arma::solve_opts::likely_sympd);
  
  return u + X_scaled_trans*w;
}

// wrapper
arma::vec sample_coef_blm_known1var(const arma::vec& y, const arma::mat& X,
                                    const arma::vec& coef_prior_mean,
                                    const arma::vec& coef_prior_variance,
                                    const bool huge) {
  arma::vec coef_draw;
  if (huge) {
    coef_draw = sample_coef_blm_known1var_large_p_small_n(y, X, coef_prior_mean,
                                                          coef_prior_variance);
  } else {
    try {
      coef_draw = sample_coef_blm_known1var_trad(y, X, coef_prior_mean,
                                                 coef_prior_variance);
    } catch (...) {
      coef_draw = sample_coef_blm_known1var_qr(y, X, coef_prior_mean,
                                               coef_prior_variance);
    }
  }
  return coef_draw;
}