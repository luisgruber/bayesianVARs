#include <RcppArmadillo.h>

// linear regression, gaussian error term with known conditional variance of 1
// multivariate normal prior on coefs
// prior covariance of coefficients is assumed to be diagonal
// quite robust algorithm, relying on forward and backward substitution of 
//  cholesky factors instead of naively inverting the full conditional posterior
//  covariance matrix
arma::vec sample_coef_blm_known1var(const arma::vec& Xy, const arma::mat& XX,
                                    const arma::vec& coef_prior_mean,
                                    const arma::vec& coef_prior_precision){
  arma::vec coef_draw;
  const int p = XX.n_cols;
  Rcpp::NumericVector rcpp_normrandvec = Rcpp::rnorm(p, 0.0, 1.0);
  arma::vec normrandvec(rcpp_normrandvec.begin(), rcpp_normrandvec.length(), false);
  arma::vec Xy_plus = Xy+coef_prior_mean%coef_prior_precision;
  try {
    // see, Appendix A.1.1.2 in https://doi.org/10.1016/j.jeconom.2018.11.006
    // this alternative computation of the full conditional posterior precision
    // matrix is numerically more robust if some elements of coef_prior_precision
    // are large (shrinkage prior)
    arma::vec scales = arma::sqrt(coef_prior_precision);
    arma::vec add(p, arma::fill::ones);
    arma::uvec idx = arma::find( scales == 0);
    scales.elem(idx).ones();
    add.elem(idx).zeros();
    // const arma::vec coef_prior_precision_sqrt = arma::sqrt(coef_prior_precision);
    arma::mat XX_scaled_plus1 = XX;
    XX_scaled_plus1.each_col() /= scales;
    XX_scaled_plus1.each_row() /= scales.as_row();
    XX_scaled_plus1.diag() += add;//1.0;
    arma::mat XX_scaled_plus1_chol = arma::chol(XX_scaled_plus1);
    arma::vec coef_post_mean = Xy_plus/scales;
    coef_post_mean = arma::solve(arma::trimatl(XX_scaled_plus1_chol.t()),
                                 coef_post_mean, arma::solve_opts::no_approx);
    coef_post_mean = arma::solve(arma::trimatu(XX_scaled_plus1_chol),
                                 coef_post_mean, arma::solve_opts::no_approx);
    coef_post_mean /= scales;
    coef_draw = arma::solve(arma::trimatu(XX_scaled_plus1_chol), normrandvec,
                            arma::solve_opts::no_approx);
    coef_draw /= scales;
    coef_draw += coef_post_mean;
  } catch (...) {
    // fall back to traditional approach
    arma::mat coef_post_precision = XX;
    coef_post_precision.diag() += coef_prior_precision;
    arma::mat coef_post_precision_chol = arma::chol(coef_post_precision);
    arma::vec coef_post_mean = arma::solve(arma::trimatl(coef_post_precision_chol.t()),
                                           Xy_plus, arma::solve_opts::no_approx);
    coef_post_mean = arma::solve(arma::trimatu(coef_post_precision_chol),
                                 coef_post_mean, arma::solve_opts::no_approx);
    coef_draw = coef_post_mean + arma::solve(arma::trimatu(coef_post_precision_chol),
                                             normrandvec, arma::solve_opts::no_approx);
  }
  return coef_draw;
}

// linear regression, gaussian error term with known conditional variance of 1
// multivariate normal prior on coefs
// prior covariance of coefficients is assumed to be diagonal
// variant of the algorithm of doi: 10.1093/biomet/asw042 that can handle a
// non-zero prior mean
arma::vec sample_coef_blm_known1var_large_p_small_n(const arma::vec& y, const arma::mat& X,
                                                    const arma::vec& coef_prior_mean,
                                                    const arma::vec& coef_prior_precision){
  arma::uword p = X.n_cols;
  arma::uword n = X.n_rows;
  arma::vec u(p); u.imbue(R::norm_rand);
  u /= arma::sqrt(coef_prior_precision);
  u += coef_prior_mean;
  arma::vec delta(n); delta.imbue(R::norm_rand);
  arma::vec v = X*u + delta;
  arma::mat X_scaled_trans = X.t(); //X.each_row()/coef_prior_precision.as_row();
  X_scaled_trans.each_col() /= coef_prior_precision;
  arma::mat tosolve = X*X_scaled_trans;//X_scaled*X.t();
  tosolve.diag() += 1.0;
  arma::vec w = arma::solve(tosolve, y-v, arma::solve_opts::likely_sympd);
  
  return u + X_scaled_trans*w;
}