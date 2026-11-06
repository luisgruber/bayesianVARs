#ifndef UNIVARIATE_REGRESSION_H
#define UNIVARIATE_REGRESSION_H

arma::vec sample_coef_blm_known1var(const arma::vec& Xy, const arma::mat& XX,
                                    const arma::vec& coef_prior_mean,
                                    const arma::vec& coef_prior_variance);

arma::vec sample_coef_blm_known1var_qr(const arma::vec& y, const arma::mat& X,
                                       const arma::vec& coef_prior_mean,
                                       const arma::vec& coef_prior_variance);

arma::vec sample_coef_blm_known1var_large_p_small_n(const arma::vec& y, const arma::mat& X,
                                                    const arma::vec& coef_prior_mean,
                                                    const arma::vec& coef_prior_variance);
#endif