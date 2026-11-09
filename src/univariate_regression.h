#ifndef UNIVARIATE_REGRESSION_H
#define UNIVARIATE_REGRESSION_H

arma::vec sample_coef_blm_known1var(const arma::vec& y, const arma::mat& X,
                                    const arma::vec& coef_prior_mean,
                                    const arma::vec& coef_prior_variance,
                                    const bool huge);
#endif