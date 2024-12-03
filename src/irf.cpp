#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

inline void shift_and_insert(
	arma::mat& X, //the columns of X should be y1,y2,y3, y1.l1,y2.l1,y3.l1,...,1
	const arma::mat& new_y //what to insert in y1,y2,y3
) {
	for (uword i = X.n_cols-2; new_y.n_cols <= i; i--) {
		X.col(i) = X.col(i-new_y.n_cols);
	}
	X.head_cols(new_y.n_cols) = new_y;
}

// [[Rcpp::export]]
arma::cube irf_cpp(
	const arma::cube& coefficients, //rows: lagged variables + intercept, columns: variables, slices: draws
	const arma::cube& factor_loadings, //rows: variables, columns: factors, slices: draws
	const arma::colvec& f_shock, //columns: how much each of the factors is shocked
	const arma::uword ahead //how far to predict ahead
) {
	const uword n_variables = coefficients.n_cols;
	const uword n_posterior_draws = coefficients.n_slices;
	arma::cube irf(ahead+1, n_variables, n_posterior_draws, arma::fill::none);
	
	// trace out n_posterior_draws paths
	arma::mat current_predictors(n_posterior_draws, coefficients.n_rows, arma::fill::zeros);
	// all paths start with the some shock at t=0 to the variables
	// the shock is to the factors, so the uncertainty of the factor loadings
	// translates to uncertainty to the shock to the variables
	for (uword r = 0; r < n_posterior_draws; r++) {
		const rowvec y_shock = (factor_loadings.slice(r) * f_shock).t();
		irf.slice(r).row(0) = y_shock;
		current_predictors.head_cols(n_variables).row(r) = y_shock;
	}
	
	for (uword t = 1; t <= ahead; t++) {
		for(uword r = 0; r < n_posterior_draws; r++) {
			irf.slice(r).row(t) = current_predictors.row(r) * coefficients.slice(r);
		}
		// shift everything and make predictions the new predictors at lag zero
		shift_and_insert(current_predictors, irf.row_as_mat(t));
	}
	
	return irf;
}
