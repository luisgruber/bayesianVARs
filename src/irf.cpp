#include "utilities_cpp.h"

inline arma::cube Sigma_chol_draws (
	const uword n_variables,
	const uword n_posterior_draws,
	const arma::cube& factor_loadings, //rows: variables, columns: factors, slices: draws
	const arma::mat& U_vecs, //rows: entries of a (variables x variables)-upper triagonal matrix with ones on the diagonals, cols: draws
	const arma::mat& logvar_T
) {
	const bool is_factor_model = factor_loadings.n_cols > 0;
	
	arma::cube ret(n_variables, n_variables, n_posterior_draws, arma::fill::none);
	for (uword r = 0; r < n_posterior_draws; r++) {
		// obtain the cholesky factorization of the covariance matrix
		arma::mat Sigma_unused, Sigma_chol, facload_mat;
		arma::vec u_vec;
		if(is_factor_model){
			facload_mat = factor_loadings.slice(r);
		} else {
			u_vec = U_vecs.col(r);
		}
		build_sigma(Sigma_unused, Sigma_chol, is_factor_model, facload_mat,
				logvar_T.col(r).as_row(), factor_loadings.n_cols,
				n_variables, u_vec, true);

		ret.slice(r) = Sigma_chol;
	}
	
	return ret;
}

// [[Rcpp::export]]
Rcpp::List compute_parameter_transformations (
	const arma::cube& reduced_coefficients, //rows: lagged variables + intercept, columns: variables, slices: draws
	const arma::cube& factor_loadings, //rows: variables, columns: factors, slices: draws
	const arma::mat& U_vecs, //rows: entries of a (variables x variables)-upper triagonal matrix with ones on the diagonals, cols: draws
	const arma::mat& logvar_T,
	
	const bool include_B0_inv = true,
	const bool include_B0 = false,
	const bool include_structural_coeff = false,
	const bool include_long_run_ir = false
) {
	const uword n_variables = reduced_coefficients.n_cols;
	const uword n_posterior_draws = reduced_coefficients.n_slices;
	
	const arma::cube Sigma_chol = Sigma_chol_draws(
		n_variables,
		n_posterior_draws,
		factor_loadings,
		U_vecs,
		logvar_T
	);
	
	Rcpp::List ret = Rcpp::List::create();
	if (include_B0_inv) ret["B0_inv"] = Sigma_chol;
  	
  	if (include_B0 || include_structural_coeff || include_long_run_ir) {
  		arma::cube B0(n_variables, n_variables, n_posterior_draws, arma::fill::none);
  		for (uword r = 0; r < n_posterior_draws; r++) {
			B0.slice(r) = arma::inv(arma::trimatu(Sigma_chol.slice(r)));
		}
		if (include_B0) ret["B0"] = B0;
		
		if (include_structural_coeff || include_long_run_ir) {
			arma::cube structural_coeff(reduced_coefficients);
			for (uword r = 0; r < n_posterior_draws; r++) {
				structural_coeff.slice(r) *= B0.slice(r);
			}
			if (include_structural_coeff) ret["structural_coeff"] = structural_coeff;
			
			if (include_long_run_ir) {
				arma::cube IR_inf(n_variables, n_variables, n_posterior_draws, arma::fill::none);
				for (uword r = 0; r < n_posterior_draws; r++) {
					arma::mat sum_of_lags(n_variables, n_variables);
					for (uword l = 0; l + n_variables < structural_coeff.n_rows; l += n_variables) {
						sum_of_lags += structural_coeff.slice(r).rows(l, l + n_variables - 1);
					}
					IR_inf.slice(r) = arma::inv((B0.slice(r) - sum_of_lags).t());
				}
				ret["IR_inf"] = IR_inf;
			}
		}
  	}
	
	return ret;
}

// [[Rcpp::export]]
arma::cube find_rotation_cpp(
	const arma::field<arma::cube>& parameter_transformations, //each field element: rows: transformation size, cols: variables, slices: draws
	const arma::field<arma::cube>& restrictions //each: rows: number of restrictions, cols: transformation size, slices: variables
) {
	const uword n_variables = parameter_transformations(0).n_cols;
	const uword n_posterior_draws = parameter_transformations(0).n_slices;
	arma::cube rotation(n_variables, n_variables, n_posterior_draws, arma::fill::none);

	for (uword r = 0; r < n_posterior_draws; r++) {
		for (uword j = 0; j < n_variables; j++) {
			arma::mat Q(rotation.slice(r).head_cols(j).t());
			for (uword i = 0; i < parameter_transformations.n_elem; i++) {
				Q.insert_rows(Q.n_rows, restrictions(i).slice(j) * parameter_transformations(i).slice(r));
			}
			arma::mat nullspaceQ = arma::null(Q);
			if (nullspaceQ.n_cols == 0) {
				throw std::logic_error("Could not satisfy restrictions");
			}
			rotation.slice(r).col(j) = nullspaceQ.col(0);
		}
	}
	return rotation;
}

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
	const arma::mat& U_vecs, //rows: entries of a (variables x variables)-upper triagonal matrix with ones on the diagonals, cols: draws
	const arma::colvec& shock, //columns: how much each of the factors is shocked
	const arma::uword ahead //how far to predict ahead
) {
	const uword n_variables = coefficients.n_cols;
	const uword n_posterior_draws = coefficients.n_slices;
	const bool is_factor_model = factor_loadings.n_cols > 0;

	arma::cube irf(ahead+1, n_variables, n_posterior_draws, arma::fill::none);
	
	// trace out n_posterior_draws paths
	arma::mat current_predictors(n_posterior_draws, coefficients.n_rows, arma::fill::zeros);
	// all paths start with the some shock at t=0 to the variables
	const arma::uvec upper_indices = arma::trimatu_ind(arma::size(n_variables, n_variables), 1);
	for (uword r = 0; r < n_posterior_draws; r++) {
		rowvec y_shock;
		if (is_factor_model) {
			y_shock = (factor_loadings.slice(r) * shock).t();
		}
		else {
		    arma::mat U(n_variables, n_variables, arma::fill::eye);
			U(upper_indices) = U_vecs.col(r);
			U = U.t();
			y_shock = solve(arma::trimatl(U), shock).t();
		}

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
