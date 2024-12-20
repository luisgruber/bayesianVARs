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
	
	const bool include_B0_inv_t = true,
	const bool include_B0 = false,
	const bool include_structural_coeff = false,
	const bool include_long_run_ir = false
) {
	const uword n_variables = reduced_coefficients.n_cols;
	const uword n_posterior_draws = reduced_coefficients.n_slices;
	
	const arma::cube B0_inv_t = Sigma_chol_draws(
		n_variables,
		n_posterior_draws,
		factor_loadings,
		U_vecs,
		logvar_T
	);
	
	Rcpp::List ret = Rcpp::List::create();
	if (include_B0_inv_t) ret.push_back(B0_inv_t, "B0_inv_t");
  	
  	if (include_B0 || include_structural_coeff || include_long_run_ir) {
  		arma::cube B0(n_variables, n_variables, n_posterior_draws, arma::fill::none);
  		for (uword r = 0; r < n_posterior_draws; r++) {
			B0.slice(r) = arma::inv(arma::trimatl(B0_inv_t.slice(r).t()));
		}
		if (include_B0) ret.push_back(B0, "B0");
		
		if (include_structural_coeff || include_long_run_ir) {
			arma::cube structural_coeff(reduced_coefficients);
			for (uword r = 0; r < n_posterior_draws; r++) {
				structural_coeff.slice(r) *= B0.slice(r);
			}
			if (include_structural_coeff) ret.push_back(structural_coeff, "structural_coeff");
			
			if (include_long_run_ir) {
				arma::cube IR_inf(n_variables, n_variables, n_posterior_draws, arma::fill::none);
				for (uword r = 0; r < n_posterior_draws; r++) {
					arma::mat sum_of_lags(n_variables, n_variables);
					for (uword l = 0; l + n_variables < structural_coeff.n_rows; l += n_variables) {
						sum_of_lags += structural_coeff.slice(r).rows(l, l + n_variables - 1);
					}
					IR_inf.slice(r) = arma::inv((B0.slice(r) - sum_of_lags).t());
				}
				ret.push_back(IR_inf, "IR_inf");
			}
		}
  	}
	
	return ret;
}

uword count(const Rcpp::LogicalVector& vec) {
	// vec can contain the values TRUE, FALSE and notably NA_LOGICAL
	uword counter = 0;
	for(int i = 0; i < vec.size(); i++) {
		if (vec[i] == TRUE) counter++;
	}
	return counter;
}

arma::mat construct_zero_restriction(const NumericMatrix::ConstColumn& spec) {
	const LogicalVector spec_is_zero = (spec == 0);
	const uword n_zero_restrictions = count(spec_is_zero);

	mat zero_restriction_matrix(n_zero_restrictions, spec.size());
	uword i = 0;
	for (uword j = 0; j < spec_is_zero.size(); j++) {
		if (spec_is_zero[j] == TRUE) {
			zero_restriction_matrix(i++, j) = 1;
		}
	}
	return zero_restriction_matrix;
}

vec calculate_sign_restriction_scores(
	const NumericMatrix::ConstColumn& spec, //rows: transformation size
	const mat& rotated_params //rows: transformation size, cols: candidates for p_j
) {
	vec scores(rotated_params.n_cols);
	for (uword i = 0; i < rotated_params.n_rows && !NumericMatrix::is_na(spec[i]); i++) {
		for (uword j = 0; j < rotated_params.n_cols; j++) {
			scores[j] += std::copysign(spec[i], rotated_params(i, j) * spec[i]);
			// if the rotated params and its sign specification spec[i] have the same sign
			// the score increases by spec[i] otherwise the score will decrease by that amount.
			// if the score is negative overall, we will choose -p_j instead of p_j.
		}
	}
	return scores;
}


// [[Rcpp::export]]
arma::cube find_rotation_cpp(
	const arma::field<arma::cube>& parameter_transformations, //each field element: rows: transformation size, cols: variables, slices: draws
	const arma::field<Rcpp::NumericMatrix>& restriction_specs, //each field element: rows: transformation size, cols: variables
	const double tolerance = 0.0
) {
	//algorithm from RUBIO-RAMÍREZ ET AL. (doi: 10.1111/j.1467-937X.2009.00578.x)

	if (restriction_specs.n_elem != parameter_transformations.n_elem) {
		throw std::logic_error("Number of restrictions does not match number of parameter transformations.");
	}

	const uword n_variables = parameter_transformations(0).n_cols;
	const uword n_posterior_draws = parameter_transformations(0).n_slices;
	arma::cube rotation(n_variables, n_variables, n_posterior_draws, arma::fill::none);

	//field rows: tranformations, field cols: cols of the transformation
	//each field element: rows: number of restrictions, cols: transformation size
	arma::field<arma::mat> zero_restrictions(restriction_specs.n_elem, n_variables);
	uword n_sign_restrictions = 0;
	for (uword i = 0; i < restriction_specs.n_elem; i++) {
		for (uword j = 0; j < n_variables; j++) {
			const NumericMatrix::ConstColumn column_restriction_spec = restriction_specs(i).column(j);
			zero_restrictions(i, j) = construct_zero_restriction(column_restriction_spec);
			n_sign_restrictions += count(column_restriction_spec != 0);
		}
	}

	for (uword r = 0; r < n_posterior_draws; r++) {
		for (uword j = 0; j < n_variables; j++) {
			arma::mat Q_tilde(rotation.slice(r).head_cols(j).t());
			for (uword i = 0; i < parameter_transformations.n_elem; i++) {
				Q_tilde.insert_rows(0, zero_restrictions(i, j) * parameter_transformations(i).slice(r));
			}
			arma::mat nullspace_Q_tilde = arma::null(Q_tilde, tolerance);
			if (nullspace_Q_tilde.n_cols == 0) {
				throw std::logic_error("Could not satisfy restrictions. Increase the tolerance for approximate results.");
			}

			colvec p_j;
			if (n_sign_restrictions > 0) {
				//find the vector in the nullspace of Q which scores best in the sign restrictions
				vec sign_restriction_scores(nullspace_Q_tilde.n_cols, arma::fill::zeros);
				for (uword i = 0; i < parameter_transformations.n_elem; i++) {
					const NumericMatrix::ConstColumn column_restriction_spec = restriction_specs(i).column(j);
					const mat rotated_params = parameter_transformations(i).slice(r) * nullspace_Q_tilde;
					sign_restriction_scores += calculate_sign_restriction_scores(column_restriction_spec, rotated_params);
				}
				uword index_of_best_score = abs(sign_restriction_scores).index_max();
				p_j = nullspace_Q_tilde.col(index_of_best_score);
				if (sign_restriction_scores[index_of_best_score] < 0) {
					p_j = -p_j;
				}
			}
			else {
				//any vector from the null space is fine
				//vector with corresponding to the smallest singular value should be the last one
				//however this is an not guranteed by the public armadillo API!
				p_j = nullspace_Q_tilde.col(nullspace_Q_tilde.n_cols - 1);
			}

			rotation.slice(r).col(j) = p_j;
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
	const arma::mat& logvar_t, //rows: log variances, cols: draws
	const arma::colvec& shock, //columns: how much each of the factors is shocked
	const arma::uword ahead, //how far to predict ahead
	const Rcpp::Nullable<Rcpp::NumericMatrix> rotation_ = R_NilValue //rows: variables, cols: dim shock, slices: draws
) {
	const uword n_variables = coefficients.n_cols;
	const uword n_posterior_draws = coefficients.n_slices;
	const bool is_factor_model = factor_loadings.n_cols > 0;
	arma::cube rotation;
	if (rotation_.isNotNull()) {
		rotation = Rcpp::as<arma::cube>(rotation_);
	}

	arma::cube irf(ahead+1, n_variables, n_posterior_draws, arma::fill::none);
	
	// trace out n_posterior_draws paths
	arma::mat current_predictors(n_posterior_draws, coefficients.n_rows, arma::fill::zeros);
	// all paths start with the some shock at t=0 to the variables
	const arma::uvec upper_indices = arma::trimatu_ind(arma::size(n_variables, n_variables), 1);
	for (uword r = 0; r < n_posterior_draws; r++) {
		arma::rowvec y_shock;
		arma::colvec rotated_shock = rotation.n_slices > 0 ? rotation.slice(r) * shock : shock;

		if (is_factor_model) {
			y_shock = (factor_loadings.slice(r) * rotated_shock).t();
		}
		else {
		    arma::mat U(n_variables, n_variables, arma::fill::eye);
			U(upper_indices) = U_vecs.col(r);

			arma::vec sqrt_D_t = arma::exp(logvar_t.col(r) / 2.0);
			y_shock = solve(arma::trimatl(U.t()), sqrt_D_t % rotated_shock).t();
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
