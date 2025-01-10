#include "utilities_cpp.h"

inline cube Sigma_chol_draws (
	const uword n_variables,
	const uword n_posterior_draws,
	const cube& factor_loadings, //rows: variables, columns: factors, slices: draws
	const mat& U_vecs, //rows: entries of a (variables x variables)-upper triagonal matrix with ones on the diagonals, cols: draws
	const mat& logvar_T
) {
	const bool is_factor_model = factor_loadings.n_cols > 0;
	
	cube ret(n_variables, n_variables, n_posterior_draws, fill::none);
	for (uword r = 0; r < n_posterior_draws; r++) {
		// obtain the cholesky factorization of the covariance matrix
		mat Sigma_unused, Sigma_chol, facload_mat;
		vec u_vec;
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
	const arma::mat& logvar_T, //cols: draws
	
	const bool include_facload = false,
	const bool include_B0_inv_t = false,
	const bool include_B0 = false,
	const bool include_structural_coeff = false,
	const bool include_long_run_ir = false
) {
	const uword n_variables = reduced_coefficients.n_cols;
	const uword n_posterior_draws = reduced_coefficients.n_slices;
	const uword n_factors = factor_loadings.n_cols;
	
	Rcpp::List ret = Rcpp::List::create();
	
	if (include_facload) {
		if (factor_loadings.n_elem == 0)
			throw std::logic_error("facload can only be included for the factor model");
		cube factor_loadings_T (factor_loadings);
		for (uword r = 0; r < n_posterior_draws; r++) {
			factor_loadings_T.slice(r).each_row() %= exp(logvar_T.col(r).head(n_factors)/2).t();
		}
		ret.push_back(factor_loadings_T, "facload");
	}

	if (!(include_B0_inv_t || include_B0 || include_structural_coeff || include_long_run_ir)) return ret;
	const cube B0_inv_t = Sigma_chol_draws(
		n_variables,
		n_posterior_draws,
		factor_loadings,
		U_vecs,
		logvar_T
	);
	if (include_B0_inv_t) ret.push_back(B0_inv_t, "B0_inv_t");

	if (!(include_B0 || include_structural_coeff || include_long_run_ir)) return ret;
	cube B0(n_variables, n_variables, n_posterior_draws, fill::none);
	for (uword r = 0; r < n_posterior_draws; r++) {
		B0.slice(r) = inv(trimatl(B0_inv_t.slice(r).t()));
	}
	if (include_B0) ret.push_back(B0, "B0");

	if (!(include_structural_coeff || include_long_run_ir)) return ret;
	cube structural_coeff(reduced_coefficients);
	for (uword r = 0; r < n_posterior_draws; r++) {
		structural_coeff.slice(r) *= B0.slice(r);
	}
	if (include_structural_coeff) ret.push_back(structural_coeff, "structural_coeff");

	if (!include_long_run_ir) return ret;
	cube IR_inf(n_variables, n_variables, n_posterior_draws, fill::none);
	for (uword r = 0; r < n_posterior_draws; r++) {
		mat sum_of_lags(n_variables, n_variables);
		for (uword l = 0; l + n_variables < structural_coeff.n_rows; l += n_variables) {
			sum_of_lags += structural_coeff.slice(r).rows(l, l + n_variables - 1);
		}
		IR_inf.slice(r) = inv((B0.slice(r) - sum_of_lags).t());
	}
	ret.push_back(IR_inf, "IR_inf");
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

mat construct_zero_restriction(const NumericMatrix::ConstColumn& spec) {
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

mat construct_sign_restriction (const NumericMatrix::ConstColumn& spec) {
	const uword n_sign_restrictions = count(spec != 0);

	mat sign_restriction_matrix(n_sign_restrictions, spec.size());
	uword i = 0;
	for (R_xlen_t j = 0; j < spec.size(); j++) {
		if ((spec[j] > 0) == TRUE) {
			sign_restriction_matrix(i++, j) = 1;
		}
		else if ((spec[j] < 0) == TRUE) {
			sign_restriction_matrix(i++, j) = -1;
		}
	}
	return sign_restriction_matrix;
}

// [[Rcpp::export]]
arma::cube find_rotation_cpp(
	const arma::field<arma::cube>& parameter_transformations, //each field element: rows: transformation size, cols: variables, slices: draws
	const arma::field<Rcpp::NumericMatrix>& restriction_specs, //each field element: rows: transformation size, cols: variables
	const double tolerance = 0.0,
	const double sign_epsilon = 1e-6
) {
	//algorithm from RUBIO-RAMÍREZ ET AL. (doi: 10.1111/j.1467-937X.2009.00578.x)

	if (restriction_specs.n_elem != parameter_transformations.n_elem) {
		throw std::logic_error("Number of restrictions does not match number of parameter transformations.");
	}
	if (!(sign_epsilon > 0)) {
		throw std::logic_error("sign_epsilon must be a positive number");
	}

	const uword n_variables = parameter_transformations(0).n_cols;
	const uword n_posterior_draws = parameter_transformations(0).n_slices;
	cube rotation(n_variables, n_variables, n_posterior_draws, fill::none);

	//field rows: tranformations, field cols: cols of the transformation
	//each field element: rows: number of restrictions, cols: transformation size
	field<mat> zero_restrictions(restriction_specs.n_elem, n_variables);
	field<mat> sign_restrictions(restriction_specs.n_elem, n_variables);
	uvec n_sign_restrictions(n_variables, fill::zeros);
	for (uword i = 0; i < restriction_specs.n_elem; i++) {
		for (uword j = 0; j < n_variables; j++) {
			const NumericMatrix::ConstColumn column_restriction_spec = restriction_specs(i).column(j);
			zero_restrictions(i, j) = construct_zero_restriction(column_restriction_spec);
			sign_restrictions(i, j) = construct_sign_restriction(column_restriction_spec);
			n_sign_restrictions(j) += sign_restrictions(i, j).n_rows;
		}
	}

	for (uword r = 0; r < n_posterior_draws; r++) {
		for (uword j = 0; j < n_variables; j++) {
			colvec p_j; // find the j-th column of the rotation matrix
			
			mat Q_zero(rotation.slice(r).head_cols(j).t());
			for (uword i = 0; i < parameter_transformations.n_elem; i++) {
				Q_zero.insert_rows(0, zero_restrictions(i, j) * parameter_transformations(i).slice(r));
			}
			mat nullspace_Q_zero;
			if (Q_zero.n_rows == 0) {
				nullspace_Q_zero = mat(n_variables, n_variables, fill::eye);
			} else {
				nullspace_Q_zero = null(Q_zero, tolerance);
				if (nullspace_Q_zero.n_cols == 0) {
					throw std::logic_error("Could not satisfy zero restrictions. Increase the tolerance for approximate results.");
				}
			}

			if (n_sign_restrictions(j) > 0) {
				//find the vector in the nullspace of Q_zero which satisfies the sign restrictions
				mat Q_sign(0, n_variables);
				for (uword i = 0; i < parameter_transformations.n_elem; i++) {
					Q_sign.insert_rows(0, sign_restrictions(i, j) * parameter_transformations(i).slice(r));
				}
				const vec small_positive_vector(Q_sign.n_rows, fill::value(sign_epsilon));
				p_j = nullspace_Q_zero * solve(Q_sign * nullspace_Q_zero, small_positive_vector);
				p_j = normalise(p_j);
			}
			else {
				//any vector from the null space is fine
				//vector with corresponding to the smallest singular value should be the last one
				//however this is an not guranteed by the public armadillo API!
				p_j = nullspace_Q_zero.col(nullspace_Q_zero.n_cols - 1);	
			}

			rotation.slice(r).col(j) = p_j;
		}
	}
	return rotation;
}

inline void shift_and_insert(
	mat& X, //the columns of X should be y1,y2,y3, y1.l1,y2.l1,y3.l1,...,1
	const mat& new_y //what to insert in y1,y2,y3
) {
	for (uword i = X.n_cols-2; new_y.n_cols <= i; i--) {
		X.col(i) = X.col(i-new_y.n_cols);
	}
	X.head_cols(new_y.n_cols) = new_y;
}

// [[Rcpp::export]]
arma::field<arma::cube> irf_cpp(
	const arma::cube& coefficients, //rows: lagged variables + intercept, columns: variables, slices: draws
	const arma::cube& factor_loadings, //rows: variables, columns: factors, slices: draws
	const arma::mat& U_vecs, //rows: entries of a (variables x variables)-upper triagonal matrix with ones on the diagonals, cols: draws
	const arma::mat& logvar_t, //rows: log variances, cols: draws
	const arma::mat& shocks, //rows: dim shock, cols: shocks
	const arma::uword ahead, //how far to predict ahead
	const Rcpp::Nullable<Rcpp::NumericMatrix> rotation_ = R_NilValue //rows: variables, cols: dim shock, slices: draws
) {
	const uword n_shocks = shocks.n_cols;
	const uword n_variables = coefficients.n_cols;
	const uword n_posterior_draws = coefficients.n_slices;
	const bool is_factor_model = factor_loadings.n_cols > 0;
	const uvec upper_indices = trimatu_ind(size(n_variables, n_variables), 1);
	cube rotation;
	if (rotation_.isNotNull()) {
		rotation = Rcpp::as<cube>(rotation_);
	}

	field<cube> ret(n_posterior_draws);
	for (uword r = 0; r < n_posterior_draws; r++) {
		cube irf(n_shocks, n_variables, ahead+1);
		
		//compute the responses to the shocks at t=0
		mat rotated_shocks = rotation.n_slices > 0 ? rotation.slice(r) * shocks : shocks;
		if (is_factor_model) {
			vec sqrt_V_t = exp(logvar_t.col(r).head(factor_loadings.n_cols) / 2);
			rotated_shocks.each_col() %= sqrt_V_t;
			irf.slice(0) = (factor_loadings.slice(r) * rotated_shocks).t();
		}
		else {
		    mat U(n_variables, n_variables, fill::eye);
			U(upper_indices) = U_vecs.col(r);

			vec sqrt_D_t = exp(logvar_t.col(r) / 2.0);
			rotated_shocks.each_col() %= sqrt_D_t;
			irf.slice(0) = solve(trimatl(U.t()), rotated_shocks).t();
		}
		
		// compute how the shocks propagate using the reduced form coeffs
		mat current_predictors(n_shocks, coefficients.n_rows, fill::zeros);
		current_predictors.head_cols(n_variables) = irf.slice(0); //set lag zero
		for (uword t = 1; t <= ahead; t++) {
			irf.slice(t) = current_predictors * coefficients.slice(r);
			// shift everything and make predictions the new predictors at lag zero
			shift_and_insert(current_predictors, irf.slice(t));
		}
		ret(r) = irf;
	}
	
	return ret;
}
