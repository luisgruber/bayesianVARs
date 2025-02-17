#include <memory>
#include <lp_lib.h>
#include "utilities_cpp.h"

inline cube Sigma_chol_t_draws (
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

		ret.slice(r) = Sigma_chol.t();
	}
	
	return ret;
}

// [[Rcpp::export]]
Rcpp::List compute_parameter_transformations (
	const arma::cube& reduced_coefficients, //rows: lagged variables + intercept, columns: variables, slices: draws
	const arma::cube& factor_loadings, //rows: variables, columns: factors, slices: draws
	const arma::mat& U_vecs, //rows: entries of a (variables x variables)-upper triagonal matrix with ones on the diagonals, cols: draws
	const arma::mat& logvar_T, //cols: draws
	const Rcpp::List& restrictions
) {
	const uword n_variables = reduced_coefficients.n_cols;
	const uword n_posterior_draws = reduced_coefficients.n_slices;
	const uword n_factors = factor_loadings.n_cols;
	
	const bool include_facload = restrictions.containsElementNamed("facload");
	const bool include_B0_inv_t = restrictions.containsElementNamed("B0_inv_t");
	const bool include_B0 = restrictions.containsElementNamed("B0");
	const bool include_structural_coeff = restrictions.containsElementNamed("structural_coeff");
	const bool include_long_run_ir = restrictions.containsElementNamed("restrictions_long_run_ir");
	
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
	const cube B0_inv_t = Sigma_chol_t_draws(
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
		B0.slice(r) = inv(trimatu(B0_inv_t.slice(r).t()));
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

class Solver {
	public:
	virtual void add_constraints(mat& constraints, int constr_type, double rhs) = 0;
	virtual vec solve() = 0;
	virtual void recycle() = 0;
	virtual void print() = 0;
};

class LPSolver : public Solver {
	using make_lp_func_ptr= lprec*(*)(int, int);
	template<typename R, typename... T> using lp_func_ptr = R(*)(lprec*, T...);
	
	private:
	static constexpr double bigM = 10.0;
	
	//import routines from the R package "lpSolveAPI"
	make_lp_func_ptr make_lp = (make_lp_func_ptr)R_GetCCallable("lpSolveAPI", "make_lp");
	lp_func_ptr<void, int> set_verbose = (lp_func_ptr<void, int>)R_GetCCallable("lpSolveAPI", "set_verbose");
	lp_func_ptr<void> set_maxim = (lp_func_ptr<void>)R_GetCCallable("lpSolveAPI", "set_maxim");
	lp_func_ptr<bool, int, double> set_obj = (lp_func_ptr<bool, int, double>)R_GetCCallable("lpSolveAPI", "set_obj");
	lp_func_ptr<bool, int, double, double> set_bounds = (lp_func_ptr<bool, int, double, double>)R_GetCCallable("lpSolveAPI", "set_bounds");
	lp_func_ptr<bool, bool> set_add_rowmode = (lp_func_ptr<bool, bool>)R_GetCCallable("lpSolveAPI", "set_add_rowmode");
	lp_func_ptr<bool, int, double*, int*, int, double> add_constraintex = (lp_func_ptr<bool, int, double*, int*, int, double>)R_GetCCallable("lpSolveAPI", "add_constraintex");
	lp_func_ptr<int> lp_solve = (lp_func_ptr<int>)R_GetCCallable("lpSolveAPI", "solve");
	lp_func_ptr<bool, double*> get_variables = (lp_func_ptr<bool, double*>)R_GetCCallable("lpSolveAPI", "get_variables");
	lp_func_ptr<char*, int> get_statustext = (lp_func_ptr<char*, int>)R_GetCCallable("lpSolveAPI", "get_statustext");
	lp_func_ptr<void> print_lp = (lp_func_ptr<void>)R_GetCCallable("lpSolveAPI", "print_lp");
	lp_func_ptr<void> delete_lp = (lp_func_ptr<void>)R_GetCCallable("lpSolveAPI", "delete_lp");
	lp_func_ptr<bool, int, double*, int*> set_obj_fnex = (lp_func_ptr<bool, int, double*, int*>)R_GetCCallable("lpSolveAPI", "set_obj_fnex");
	lp_func_ptr<bool, int, int> resize_lp = (lp_func_ptr<bool, int, int>)R_GetCCallable("lpSolveAPI", "resize_lp");
	lp_func_ptr<bool, int, bool> set_binary = (lp_func_ptr<bool, int, bool>)R_GetCCallable("lpSolveAPI", "set_binary");
	lp_func_ptr<bool, int, char*> set_col_name = (lp_func_ptr<bool, int, char*>)R_GetCCallable("lpSolveAPI", "set_col_name");
	lp_func_ptr<void, int, int> set_presolve = NULL; //do not combine presolve and recycling!!
	lp_func_ptr<int> get_Nrows = (lp_func_ptr<int>)R_GetCCallable("lpSolveAPI", "get_Nrows");
	lp_func_ptr<int> get_Ncolumns = (lp_func_ptr<int>)R_GetCCallable("lpSolveAPI", "get_Ncolumns");
	
	lprec *lp;
	uword n; //dimension of x
	ivec x_cols;
	ivec U_cols;
	ivec y_cols;
	int n_initial_rows;

	public:
	LPSolver(const uword dim_x) : n(dim_x) {
		lp = (*make_lp)(0, 3*n);
		if (lp == NULL) throw std::bad_alloc();
		
		(*set_verbose)(lp, IMPORTANT);
		(*set_maxim)(lp);
		
		// we have variables x_i, U_i, y_i for i=1...n
		x_cols = regspace<ivec>(1, n);
		U_cols = x_cols + n;
		y_cols = U_cols + n;
		for (uword j = 0; j < n; j++) {
			std::string col_name("x");
			col_name += std::to_string(j+1);
			(*set_col_name)(lp, x_cols[j], col_name.data());
			col_name[0] = 'U';
			(*set_col_name)(lp, U_cols[j], col_name.data());
			col_name[0] = 'y';
			(*set_col_name)(lp, y_cols[j], col_name.data());
		}
		
		(*set_add_rowmode)(lp, TRUE);
				
		// maximize over sum of U
		vec ones(U_cols.n_elem, fill::ones);
		(*set_obj_fnex)(lp, U_cols.n_elem, ones.memptr(), U_cols.memptr());
		
		// set variable bounds
		for (uword j = 0; j < n; j++) {
			(*set_bounds)(lp, x_cols[j], -1, 1); //x_i in [-1, 1]
			(*set_bounds)(lp, U_cols[j], 0, 1); //U_i in [0, 1]
			(*set_binary)(lp, y_cols[j], TRUE); //y_i in {0, 1}
		}
		
		//setup initial constraints
		//x_i must not be within the interval [-U_i, U_i]
		//y_i switches if x_i is right to the interval or left to the interval
		//=> in an optimal solution, U_i is the greatest possible value of |x_i|
		std::array<double, 3> left_constraint  = { 1.0 /* x_i */, 1.0 /*U_i*/, -bigM /*y*/};
		std::array<double, 3> right_constraint = {-1.0 /* x_i */, 1.0 /*U_i*/,  bigM /*y*/};
		for (uword j = 0; j < n; j++) {
			std::array<int, 3> cols = {x_cols[j], U_cols[j], y_cols[j]};
			(*add_constraintex)(lp, cols.size(), left_constraint.data(), cols.data(), LE, 0);
			(*add_constraintex)(lp, cols.size(), right_constraint.data(), cols.data(), LE, bigM);
		}
		
		//recycling will remove all constraints added by `add_constraint`
		n_initial_rows = (*get_Nrows)(lp);
	}
	
	void add_constraints(mat& constraints, int constr_type, double rhs) override {
		constraints.each_row([&](rowvec& row) {
			const bool success = (*add_constraintex)(lp, n, row.memptr(), x_cols.memptr(), constr_type, rhs);
			if (success == FALSE) {
				throw std::runtime_error("Could not add constraints");
			}
		});
	}
	
	vec solve() override {
		(*set_add_rowmode)(lp, FALSE);
		int ret = (*lp_solve)(lp);
		if(ret != 0) {
			print();
			throw std::logic_error((*get_statustext)(lp, ret));
		};
		vec lp_vars_store(3*n);
		(*get_variables)(lp, lp_vars_store.memptr());
		return lp_vars_store.head_rows(n); // only return x
	}
	
	void recycle() override {
		// delete all non-initial constraints
		(*resize_lp)(lp, n_initial_rows, (*get_Ncolumns)(lp));
		(*set_add_rowmode)(lp, TRUE);
	}
	
	void print() override {
		(*set_add_rowmode)(lp, FALSE);
		(*print_lp)(lp);
	}
		
	~LPSolver() {
		(*delete_lp)(lp);
	}
	
};

class RandomizedSolver : public Solver {
	private:
		mat zero_restrictions;
		mat sign_restrictions;
		const uword dim_x;
		const double tol;
		const uword max_attempts;
	public:
	RandomizedSolver(const uword dim_x, const uword max_attempts, const double tol) :
		zero_restrictions(mat(0, dim_x)),
		sign_restrictions(mat(0, dim_x)),
		dim_x(dim_x), tol(tol), max_attempts(max_attempts)
	{}
	
	void add_constraints(mat& constraints, int constr_type, double rhs) override {
		if ( constr_type == EQ && rhs == 0) {
			zero_restrictions = join_vert(zero_restrictions, constraints);
		}
		else if ( constr_type == GE && rhs == 0) {
			sign_restrictions = join_vert(sign_restrictions, constraints);
		}
		else throw std::domain_error("Constraint type not implemented");
	};
	virtual vec solve() override {
		mat zero_restrictions_solution_space;
		if (zero_restrictions.n_rows == 0) {
			zero_restrictions_solution_space = eye(dim_x, dim_x);
		} else {
			zero_restrictions_solution_space = null(zero_restrictions, tol);
		}
		for (uword i = 0; i < max_attempts; i++) {
			// draw random unit-length vector
			vec v(zero_restrictions_solution_space.n_cols);
    		v.imbue(R::norm_rand);
    		v = normalise(v);
    		
    		const vec x_candidate = zero_restrictions_solution_space*v;
    		if (all(sign_restrictions * x_candidate >= 0-tol))
    			return x_candidate;
		}
		throw std::runtime_error("Max attempts reached while searching for solution");
	};
	virtual void recycle() override {
		zero_restrictions = mat(0, zero_restrictions.n_cols);
		sign_restrictions = mat(0, sign_restrictions.n_cols);
	};
	virtual void print() override {
		Rcout << "Zero restrictions:" << endl << zero_restrictions <<
		         "Sign restrictions:" << endl << sign_restrictions;
	};
};

// [[Rcpp::export]]
arma::cube find_rotation_cpp(
	const arma::field<arma::cube>& parameter_transformations, //each field element: rows: transformation size, cols: variables, slices: draws
	const arma::field<Rcpp::NumericMatrix>& restriction_specs, //each field element: rows: transformation size, cols: variables
	const std::string& solver_option = "randomized", // "randomized" or "lp"
	const arma::uword randomized_max_attempts = 100, // ignored when solver is "lp"
	const double tol = 1e-6
) {
	//algorithm from RUBIO-RAMÍREZ ET AL. (doi: 10.1111/j.1467-937X.2009.00578.x)
	if (restriction_specs.n_elem != parameter_transformations.n_elem) {
		throw std::logic_error("Number of restrictions does not match number of parameter transformations.");
	}

	const uword n_variables = parameter_transformations(0).n_cols;
	const uword n_posterior_draws = parameter_transformations(0).n_slices;
	cube rotation(n_variables, n_variables, n_posterior_draws, fill::none);

	std::unique_ptr<Solver> solver;
	if (solver_option == "randomized") {
		solver.reset(new RandomizedSolver(n_variables, randomized_max_attempts, tol));
	}
	else if (solver_option == "lp") {
		solver.reset(new LPSolver(n_variables));
	}
	else throw std::domain_error("Unknown solver option");
	
	//field rows: tranformations, field cols: cols of the transformation
	//each field element: rows: number of restrictions, cols: transformation size
	field<mat> zero_restrictions(restriction_specs.n_elem, n_variables);
	field<mat> sign_restrictions(restriction_specs.n_elem, n_variables);
	uvec n_zero_restrictions(n_variables, fill::zeros);
	uvec n_sign_restrictions(n_variables, fill::zeros);
	for (uword i = 0; i < restriction_specs.n_elem; i++) {
		for (uword j = 0; j < n_variables; j++) {
			const NumericMatrix::ConstColumn column_restriction_spec = restriction_specs(i).column(j);
			zero_restrictions(i, j) = construct_zero_restriction(column_restriction_spec);
			sign_restrictions(i, j) = construct_sign_restriction(column_restriction_spec);
			n_zero_restrictions(j) += zero_restrictions(i, j).n_rows; //rank = n_rows by construction!
			n_sign_restrictions(j) += sign_restrictions(i, j).n_rows; //rank = n_rows by construction!
		}
	}
	
	//iterate over columns in order or descending rank.
	//since each column must be orthogonal to the ones that came before,
	//we start with the column with the most restrictions	
	uvec col_order = sort_index(2 * n_zero_restrictions + n_sign_restrictions, "descend");
	for (uword r = 0; r < n_posterior_draws; r++) {
		if (r % 512 == 0) Rcpp::checkUserInterrupt();
		for (uword j_index = 0; j_index < n_variables; solver->recycle(), j_index++) {
			const uword j = col_order[j_index];
			
			// the column j of the rotation matrix must be orthogonal to columns that came before
			mat previous_columns = rotation.slice(r).cols(col_order.head(j_index)).t();
			solver->add_constraints(previous_columns, EQ, 0);
			
			if (n_zero_restrictions(j) > 0) {
				for (uword i = 0; i < parameter_transformations.n_elem; i++) {
					mat zero_constraints = zero_restrictions(i, j) * parameter_transformations(i).slice(r);
					solver->add_constraints(zero_constraints, EQ, 0);
				}
			}
			
			if (n_sign_restrictions(j) > 0) {
				for (uword i = 0; i < parameter_transformations.n_elem; i++) {
					mat sign_constraints = sign_restrictions(i, j) * parameter_transformations(i).slice(r);
					solver->add_constraints(sign_constraints, GE, 0);
				}
			}
			
			const vec p_j = solver->solve().clean(tol);
			if (p_j.is_zero()) {
				//zero was the optimal solution
				Rcerr << "Cannot satisfy restrictions for posterior sample: #" << r
					<< ", column:" << j+1
					<< " (" << j_index+1 << "-th in order of rank)" << endl;
				solver->print();
				throw std::logic_error("Could not satisfy restrictions");
			}
			rotation.slice(r).col(j) = normalise(p_j);
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
		cube irf(n_variables, n_shocks, ahead+1);
		
		//compute the responses to the shocks at t=0
		mat rotated_shocks = rotation.n_slices > 0 ? rotation.slice(r) * shocks : shocks;
		if (is_factor_model) {
			vec sqrt_V_t = exp(logvar_t.col(r).head(factor_loadings.n_cols) / 2);
			rotated_shocks.each_col() %= sqrt_V_t;
			irf.slice(0) = factor_loadings.slice(r) * rotated_shocks;
		}
		else {
		    mat U(n_variables, n_variables, fill::eye);
			U(upper_indices) = U_vecs.col(r);

			vec sqrt_D_t = exp(logvar_t.col(r) / 2.0);
			rotated_shocks.each_col() %= sqrt_D_t;
			irf.slice(0) = solve(trimatl(U.t()), rotated_shocks);
		}
		
		// compute how the shocks propagate using the reduced form coeffs
		mat current_predictors(n_shocks, coefficients.n_rows, fill::zeros);
		current_predictors.head_cols(n_variables) = irf.slice(0).t(); //set lag zero
		for (uword t = 1; t <= ahead; t++) {
			mat new_predictiors = current_predictors * coefficients.slice(r);
			irf.slice(t) = new_predictiors.t();
			// shift everything and make predictions the new predictors at lag zero
			shift_and_insert(current_predictors, new_predictiors);
		}
		ret(r) = irf;
	}
	
	return ret;
}

// [[Rcpp::export]]
arma::ivec irf_bayes_optimal_order(arma::field<arma::cube>& irf) {
	const uword n_posterior_draws = irf.n_elem;
	vec losses(n_posterior_draws, fill::zeros);
	for (uword r = 0; r < n_posterior_draws; r++) {
		for (uword r_other = 0; r_other < r; r_other++) {
			const double loss = accu(abs(irf(r) - irf(r_other)));
			losses(r) += loss;
			losses(r_other) += loss;
		}
	}
	return conv_to<ivec>::from(sort_index(losses, "ascend"));
}

// [[Rcpp::export]]
arma::cube irf_from_true_parameters(
	arma::mat true_structural_matrix,
	arma::mat true_reduced_coeff,
	arma::uword ahead
) {
	const uword n_variables = true_structural_matrix.n_rows;
	const uword n_shocks = true_structural_matrix.n_cols;
	
	cube irf(n_variables, n_shocks, ahead+1);
	irf.slice(0) = inv(true_structural_matrix).t();
	
	mat current_predictors(n_shocks, true_reduced_coeff.n_rows, fill::zeros);
	current_predictors.head_cols(n_variables) = irf.slice(0).t();
	for (uword t = 1; t <= ahead; t++) {
		mat new_predictiors = current_predictors * true_reduced_coeff;
		irf.slice(t) = new_predictiors.t();
		shift_and_insert(current_predictors, new_predictiors);
	}
	return irf;
}
