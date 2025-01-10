#' Impulse response functions
#'
#' Effect of a shock on the factors or variables over time
#' 
#' If a factor model was used, then the number of shocks is equal to the number
#' of factors.
#'
#' If the Cholesky model was used, then the number of shocks is equal to the
#' number of variables / time series.
#'
#' @param x An object of type `bayesianVARs_bvar`.
#'
#' @return Returns a `bayesianVARs_irf` object
#'
#' @examples
#' train_data <- 100 * usmacro_growth[,c("GDPC1", "PCECC96", "GPDIC1", "AWHMAN", "GDPCTPI", "CES2000000008x", "FEDFUNDS", "GS10", "EXUSUKx", "S&P 500")]
#' prior_sigma <- specify_prior_sigma(data=train_data, type="cholesky")
#' mod <- bvar(train_data, lags=2L, draws=2000, prior_sigma=prior_sigma) 
#' plot(irf(mod))
#'
#' @export
irf <- function(x, ahead=8, rotation=NULL, shocks=NULL) {
	n_variables <- ncol(x$PHI)
	n_posterior_draws <- dim(x$PHI)[3];
	
	if (is.null(shocks)) {
		shocks <- diag(
			if (x$sigma_type == "factor") ncol(x$facload)
			else if (x$sigma_type == "cholesky") ncol(x$Y)
			else stop("unknown model")
		)
	}
	n_shocks <- ncol(shocks)

	ret <- unlist(irf_cpp(
		x$PHI,
		x$facload,
		x$U,
		x$logvar[nrow(x$logvar),,], #most recent log volatility
		shocks=shocks,
		ahead,
		rotation
	))
	dim(ret) <- c(n_variables, n_shocks, 1+ahead, n_posterior_draws)
	dimnames(ret) <- list(colnames(x$Y), paste("shock",1:n_shocks), paste0("t=",0:ahead))
	class(ret) <- "bayesianVARs_irf"
	ret
}

#' Identifying restrictions for the structural parameters
#'
#' @export
#' @param restrictions_* A matrix. The entries have following meaning
#' \describe{
#'	\item{NA:}{an unrestricted entry.}
#'	\item{0:}{The entry at this position is restricted to be zero (i.e. an exclusion restriction).}
#'	\item{A positive number:}{The sign of this entry should be positive.}
#'	\item{A negative number:}{The sign of this entry should be negative.}
#' }
#'
#' @examples
#' train_data <- 100 * usmacro_growth[,c("GDPC1", "GPDICTPI", "GS1", "M2REAL", "CPIAUCSL")]
#' prior_sigma <- specify_prior_sigma(data=train_data, type="cholesky")
#' x <- bvar(train_data, lags=2L, draws=10000, prior_sigma=prior_sigma)
#' 
#' restrictions_B0 <- cbind(
#' 	c(1, 0, 0, 0, 0),
#' 	c(NA, 1, 0, 0, 0),
#' 	c(NA, 0, NA, NA, NA),
#' 	c(NA, NA, NA, NA, NA),
#' 	c(NA, NA, NA, NA, NA)
#' )
#' 
#' restrictions_IR_inf <- cbind(
#' 	c(NA, NA, NA, NA, NA),
#' 	c(NA, NA, NA, NA, NA),
#' 	c(0, NA, NA, NA, NA),
#' 	c(NA, NA, NA, NA, NA),
#' 	c(NA, NA, NA, NA, NA)
#' )
#' 
#' rotation <- find_rotation(x,
#' 	restrictions_B0 = restrictions_B0,
#' 	restrictions_long_run_ir = restrictions_IR_inf
#' )
#' plot(irf(x, rotation=rotation))
find_rotation <- function(
	x,
	restrictions_facload = NULL,
	restrictions_B0_inv_t = NULL,
	restrictions_B0 = NULL,
	restrictions_structural_coeff = NULL,
	restrictions_long_run_ir = NULL,
	tolerance = 0.0
) {
	if (x$sigma_type == "factor") {
		forbidden_restrictions <- c(
			!is.null(restrictions_B0_inv_t),
			!is.null(restrictions_B0),
			!is.null(restrictions_structural_coeff),
			!is.null(restrictions_long_run_ir)
		)
		if (any(forbidden_restrictions))
			stop("For factor models, only the factor loadings can be restricted")
	}
	else if (x$sigma_type == "cholesky") {
		if (!is.null(restrictions_facload))
			stop("Factor loadings cannot be restricted for a Cholesky model")
	}
	else {
		stop("Unkown model type")
	}
	
	parameter_transformations <- compute_parameter_transformations(
		x$PHI,
		x$facload,
		x$U,
		x$logvar[nrow(x$logvar),,], #most recent log volatility,
		include_facload = !is.null(restrictions_facload),
		include_B0_inv_t = !is.null(restrictions_B0_inv_t),
		include_B0 = !is.null(restrictions_B0),
		include_structural_coeff = !is.null(restrictions_structural_coeff),
		include_long_run_ir = !is.null(restrictions_long_run_ir)
	)
	# order must match the output of `compute_parameter_transformations`!
	restrictions <- list(
		"facload" = restrictions_facload,
		"B0_inv_t" = restrictions_B0_inv_t,
		"B0" = restrictions_B0,
	    "structural_coeff" = restrictions_structural_coeff,
	    "restrictions_long_run_ir" = restrictions_long_run_ir
	)
	
	if (all(lengths(parameter_transformations) == 0) || all(lengths(restrictions) == 0)) {
		stop("None of the restrictions apply to this model")
	}

	find_rotation_cpp(
		parameter_transformations = parameter_transformations,
		restriction_specs = restrictions[lengths(restrictions) > 0],
		tolerance = tolerance
	)
}
