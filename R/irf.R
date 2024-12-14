#' Impulse response functions
#'
#' Effect of a shock on the factors or variables over time
#'
#' @param x An object of type `bayesianVARs_bvar`.
#' @param shock A vector of shocks. If a factor model was used, then `shock` is
#'		applied to the factors, and its dimension must be the number of factors.
#'
#'		If the Cholesky model was used, then the dimension of `shock` is expected
#'		to be equals the number of equations. Note that the contemporaneous effects
#'		of a shock in the Cholesky model depend on the ordering of equations.
#'
#' @return Returns a `bayesianVARs_irf` object
#'
#' @examples
#' train_data <- 100 * usmacro_growth[,c("GDPC1", "PCECC96", "GPDIC1", "AWHMAN", "GDPCTPI", "CES2000000008x", "FEDFUNDS", "GS10", "EXUSUKx", "S&P 500")]
#' prior_sigma <- specify_prior_sigma(data=train_data, type="cholesky")
#' mod <- bvar(train_data, lags=2L, draws=2000, prior_sigma=prior_sigma)
#' ir <- irf(mod, shock=c(0,0,0,0,0,0,1,0,0,0))
#' plot(ir, n_col=2)
#'
#' @export
irf <- function(x, shock, ahead=8, rotation=NULL) {
	if (x$sigma_type == "factor") {
		number_of_factors <- dim(x$facload)[2]
		if (length(shock) != number_of_factors) {
			stop("the shock vector must have dimension equals to the number of factors")
		}
	}
	else if (x$sigma_type == "cholesky") {
		if (length(shock) != ncol(x$PHI)) {
			stop("the shock vector must have dimension equals to the number of variables")
		}
	}
	else {
		stop("unknown model")
	}
		
	ret <- irf_cpp(
		x$PHI,
		x$facload,
		x$U,
		x$logvar[nrow(x$logvar),,], #most recent log volatility
		shock,
		ahead,
		rotation
	)

	colnames(ret) <- colnames(x$Y)
	class(ret) <- "bayesianVARs_irf"
	ret
}

specify_zero_restrictions <- function(spec) {
	ret <- array(0, dim=c(nrow(spec), nrow(spec), ncol(spec)))
	for (j in seq_len(ncol(spec))) {
		ret[,,j] = diag(spec[,j] + 1)
	}
	ret[is.na(ret)] <- 0
	ret
}

find_rotation <- function(
	x,
	restrictions_B0_inv_t = NULL,
	restrictions_B0 = NULL,
	restrictions_structural_coeff = NULL,
	restrictions_long_run_ir = NULL
) {
	parameter_transformations <- compute_parameter_transformations(
		x$PHI,
		x$facload,
		x$U,
		x$logvar[nrow(x$logvar),,], #most recent log volatility
		include_B0_inv_t = !is.null(restrictions_B0_inv_t),
		include_B0 = !is.null(restrictions_B0),
		include_structural_coeff = !is.null(restrictions_structural_coeff),
		include_long_run_ir = !is.null(restrictions_long_run_ir)
	)
	# order must match the output of `compute_parameter_transformations`!
	restrictions <- list("B0_inv_t" = restrictions_B0_inv_t,
			 "B0" = restrictions_B0,
		     "structural_coeff" = restrictions_structural_coeff,
		     "restrictions_long_run_ir" = restrictions_long_run_ir)
	find_rotation_cpp(
		parameter_transformations,
		restrictions = restrictions[lengths(restrictions) > 0]
	)
}
