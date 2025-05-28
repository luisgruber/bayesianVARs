dim_shocks <- function(x) {
	if (x$sigma_type == "factor") ncol(x$facload)
	else if (x$sigma_type == "cholesky") ncol(x$Y)
	else stop("unknown model")
}

mspe_decomposition <- function(ir) {
	# cf. Structural Vector Autoregressive Analysis, Kilian & Lütkepohl, section 4.2

	# cumsum over time (-> dimension 3 of ir)
	# result has dimensions: time * variables * shocks * draws
	mspe_h <- apply(ir^2, MARGIN=c(1,2,4), FUN=cumsum)

	# proportions of the mspe_h over the shocks (-> dimension 3 of mspe_h)
	# result has dimensions: time * variables * shocks * draws
	ret <- proportions(mspe_h, margin=c(1,2,4))

	# use the same order of dimensions as `irf` for consistency
	aperm(ret, c(2, 3, 1, 4))
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
specify_structural_restrictions <- function(
	x,
	restrictions_facload = NULL,
	restrictions_B0_inv_t = NULL,
	restrictions_B0 = NULL,
	restrictions_structural_coeff = NULL,
	restrictions_long_run_ir = NULL
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
	# order must match the output of `compute_parameter_transformations`!
	list(
		"facload" = restrictions_facload,
		"B0_inv_t" = restrictions_B0_inv_t,
		"B0" = restrictions_B0,
	    "structural_coeff" = restrictions_structural_coeff,
	    "restrictions_long_run_ir" = restrictions_long_run_ir
	)
}

find_rotation <- function(
	x,
	structural_restrictions,
	t = nrow(x$logvar),
	
	solver = "randomized",
	randomized_max_rotations_per_sample = 2
) {
	stopifnot(length(t) == 1)
	stopifnot(solver %in% c("randomized", "lp"))
	stopifnot(randomized_max_rotations_per_sample > 0)
	if (all(lengths(structural_restrictions) == 0)) {
		stop("No restrictions specified")
	}
	structural_restrictions <- structural_restrictions[lengths(structural_restrictions) > 0]
    	
	parameter_transformations <- compute_parameter_transformations(
		x$PHI,
		x$facload,
		x$U,
		x$logvar[t,,],
		structural_restrictions
	)
	
	find_rotation_cpp(
		parameter_transformations = parameter_transformations,
		restriction_specs = structural_restrictions,
		solver_option = solver,
		randomized_max_rotations_per_sample = randomized_max_rotations_per_sample
	)
}

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
irf <- function(x, ahead=8, structural_restrictions=NULL, shocks=NULL, hairy=FALSE, t = nrow(x$logvar), ...) {
	n_variables <- ncol(x$PHI)
	n_posterior_draws <- dim(x$PHI)[3]
	
	if (is.null(shocks)) {
		shocks <- diag(dim_shocks(x))
	}
	n_shocks <- ncol(shocks)
	
	if (length(t) != 1 || t > nrow(x$logvar)) {
		stop("t must be a single valid timepoint. This needs sv_keep = \"all\" when calling bvar.")
	}
	
	rotation <- NULL
	rotation_sample_map <- 1:n_posterior_draws
	if (!is.null(structural_restrictions)) {
		result <- find_rotation(x, structural_restrictions, t=t, ...)
		rotation <- result$rotation
		rotation_sample_map <- result$rotation_sample_map
	}
	if (length(rotation_sample_map) == 0) {
		stop("The restrictions could not be satisfied")
	}
	if (length(rotation_sample_map) < x$config$draws * 0.25) {
		warning(paste("There are", length(rotation_sample_map), "samples remaining for IRF calculation"))
	}
	
	factor_loadings <- x$facload
	if (length(factor_loadings) > 0) {
		factor_loadings <- factor_loadings[,,rotation_sample_map]
	}
	
	U_vecs <- x$U
	if (length(U_vecs) > 0) {
		U_vecs <- U_vecs[,rotation_sample_map]
	}
	
	ret <- irf_cpp(
		coefficients=x$PHI[,,rotation_sample_map],
		factor_loadings=factor_loadings,
		U_vecs=U_vecs,
		logvar_t=x$logvar[t,,rotation_sample_map],
		shocks=shocks,
		ahead=ahead,
		rotation_=rotation
	)
	
	hair_order <- NULL
	if (hairy) {
		hair_order <- irf_bayes_optimal_order(ret) + 1 #R indices are one based!
		dim(hair_order) <- NULL
	}
	
	ret <- unlist(ret)
	dim(ret) <- c(n_variables, n_shocks, 1+ahead, length(rotation_sample_map))
	dimnames(ret) <- list(colnames(x$Y), paste("shock",1:n_shocks), paste0("t=",0:ahead))
	class(ret) <- "bayesianVARs_irf"
	attr(ret, "hair_order") <- hair_order
	ret
}

estimate_structural_shocks <- function(x, structural_restrictions, ...) {
	stop("not yet implemented")
	n_shocks <- dim_shocks(x)
	n_posterior_draws <- dim(x$PHI)[3]
	
	reduced_shocks <- residuals(x)
	
	structural_shocks_dim <- dim(reduced_shocks)
	structural_shocks_dim[2] <- n_shocks
	structural_shocks_dimnames <- dimnames(reduced_shocks)
	structural_shocks_dimnames[2] <- NULL
	structural_shocks <- array(dim=structural_shocks_dim, dimnames=structural_shocks_dimnames)
	for (t in seq_len(nrow(reduced_shocks))) {
		rot <- find_rotation(x, structural_restrictions, t=t, ...)
		B0 <- compute_parameter_transformations(
			x$PHI,
			x$facload,
			x$U,
			x$logvar[t,,],
			list("B0" = TRUE)
		)[["B0"]]
		for (r in seq_len(n_posterior_draws)) {
			structural_shocks[t,,r] <- reduced_shocks[t,,r] %*% B0[,,r] %*% rot[,,r]
		}
	}
	
	structural_shocks
}
