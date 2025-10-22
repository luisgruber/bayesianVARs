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

#' Set identifying restrictions for the structural VAR parameters.
#'
#' @param x An object of type `bayesianVARs_bvar`.
#'
#' @param restrictions_facload an M times r matrix of restrictions on the factor loadings.
#' This is equivalent to the instantaneous effects of the factor shocks.
#' Can only be used, if the factor decomposition for sigma was specified.
#'
#' @param restrictions_B0_inv_t an M times M matrix of restrictions on the instantaneous
#' effects of the structural shocks. Columns correspond to shocks, rows to variables.
#' Can only be used, if the cholesky decomposition for sigma was specified.
#'
#' @param restrictions_B0 an M times M matrix of restrictions on the structural matrix.
#' Can only be used, if the cholesky decomposition for sigma was specified.
#'
#' @param restrictions_structural_coeff a matrix of restrictions on the structural VAR coefficients.
#' Its size should match the dimensions of \code{x$PHI}.
#' Can only be used, if the cholesky decomposition for sigma was specified.
#'
#' @param restrictions_long_run_ir a matrix of restrictions on the long run impulse responses.
#' The long run impulse responses are the sum of the impulse responses summed over all time horizons.
#' Restrictions on the long run impulse responses can be specified for both the factor and the
#' cholesky decomposition of sigma. In the case of a factor decomposition its size is expected to be
#' M times r. In the case of a cholesky decomposition the size must be M times M.
#'
#' @details All \code{restrictions_*} entries have following meaning
#' \describe{
#'	\item{NA:}{an unrestricted entry.}
#'	\item{0:}{The entry at this position is restricted to be zero (i.e. an exclusion restriction).}
#'	\item{A positive number:}{The sign of this entry should be positive.}
#'	\item{A negative number:}{The sign of this entry should be negative.}
#' }
#'
#' The structural VAR(p) model is of the following form: \deqn{\boldsymbol{y}^\prime_t \boldsymbol{B}_0
#' = \boldsymbol{x}^\prime_t\boldsymbol{\Phi} \boldsymbol{B}_0 + \boldsymbol{\omega}^\prime_t}
#'
#' @examples
#'train_data <- 100 * usmacro_growth[,c("GDPC1", "GDPCTPI", "GS1", "M2REAL", "CPIAUCSL")]
#'prior_sigma <- specify_prior_sigma(M=ncol(train_data), type="cholesky", cholesky_heteroscedastic=FALSE)
#'mod <- bvar(train_data, lags=5L, draws=3000, prior_sigma=prior_sigma)
#'
#'structural_restrictions <- specify_structural_restrictions(
#'  mod,
#'  restrictions_B0=rbind(
#'    c(1 ,NA,0 ,NA,NA),
#'    c(0 ,1 ,0 ,NA,NA),
#'    c(0 ,NA,1 ,NA,NA),
#'    c(0 ,0 ,NA,1 ,NA),
#'    c(0 ,0 ,0 ,0 ,1 )
#'  )
#')
#'irf_structural <- irf(
#'  mod, ahead=8,
#'  structural_restrictions=structural_restrictions
#')
#'plot(irf_structural)
#'
#' @seealso [`irf`], [`extract_B0`], [`specify_prior_sigma`]
#' @author Stefan Haan \email{sthaan@edu.aau.at}
#' @export
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
			!is.null(restrictions_structural_coeff)
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
#' Effect of the structural (factor) shocks over time.
#' 
#' If a factor model was used, then the number of shocks is equal to the number
#' of factors.
#'
#' If the Cholesky model was used, then the number of shocks is equal to the
#' number of variables.
#'
#' @param x An object of type `bayesianVARs_bvar`.
#' @param ahead maximum number of time horizons.
#' @param structural_restrictions an object generated by [`specify_structural_restrictions`]. If specified,
#' the IRFs to the structural shocks identified by the restrictions given in this argument will be calculated.
#' If not specified, the orthogonal IRFs based on the Cholesky decomposition will be calculated.
#' Note that the orthogonal IRFs depend on the ordering of the variables and do not necessarily
#' correspond to shocks with a meaningful interpretation.
#'
#' @param shocks an matrix with r rows, where r is the number of shocks (see 'Details'). Each column specifies a shock.
#' Default: \code{diag(r)}, will calculate the responses to all structural (or factor) shocks of one standard deviation.
#' @param hairy set to \code{TRUE} in order to plot each path seperately.
#' To show valid quantiles, an Bayes optimal order of the posterior samples will be calculated which can take a long time
#' even for moderately many samples. Default: \code{FALSE}.
#' @param ... Following expert arguments can be specified:
#' \describe{
#'  \item{solver:}{\code{"randomized"} or \code{"lp"}. If some columns have more than one sign restriction, \code{"lp"}
#'  might find a solution, even if \code{"randomized"} is unable to. However \code{"lp"} can produce artifically narrow confidence bands
#'  which do not properly reflect the uncertainty in the identification scheme. Default: \code{"randomized"}}
#'  \item{randomized_max_rotations_per_sample:}{if using the \code{"randomized"} solver, how many rotations are drawn for each sample of
#'  the reduced form parameters in \code{x}. Default: \code{2}.}
#' }
#'
#' @return Returns a `bayesianVARs_irf` object.
#'
#' @examples
#'train_data <- 100 * usmacro_growth[,c("GDPC1", "GDPCTPI", "GS1", "M2REAL", "CPIAUCSL")]
#'prior_sigma <- specify_prior_sigma(M=ncol(train_data), type="cholesky", cholesky_heteroscedastic=FALSE)
#'mod <- bvar(train_data, lags=5L, draws=3000, prior_sigma=prior_sigma)
#'
#'structural_restrictions <- specify_structural_restrictions(
#'  mod,
#'  restrictions_B0=rbind(
#'    c(1 ,NA,0 ,NA,NA),
#'    c(0 ,1 ,0 ,NA,NA),
#'    c(0 ,NA,1 ,NA,NA),
#'    c(0 ,0 ,NA,1 ,NA),
#'    c(0 ,0 ,0 ,0 ,1 )
#'  )
#')
#'irf_structural <- irf(
#'  mod, ahead=8,
#'  structural_restrictions=structural_restrictions
#')
#'plot(irf_structural)
#'
#' @export
#' @seealso [`specify_structural_restrictions`], [`extract_B0`]
#' @references Arias, J. and Rubio-Ramírez, J. and Waggoner, D. (2014).
#'  Inference Based on SVARs Identified with Sign and Zero Restrictions: Theory and Applications.
#'  \emph{FRB Atlanta Working Paper Series}, \doi{10.2139/ssrn.2580264}.
#'
#' @author Stefan Haan \email{sthaan@edu.aau.at}
irf <- function(x, ahead=8, structural_restrictions=NULL, shocks=NULL, hairy=FALSE, ...) {
	t = nrow(x$logvar)
	n_variables <- ncol(x$PHI)
	n_posterior_draws <- dim(x$PHI)[3]
	
	r <- dim_shocks(x)
	if (is.null(shocks)) {
		shocks <- diag(r)
	}
	else {
	  if (nrow(shocks) != r) {
	    stop(paste("'shocks' must be a matrix with", r, "rows."))
	  }
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
		factor_loadings <- factor_loadings[,,rotation_sample_map, drop=FALSE]
	}
	
	U_vecs <- x$U
	if (length(U_vecs) > 0) {
		U_vecs <- U_vecs[,rotation_sample_map]
	}
	
	ret <- irf_cpp(
		coefficients=x$PHI[,,rotation_sample_map, drop=FALSE],
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
	dimnames(ret) <- list(
		variables=colnames(x$Y),
		shocks=colnames(shocks),
		horizon=paste0("t=",0:ahead),
		samples=NULL
	)
	class(ret) <- "bayesianVARs_irf"
	attr(ret, "hair_order") <- hair_order
	ret
}

#' @export
print.bayesianVARs_irf <- function(x) {
	d <- dim(x)
	cat(
		"Impulse responses of",
		d[1], "variables to",
		d[2], "shocks over",
		d[3], "periods. \n"
	)
}


#' Retrieve the structural parameter \eqn{\boldsymbol{B}_0} samples from an IRF object.
#'
#' @param x a `bayesianVARs_irf` object
#' @export
#' @seealso [`specify_structural_restrictions`]
#' @author Stefan Haan \email{sthaan@edu.aau.at}
extract_B0 <- function(x) {
	if (ncol(x) != nrow(x)) {
		stop("IRFs must have an equal number of variables and shocks")
	}
	n_posterior_draws <- dim(x)[4]
	B0_samples <- array(
		NA,
		dim=c(nrow(x), ncol(x), n_posterior_draws),
		dimnames=list(variables=rownames(x), shocks=colnames(x), samples=NULL)
	)
	for (r in seq_len(n_posterior_draws)) {
		B0_samples[,,r] <- t(solve(x[,,1,r]))
	}
	B0_samples
}
