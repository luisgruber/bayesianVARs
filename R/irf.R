#' Impulse response functions
#'
#' Effect of a shock on the variables over time
#'
#' @param x An object of type `bayesianVARs_bvar`
#'
#' @export
irf <- function(x, factor_shock, ahead=8) {
	if (x$sigma_type != "factor") {
		stop("impulse reponse functions are only available for factor models")
	}
	
	number_of_factors <- dim(x$facload)[2]
	if (length(factor_shock) != number_of_factors) {
		stop("the shock vector to the factors must have dimension equals to the number of factors")
	}

	ret <- irf_cpp(
		x$PHI,
		x$facload,
		factor_shock,
		ahead
	)

	colnames(ret) <- colnames(x$Y)
	class(ret) <- "bayesianVARs_irf"
	ret
}
