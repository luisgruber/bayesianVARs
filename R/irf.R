#' Impulse response functions
#'
#' Effect of a shock on the variables over time
#'
#' @param x An object of type `bayesianVARs_bvar`
#'
#' @export
irf <- function(x, shock, ahead=8) {
	if (x$sigma_type == "factor") {
		number_of_factors <- dim(x$facload)[2]
		if (length(shock) != number_of_factors) {
			stop("the shock vector must have dimension equals to the number of factors")
		}
	}
	else if (x$sigma_type == "cholesky") {
		if (length(shock) != ncol(mod$PHI)) {
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
		shock,
		ahead
	)

	colnames(ret) <- colnames(x$Y)
	class(ret) <- "bayesianVARs_irf"
	ret
}
