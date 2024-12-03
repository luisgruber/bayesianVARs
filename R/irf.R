# x should be of class bayesianVARs_bvar
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
	dimnames(ret)[[2]] <- dimnames(x$Y)[[2]]
	ret
}
