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

# obtain_restrictable_matrices(
#		x$PHI,
#		x$facload,
#		x$U,
#		x$logvar[nrow(x$logvar),,] #most recent log volatility
# )
