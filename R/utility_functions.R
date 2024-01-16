
# Methods for bayesianVARs_bvar objects -----------------------------------

#' Pretty printing of a bvar object
#'
#' @param x Object of class `bayesianVARs_bvar`, usually resulting from a call
#' of [`bvar()`].
#' @param ... Ignored.
#'
#' @return Returns \code{x} invisibly.
#'
#' @export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Estimate a model
#' mod <- bvar(data, sv_keep = "all", quiet = TRUE)
#'
#' # Print model
#' mod
#'
print.bayesianVARs_bvar <- function(x, ...) {
  M <- ncol(x$Y)
  Tobs <- nrow(x$Y)
  Traw <- nrow(x$Yraw)

  if(inherits(x, "bayesianVARs_bvar_stable")){
    cat(paste("\nFitted bayesianVARs_bvar_stable object with\n",
              " -", formatC(dim(x$PHI)[3], width = 7), "remaining stable draws\n"))
  }
  cat(paste("\nFitted bayesianVARs_bvar object with\n",
            " -", formatC(M, width = 7), "series\n",
            " -", formatC(x$lags, width = 7), "lag(s)\n",
            " -", formatC(Tobs, width = 7), "used observations\n",
            " -", formatC(Traw, width = 7), "total observations\n",
            " -", formatC(x$config$draws, width = 7), "MCMC draws\n",
            " -", formatC(x$config$thin, width = 7), "thinning\n",
            " -", formatC(x$config$burnin, width = 7), "burn-in\n\n"))
  invisible(x)
}

#' Summary method for bayesianVARs_bvar objects
#'
#' Summary method for `bayesianVARs_bvar` objects.
#'
#' @param object A `bayesianVARs_bvar` object obtained via [`bvar()`].
#' @param quantiles numeric vector which quantiles to compute.
#' @param ... Currently ignored!
#'
#' @return An object of type `summary.bayesianVARs_bvar`.
#' @export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Estimate model
#' mod <- bvar(data, quiet = TRUE)
#'
#' # Summary
#' sum <- summary(mod)
summary.bayesianVARs_bvar <- function(object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),...){
  PHImedian <- apply(object$PHI, 1:2, stats::median)
  PHIquantiles <- apply(object$PHI, 1:2, stats::quantile, quantiles)
  PHIiqr <- apply(object$PHI, 1:2, stats::IQR)
  if(object$sigma_type == "cholesky"){
    Umedian <- Uiqr <- diag(1, nrow = ncol(object$Y))
    colnames(Umedian) <- rownames(Umedian) <-
      colnames(Uiqr) <- rownames(Uiqr) <- colnames(object$Y)
    Umedian[upper.tri(Umedian)] <- apply(object$U, 1, stats::median)
    Uquantiles <- apply(object$U, 1, stats::quantile, quantiles)
    Uiqr[upper.tri(Uiqr)] <- apply(object$U, 1, stats::IQR)
  }

  out <- list(PHImedian = PHImedian,
              PHIquantiles = PHIquantiles,
              PHIiqr = PHIiqr
  )
  if(object$sigma_type == "cholesky"){
    out$Umedian <- Umedian
    out$Uquantiles <- Uquantiles
    out$Uiqr <- Uiqr
  }
  out$sigma_type <- object$sigma_type

  class(out) <- "summary.bayesianVARs_bvar"
  out
}

#' Print method for summary.bayesianVARs_bvar objects
#'
#' Print method for `summary.bayesianVARs_bvar` objects.
#'
#' @param x A `summary.bayesianVARs_bvar` object obtained via
#'   [`summary.bayesianVARs_bvar()`].
#' @param ... Currently ignored!
#'
#' @return Returns `x` invisibly!
#' @export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Estimate model
#' mod <- bvar(data, quiet = TRUE)
#'
#' # Print summary
#' summary(mod)
print.summary.bayesianVARs_bvar <- function(x, ...){
  digits <- max(3, getOption("digits") - 3)
  cat("\nPosterior median of reduced-form coefficients:\n")
  print(x$PHImedian, digits = digits)
  cat("\nPosterior interquartile range of of reduced-form coefficients:\n")
  print(x$PHIiqr, digits = digits)
  if(x$sigma_type == "cholesky"){
    cat("\nPosterior median of contemporaneous coefficients:\n")
    print(as.table(x$Umedian - diag(nrow(x$Umedian))), digits = digits, zero.print = "-")
    cat("\nPosterior interquartile range of contemporaneous coefficients:\n")
    print(as.table(x$Uiqr- diag(nrow(x$Uiqr))), digits = digits, zero.print = "-")
  }
  invisible(x)
}


#' @name coef
#' @title Extract VAR coefficients
#' @description Extracts posterior draws of the VAR coefficients from a VAR
#'   model estimated with [`bvar()`].
#'
#' @param object A `bayesianVARs_bvar` object obtained from [`bvar()`].
#' @param ... Currently ignored.
#'
#' @return Returns a numeric array of dimension \eqn{M \times K \times draws},
#'   where M is the number of time-series, K is the number of covariates per
#'   equation (including the intercept) and draws is the number of stored
#'   posterior draws.
#' @seealso [`summary.bayesianVARs_draws()`], [`vcov.bayesianVARs_bvar()`].
#' @export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Estimate a model
#' mod <- bvar(data, sv_keep = "all", quiet = TRUE)
#'
#' # Extract posterior draws of VAR coefficients
#' bvar_coefs <- coef(mod)
coef.bayesianVARs_bvar <- function(object, ...){
  ret <- object$PHI
  class(ret) <- c("bayesianVARs_coef", "bayesianVARs_draws")
  ret
}

#' Extract posterior draws of the (time-varying) variance-covariance matrix for
#' a VAR model
#'
#' Returns the posterior draws of the possibly time-varying variance-covariance
#' matrix of a VAR estimated via [`bvar()`]. Returns the full paths if
#' `sv_keep="all"` when calling [`bvar()`]. Otherwise, the draws of the
#' variance-covariance matrix for the last observation are returned, only.
#'
#' @param object An object of class `bayesianVARs_bvar` obtained via [`bvar()`].
#' @param ... Currently ignored.
#'
#' @return An array of class `bayesianVARs_draws` of dimension \eqn{T \times M
#'   \times M \times draws}, where \eqn{T} is the number of observations,
#'   \eqn{M} the number of time-series and \eqn{draws} the number of stored
#'   posterior draws.
#' @seealso [`summary.bayesianVARs_draws`], [`coef.bayesianVARs_bvar()`].
#' @export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#' # Estimate a model
#' mod <- bvar(data, sv_keep = "all", quiet = TRUE)
#'
#' # Extract posterior draws of the variance-covariance matrix
#' bvar_vcov <- vcov(mod)
vcov.bayesianVARs_bvar <- function(object, ...){
  M <- ncol(object$Y)
  factors <- dim(object$facload)[2]
  Tobs <- nrow(object$logvar)
  nsave <- dim(object$PHI)[3]

  dates <- NULL
  if(Tobs==1L){
    dates <- nrow(object$Y)
  }else if(Tobs == nrow(object$Y)){
    dates <- 1:nrow(object$Y)
  }

  if(!is.null(rownames(object$Y)) & !is.null(dates)){
    dates <- tryCatch(as.Date(rownames(object$Yraw)[dates]),
                      error = function(e) dates)
  }

  out <- vcov_cpp(object$sigma_type == "factor",
                  object$facload,
                  object$logvar,
                  object$U,
                  M,
                  factors)
  out <- array(out, c(Tobs, M, M, nsave))
  dimnames(out) <- list(as.character(dates), colnames(object$Y),
                        colnames(object$Y), NULL)
  class(out) <- "bayesianVARs_draws"
  out
}

#' Summary statistics for bayesianVARs posterior draws.
#'
#' @param object An object of class `bayesianVARs_draws` usually obtained through
#'   extractors like [`coef.bayesianVARs_bvar()`] and
#'   [`vcov.bayesianVARs_bvar()`].
#' @param quantiles A vector of quantiles to evaluate.
#' @param ... Currently ignored.
#'
#' @return A list object of class `bayesianVARs_draws_summary` holding
#' * `mean`: Vector or matrix containing the posterior mean.
#' * `sd`: Vector or matrix containing the posterior standard deviation .
#' * `quantiles`: Array containing the posterior quantiles.
#' @seealso Available extractors: [`coef.bayesianVARs_bvar()`],
#'   [`vcov.bayesianVARs_bvar()`].
#' @export
#'
#' @examples
#'
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#' # Estimate a model
#' mod <- bvar(data, sv_keep = "all", quiet = TRUE)
#'
#' # Extract posterior draws of VAR coefficients
#' bvar_coefs <- coef(mod)
#'
#' # Compute summary statistics
#' summary_stats <- summary(bvar_coefs)
#'
#' # Compute summary statistics of VAR coefficients without using coef()
#' summary_stats <- summary(mod$PHI)
#'
#' # Test which list elements of 'mod' are of class 'bayesianVARs_draws'.
#' names(mod)[sapply(names(mod), function(x) inherits(mod[[x]], "bayesianVARs_draws"))]
#'
summary.bayesianVARs_draws <- function(object, quantiles = c(0.25, 0.5, 0.75),
                                       ...){

  dim_object <- dim(object)
  dim_length <- length(dim_object)
  object_mean <- apply(object, 1:(dim_length-1), mean)
  object_sd <- apply(object, 1:(dim_length-1), sd)
  object_quantiles <- apply(object, 1:(dim_length-1), quantile, quantiles)

  out <- list(
    mean = object_mean,
    sd = object_sd,
    quantiles = object_quantiles
  )
  class(out) <- "bayesianVARs_draws_summary"
  return(out)
}

#' Simulate fitted/predicted historical values for an estimated VAR model
#'
#' Simulates the fitted/predicted (in-sample) values for an estimated VAR model.
#'
#' @param object A `bayesianVARs_bvar` object estimated via [bvar()].
#' @param error_term logical indicating whether to include the error term or
#'   not.
#' @param ... Currently ignored.
#'
#' @return An object of class `bayesianVARs_fitted`.
#' @export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Estimate a model
#' mod <- bvar(data, sv_keep = "all", quiet = TRUE)
#'
#' # Simulate predicted historical values including the error term.
#' pred <- fitted(mod, error_term = TRUE)
#'
#' # Simulate fitted historical values not including the error term.
#' fit <- fitted(mod, error_term = FALSE)
#'
#' # Visualize
#' plot(pred)
#' plot(fit)
fitted.bayesianVARs_bvar <- function(object, error_term = TRUE, ...){

  if(dim(object$logvar)[1] != nrow(object$Y) & error_term){
    warning("Setting 'error_term=FALSE'! To calculate predicted historical
    values including the error term, the full path of logvariances is needed,
    i.e. set 'sv_keep='all'' when calling bvar()!")
    error_term <- FALSE
  }

  ret <- insample(object$X,
                  object$PHI,
                  object$U,
                  object$facload,
                  object$logvar,
                  error_term,
                  object$sigma_type == "factor")
  colnames(ret) <- colnames(object$Y)
  rownames(ret) <- rownames(object$Y)
  out <- list(
    fitted = ret,
    Yraw = object$Yraw
  )
  class(out) <- "bayesianVARs_fitted"
  out
}

# Functions for prior configuration ---------------------------------------

#'Specify prior on Sigma
#'
#'Configures prior on the variance-covariance of the VAR.
#'
#'\code{bvar} offers two different specifications for the errors: The user can
#'choose between a factor stochastic volatility structure or a cholesky
#'stochastic volatility structure. In both cases the disturbances
#'\eqn{\boldsymbol{\epsilon}_t} are assumed to follow a \eqn{M}-dimensional
#'multivariate normal distribution with zero mean and variance-covariance matrix
#'\eqn{\boldsymbol{\Sigma}_t}. In case of the cholesky specification
#'\eqn{\boldsymbol{\Sigma}_t = \boldsymbol{U}^{\prime -1} \boldsymbol{D}_t
#'\boldsymbol{U}^{-1}}, where \eqn{\boldsymbol{U}^{-1}} is upper unitriangular (with ones on
#'the diagonal). The diagonal matrix \eqn{\boldsymbol{D}_t} depends upon latent
#'log-variances, i.e. \eqn{\boldsymbol{D}_t=diag(exp(h_{1t}),\dots, exp(h_{Mt})}. The
#'log-variances follow a priori independent autoregressive processes
#'\eqn{h_{it}\sim N(\mu_i + \phi_i(h_{i,t-1}-\mu_i),\sigma_i^2)} for
#'\eqn{i=1,\dots,M}. In case of the factor structure,
#'\eqn{\boldsymbol{\Sigma}_t = \boldsymbol{\Lambda} \boldsymbol{V}_t \boldsymbol{\Lambda}^\prime +
#'\boldsymbol{G}_t}. The diagonal matrices \eqn{\boldsymbol{V}_t} and
#'\eqn{\boldsymbol{G}_t} depend upon latent log-variances, i.e.
#'\eqn{\boldsymbol{G}_t=diag(exp(h_{1t}),\dots, exp(h_{Mt})} and
#'\eqn{\boldsymbol{V}_t=diag(exp(h_{M+1,t}),\dots, exp(h_{M+r,t})}. The log-variances
#'follow a priori independent autoregressive processes \eqn{h_{it}\sim N(\mu_i +
#'\phi_i(h_{i,t-1}-\mu_i),\sigma_i^2)} for \eqn{i=1,\dots,M} and
#'\eqn{h_{M+j,t}\sim N(\phi_ih_{M+j,t-1},\sigma_{M+j}^2)} for \eqn{j=1,\dots,r}.
#'
#'@param data Optional. Data matrix (can be a time series object). Each of
#'  \eqn{M} columns is assumed to contain a single time-series of length
#'  \eqn{T}.
#'@param M positive integer indicating the number of time-series of the VAR.
#'@param type character, one of `"factor"` (the default) or `"cholesky"`,
#'  indicating which decomposition to be applied to the covariance-matrix.
#'@param factor_factors Number of latent factors to be estimated. Only required
#'  if `type="factor"`.
#'@param factor_restrict Either "upper" or "none", indicating whether the factor
#'  loadings matrix should be restricted to have zeros above the diagonal
#'  ("upper") or whether all elements should be estimated from the data
#'  ("none"). Setting \code{restrict} to "upper" often stabilizes MCMC
#'  estimation and can be important for identifying the factor loadings matrix,
#'  however, it generally is a strong prior assumption. Setting \code{restrict}
#'  to "none" is usually the preferred option if identification of the factor
#'  loadings matrix is of less concern but covariance estimation or prediction
#'  is the goal. Only required if `type="factor"`.
#'@param factor_priorfacloadtype Can be \code{"normal"}, \code{"rowwiseng"},
#'  \code{"colwiseng"}. Only required if `type="factor"`.
#' \describe{
#'  \item{\code{"normal"}: }{Normal prior. The value of \code{priorfacload}
#'                           is interpreted as the standard deviations of the
#'                           Gaussian prior distributions for the factor loadings.}
#'  \item{\code{"rowwiseng"}: }{Row-wise Normal-Gamma prior. The value of \code{priorfacload}
#'                              is interpreted as the shrinkage parameter \code{a}.}
#'  \item{\code{"colwiseng"}: }{Column-wise Normal-Gamma prior. The value of \code{priorfacload}
#'                              is interpreted as the shrinkage parameter \code{a}.}
#' }
#'  For details please see Kastner (2019).
#'@param factor_priorfacload Either a matrix of dimensions \code{M} times
#'  \code{factor_factors} with positive elements or a single number (which will
#'  be recycled accordingly). Only required if `type="factor"`. The meaning of
#'  \code{factor_priorfacload} depends on the setting of
#'  \code{factor_priorfacloadtype} and is explained there.
#'@param factor_facloadtol Minimum number that the absolute value of a factor
#'  loadings draw can take. Prevents numerical issues that can appear when
#'  strong shrinkage is enforced if chosen to be greater than zero. Only
#'  required if `type="factor"`.
#'@param factor_priorng Two-element vector with positive entries indicating the
#'  Normal-Gamma prior's hyperhyperparameters \code{c} and \code{d} (cf. Kastner
#'  (2019)). Only required if `type="factor"`.
#'@param factor_priormu Vector of length 2 denoting prior mean and standard
#'  deviation for unconditional levels of the idiosyncratic log variance
#'  processes. Only required if `type="factor"`.
#'@param factor_priorphiidi Vector of length 2, indicating the shape parameters
#'  for the Beta prior distributions of the transformed parameters
#'  \code{(phi+1)/2}, where \code{phi} denotes the persistence of the
#'  idiosyncratic log variances. Only required if `type="factor"`.
#'@param factor_priorphifac Vector of length 2, indicating the shape parameters
#'  for the Beta prior distributions of the transformed parameters
#'  \code{(phi+1)/2}, where \code{phi} denotes the persistence of the factor log
#'  variances. Only required if `type="factor"`.
#'@param factor_priorsigmaidi Vector of length \code{M} containing the prior
#'  volatilities of log variances. If \code{factor_priorsigmaidi} has exactly
#'  one element, it will be recycled for all idiosyncratic log variances. Only
#'  required if `type="factor"`.
#'@param factor_priorsigmafac Vector of length \code{factor_factors} containing
#'  the prior volatilities of log variances. If \code{factor_priorsigmafac} has
#'  exactly one element, it will be recycled for all factor log variances. Only
#'  required if `type="factor"`.
#'@param factor_priorh0idi Vector of length 1 or \code{M}, containing
#'  information about the Gaussian prior for the initial idiosyncratic log
#'  variances. Only required if `type="factor"`. If an element of
#'  \code{factor_priorh0idi} is a nonnegative number, the conditional prior of
#'  the corresponding initial log variance h0 is assumed to be Gaussian with
#'  mean 0 and standard deviation \code{factor_priorh0idi} times \eqn{sigma}. If
#'  an element of \code{factor_priorh0idi} is the string 'stationary', the prior
#'  of the corresponding initial log volatility is taken to be from the
#'  stationary distribution, i.e. h0 is assumed to be Gaussian with mean 0 and
#'  variance \eqn{sigma^2/(1-phi^2)}.
#'@param factor_priorh0fac Vector of length 1 or \code{factor_factors},
#'  containing information about the Gaussian prior for the initial factor log
#'  variances. Only required if `type="factor"`. If an element of
#'  \code{factor_priorh0fac} is a nonnegative number, the conditional prior of
#'  the corresponding initial log variance h0 is assumed to be Gaussian with
#'  mean 0 and standard deviation \code{factor_priorh0fac} times \eqn{sigma}. If
#'  an element of \code{factor_priorh0fac} is the string 'stationary', the prior
#'  of the corresponding initial log volatility is taken to be from the
#'  stationary distribution, i.e. h0 is assumed to be Gaussian with mean 0 and
#'  variance \eqn{sigma^2/(1-phi^2)}.
#'@param factor_heteroskedastic Vector of length 1, 2, or \code{M +
#'  factor_factors}, containing logical values indicating whether time-varying
#'   (\code{factor_heteroskedastic = TRUE}) or constant (\code{factor_heteroskedastic =
#'   FALSE}) variance should be estimated. If \code{factor_heteroskedastic} is of
#'  length 2 it will be recycled accordingly, whereby the first element is used
#'  for all idiosyncratic variances and the second element is used for all
#'  factor variances. Only required if `type="factor"`.
#'@param factor_priorhomoskedastic Only used if at least one element of
#'  \code{factor_heteroskedastic} is set to \code{FALSE}. In that case,
#'  \code{factor_priorhomoskedastic} must be a matrix with positive entries and
#'  dimension c(M, 2). Values in column 1 will be interpreted as the shape and
#'  values in column 2 will be interpreted as the rate parameter of the
#'  corresponding inverse gamma prior distribution of the idiosyncratic
#'  variances. Only required if `type="factor"`.
#'@param factor_interweaving The following values for interweaving the factor
#'  loadings are accepted (Only required if `type="factor"`):
#' \describe{
#'  \item{0: }{No interweaving.}
#'  \item{1: }{Shallow interweaving through the diagonal entries.}
#'  \item{2: }{Deep interweaving through the diagonal entries.}
#'  \item{3: }{Shallow interweaving through the largest absolute entries in each column.}
#'  \item{4: }{Deep interweaving through the largest absolute entries in each column.}
#' }
#'  For details please see Kastner et al. (2017). A value of 4 is the highly
#'  recommended default.
#'@param cholesky_U_prior character, one of \code{"HS"}, \code{"R2D2"}, `"NG"`,
#'  \code{"DL"}, \code{"SSVS"}, \code{"HMP"} or \code{"normal"}. Only required
#'  if `type="cholesky"`.
#'@param cholesky_U_tol Minimum number that the absolute value of an free
#'  off-diagonal element of an \eqn{U}-draw can take. Prevents numerical issues
#'  that can appear when strong shrinkage is enforced if chosen to be greater
#'  than zero. Only required if `type="cholesky"`.
#'@param cholesky_heteroscedastic single logical indicating whether time-varying
#'  (\code{cholesky_heteroscedastic = TRUE}) or constant
#'  (\code{cholesky_heteroscedastic = FALSE}) variance should be estimated. Only
#'  required if `type="cholesky"`.
#'@param cholesky_priormu Vector of length 2 denoting prior mean and standard
#'  deviation for unconditional levels of the log variance processes. Only
#'  required if `type="cholesky"`.
#'@param cholesky_priorphi Vector of length 2, indicating the shape parameters
#'  for the Beta prior distributions of the transformed parameters
#'  \code{(phi+1)/2}, where \code{phi} denotes the persistence of the log
#'  variances. Only required if `type="cholesky"`.
#'@param cholesky_priorsigma2 Vector of length 2, indicating the shape and the
#'  rate for the Gamma prior distributions on the variance of the log variance
#'  processes. (Currently only one global setting for all \eqn{M} processes is
#'  supported). Only required if `type="cholesky"`.
#'@param cholesky_priorh0 Vector of length 1 or \code{M}, containing information
#'  about the Gaussian prior for the initial idiosyncratic log variances. Only
#'  required if `type="cholesky"`. If an element of \code{cholesky_priorh0} is a
#'  nonnegative number, the conditional prior of the corresponding initial log
#'  variance h0 is assumed to be Gaussian with mean 0 and standard deviation
#'  \code{cholesky_priorh0} times \eqn{sigma}. If an element of
#'  \code{cholesky_priorh0} is the string 'stationary', the prior of the
#'  corresponding initial log volatility is taken to be from the stationary
#'  distribution, i.e. h0 is assumed to be Gaussian with mean 0 and variance
#'  \eqn{sigma^2/(1-phi^2)}.
#'@param cholesky_priorhomoscedastic Only used if
#'  \code{cholesky_heteroscedastic=FALSE}. In that case,
#'  \code{cholesky_priorhomoscedastic} must be a matrix with positive entries
#'  and dimension c(M, 2). Values in column 1 will be interpreted as the shape
#'  and values in column 2 will be interpreted as the scale parameter of the
#'  corresponding inverse gamma prior distribution of the variances. Only
#'  required if `type="cholesky"`.
#'@param cholesky_DL_a (Single) positive real number. The value is interpreted
#'  as the concentration parameter for the local scales. Smaller values enforce
#'  heavier shrinkage. A matrix of dimension `c(s,2)` specifies a discrete
#'  hyperprior, where the first column contains s support points and the second
#'  column contains the associated prior probabilities. `cholesky_DL_a` has only
#'  to be specified if `cholesky_U_prior="DL"`.
#'@param cholesky_DL_tol Minimum number that a parameter draw of one of the
#'  shrinking parameters of the Dirichlet Laplace prior can take. Prevents
#'  numerical issues that can appear when strong shrinkage is enforced if chosen
#'  to be greater than zero. `DL_tol` has only to be specified if
#'  `cholesky_U_prior="DL"`.
#'@param cholesky_R2D2_a (Single) positive real number. The value is interpreted
#'  as the concentration parameter for the local scales. Smaller values enforce
#'  heavier shrinkage. A matrix of dimension `c(s,2)` specifies a discrete
#'  hyperprior, where the first column contains s support points and the second
#'  column contains the associated prior probabilities. cholesky_R2D2_a has only
#'  to be specified if `cholesky_U_prior="R2D2"`.
#'@param cholesky_R2D2_b single positive number, where greater values indicate
#'  heavier regularization. \code{cholesky_R2D2_b} has only to be specified if
#'  \code{cholesky_U_prior="R2D2"}.
#'@param cholesky_R2D2_tol Minimum number that a parameter draw of one of the
#'  shrinking parameters of the R2D2 prior can take. Prevents numerical issues
#'  that can appear when strong shrinkage is enforced if chosen to be greater
#'  than zero. `cholesky_R2D2_tol` has only to be specified if
#'  `cholesky_U_prior="R2D2"`.
#'@param cholesky_NG_a (Single) positive real number. The value is interpreted
#'  as the concentration parameter for the local scales. Smaller values enforce
#'  heavier shrinkage. A matrix of dimension `c(s,2)` specifies a discrete
#'  hyperprior, where the first column contains s support points and the second
#'  column contains the associated prior probabilities. `cholesky_NG_a` has only
#'  to be specified if `cholesky_U_prior="NG"`.
#'@param cholesky_NG_b (Single) positive real number. The value indicates the
#'  shape parameter of the inverse gamma prior on the global scales.
#'  `cholesky_NG_b` has only to be specified if `cholesky_U_prior="NG"`.
#'@param cholesky_NG_c (Single) positive real number. The value indicates the
#'  scale parameter of the inverse gamma prior on the global scales.
#'  Expert option would be to set the scale parameter proportional to NG_a. E.g.
#'  in the case where a discrete hyperprior for NG_a is chosen, a desired
#'  proportion of let's say 0.2 is achieved by setting NG_c="0.2a" (character
#'  input!). `cholesky_NG_c` has only to be specified if
#'  `cholesky_U_prior="NG"`.
#'@param cholesky_NG_tol Minimum number that a parameter draw of one of the
#'  shrinking parameters of the normal-gamma prior can take. Prevents numerical
#'  issues that can appear when strong shrinkage is enforced if chosen to be
#'  greater than zero. `cholesky_NG_tol` has only to be specified if
#'  `cholesky_U_prior="NG"`.
#'@param cholesky_SSVS_c0 single positive number indicating the (unscaled)
#'  standard deviation of the spike component. \code{cholesky_SSVS_c0} has only
#'  to be specified if \code{choleksy_U_prior="SSVS"}.
#'   It should be that \eqn{SSVS_{c0}
#'   \ll SSVS_{c1}}!
#'@param cholesky_SSVS_c1 single positive number indicating the (unscaled)
#'  standard deviation of the slab component. \code{cholesky_SSVS_c1} has only
#'  to be specified if \code{choleksy_U_prior="SSVS"}. It should be that
#'  \eqn{SSVS_{c0} \ll SSVS_{c1}}!
#'@param cholesky_SSVS_p Either a single positive number in the range `(0,1)`
#'  indicating the (fixed) prior inclusion probability of each coefficient. Or
#'  numeric vector of length 2 with positive entries indicating the shape
#'  parameters of the Beta distribution. In that case a Beta hyperprior is
#'  placed on the prior inclusion probability. `cholesky_SSVS_p` has only to be
#'  specified if \code{choleksy_U_prior="SSVS"}.
#'@param cholesky_HMP_lambda3 numeric vector of length 2. Both entries must be
#'  positive. The first indicates the shape and the second the rate of the Gamma
#'  hyperprior on the contemporaneous coefficients. `cholesky_HMP_lambda3` has
#'  only to be specified if \code{choleksy_U_prior="HMP"}.
#'@param cholesky_normal_sds numeric vector of length \eqn{\frac{M^2-M}{2}},
#'  indicating the prior variances for the free off-diagonal elements in
#'  \eqn{U}. A single number will be recycled accordingly! Must be positive.
#'  `cholesky_normal_sds` has only to be specified if
#'  `choleksy_U_prior="normal"`.
#'@param expert_sv_offset ... Do not use!
#'@param quiet logical indicating whether informative output should be omitted.
#'@param ... Do not use!
#'
#'@seealso [specify_prior_phi()].
#'
#'@return Object of class `bayesianVARs_prior_sigma`.
#'@export
#'
#'@references Kastner, G. (2019). Sparse Bayesian Time-Varying Covariance
#'  Estimation in Many Dimensions \emph{Journal of Econometrics}, \bold{210}(1),
#'  98--115, \doi{10.1016/j.jeconom.2018.11.007}
#'
#'@references Kastner, G., Frühwirth-Schnatter, S., and Lopes, H.F. (2017).
#'  Efficient Bayesian Inference for Multivariate Factor Stochastic Volatility
#'  Models. \emph{Journal of Computational and Graphical Statistics},
#'  \bold{26}(4), 905--917, \doi{10.1080/10618600.2017.1322091}.
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # examples with stochastic volatility (heteroscedasticity) -----------------
#' # factor-decomposition with 2 factors and colwise normal-gamma prior on the loadings
#' sigma_factor_cng_sv <- specify_prior_sigma(data = data, type = "factor",
#' factor_factors = 2L, factor_priorfacloadtype = "colwiseng", factor_heteroskedastic = TRUE)
#'
#' # cholesky-decomposition with Dirichlet-Laplace prior on U
#' sigma_cholesky_dl_sv <- specify_prior_sigma(data = data, type = "cholesky",
#' cholesky_U_prior = "DL", cholesky_DL_a = 0.5, cholesky_heteroscedastic = TRUE)
#'
#' # examples without stochastic volatility (homoscedasticity) ----------------
#' # factor-decomposition with 2 factors and colwise normal-gamma prior on the loadings
#' sigma_factor_cng <- specify_prior_sigma(data = data, type = "factor",
#' factor_factors = 2L, factor_priorfacloadtype = "colwiseng",
#' factor_heteroskedastic = FALSE, factor_priorhomoskedastic = matrix(c(0.5,0.5),
#' ncol(data), 2))
#'
#' # cholesky-decomposition with Horseshoe prior on U
#' sigma_cholesky_dl <- specify_prior_sigma(data = data, type = "cholesky",
#' cholesky_U_prior = "HS", cholesky_heteroscedastic = FALSE)
#'
#'\donttest{
#' # Estimate model with your prior configuration of choice
#' mod <- bvar(data, prior_sigma = sigma_factor_cng_sv, quiet = TRUE)
#'}
specify_prior_sigma <- function(data=NULL,
                                M = ncol(data),
                                type = c("factor", "cholesky"),
                                factor_factors = 1L,
                                factor_restrict = c("none", "upper"),
                                factor_priorfacloadtype = c("rowwiseng", "colwiseng", "normal"),
                                factor_priorfacload = 0.1,
                                factor_facloadtol = 1e-18,
                                factor_priorng = c(1,1),
                                factor_priormu = c(0,10),
                                factor_priorphiidi = c(10, 3),
                                factor_priorphifac = c(10, 3),
                                factor_priorsigmaidi = 1,
                                factor_priorsigmafac = 1,
                                factor_priorh0idi = "stationary",
                                factor_priorh0fac = "stationary",
                                factor_heteroskedastic = TRUE,
                                factor_priorhomoskedastic = NA,
                                factor_interweaving = 4,
                                cholesky_U_prior = c("HS", "DL", "R2D2", "NG", "SSVS", "normal", "HMP"),
                                cholesky_U_tol = 1e-18,
                                cholesky_heteroscedastic = TRUE,
                                cholesky_priormu = c(0,100),
                                cholesky_priorphi = c(20, 1.5),
                                cholesky_priorsigma2 = c(0.5, 0.5),
                                cholesky_priorh0 = "stationary",
                                cholesky_priorhomoscedastic = as.numeric(NA),
                                cholesky_DL_a = "1/n",
                                cholesky_DL_tol = 0,
                                cholesky_R2D2_a =0.4,
                                cholesky_R2D2_b = 0.5,
                                cholesky_R2D2_tol=0,
                                cholesky_NG_a = .5,
                                cholesky_NG_b = .5,
                                cholesky_NG_c = .5,
                                cholesky_NG_tol = 0,
                                cholesky_SSVS_c0 = 0.001,
                                cholesky_SSVS_c1 = 1,
                                cholesky_SSVS_p = 0.5,
                                cholesky_HMP_lambda3 = c(0.01,0.01),
                                cholesky_normal_sds = 10,
                                expert_sv_offset = 0,
                                quiet = FALSE,
                                ...){

  if(is.null(data) & !is.numeric(M)){
    stop("\nEither provide 'data', the data for the VAR to be estimated, or 'M', the dimensionality of the data (number of time series)!\n")
  }

  if(M<2 | M%%1!=0){
    stop("\nArgument 'M' must be an integer greater or equal to 2.\n")
  }
  M <- as.integer(M)

  if(!is.null(data)){
    if(ncol(data) != M){
      warning(paste0("\nArgument 'M' does not coincide with 'ncol(data)'. Setting M=", ncol(data),"!\n"))
      M <- ncol(data)
    }
  }

  optionals <- list(...)
  if(length(optionals)>0){
    warning(paste0("\nYou provided an additional argument. Additional argument:\n", if(is.null(names(optionals))) NULL else paste0(names(optionals),"="), optionals ))
  }

  # error checks 'type'
  if(length(type)<=2L & is.character(type)){
    if(length(type)==2L) type <- type[1]
    if(!(type %in% c("factor", "cholesky"))){
      stop("type must be either 'factor' or 'cholesky'.")
    }
  }else stop("type must be either 'factor' or 'cholesky'.")

  # placeholder for cpp (cpp function expects a list with all the following elements)
  # (e.g. even if type is specified as cholesky, cpp function requires the 'factor_list')
  sv_priormu <- sv_priorphi <- sv_priorh0 <- numeric(1L)
  sv_priorsigma2 <- matrix(as.numeric(NA),1,1)
  sv_heteroscedastic <- logical(1L)
  cholesky_list <- list(
    #cholesky_heteroscedastic = logical(1L),
    cholesky_priorhomoscedastic = matrix(as.numeric(NA),1,1),
    cholesky_U_prior = character(1L),
    cholesky_U_tol = numeric(1L),
    ## GL priors
    cholesky_GL_tol = double(1L),
    cholesky_a = double(1L),
    cholesky_b = double(1L),
    cholesky_c = double(1L),
    cholesky_GT_vs = double(1L),
    cholesky_GT_priorkernel = character(1L),
    cholesky_a_vec = double(1L),
    cholesky_a_weight = double(1L),
    cholesky_norm_consts = double(1L),
    cholesky_c_vec = double(1),
    cholesky_c_rel_a = logical(1L),
    cholesky_GT_hyper = logical(1),
    #DL
    cholesky_DL_hyper = logical(1L),
    cholesky_DL_plus = logical(1L),
    #SSVS
    cholesky_SSVS_tau0 = double(1L),
    cholesky_SSVS_tau1 = double(1L),
    cholesky_SSVS_s_a = double(1L),
    cholesky_SSVS_s_b = double(1L),
    cholesky_SSVS_hyper = logical(1L),
    cholesky_SSVS_p = double(1L),
    #HM
    cholesky_lambda_3 = double(1L),
    cholesky_sv_offset = double(1L))
  factor_list <- list(factor_factors = integer(1L),
                      factor_restrinv = matrix(1L,1,1),
                      factor_ngprior = logical(1L),
                      factor_columnwise = logical(1L),
                      factor_shrinkagepriors = list(a = double(1L),
                                                    c = double(1L),
                                                    d = double(1L)),
                      factor_facloadtol = numeric(1L),
                      factor_interweaving = integer(1L),
                      #factor_heteroskedastic = logical(1L),
                      factor_priorhomoskedastic = matrix(as.numeric(NA),1,1),
                      factor_starttau2 = matrix(as.numeric(NA), 1,1)
  )

  if(type == "factor"){
    if(!quiet){
      cat("\nSince argument 'type' is specified with 'factor', all arguments starting with 'cholesky_' are being ignored.\n")
    }

    # error checks for 'factor_factors'
    if (!is.numeric(factor_factors) | factor_factors < 0) {
      stop("Argument 'factor_factors' (number of latent factor_factors) must be a single number >= 0.")
    } else {
      factor_list$factor_factors <- as.integer(factor_factors)
    }

    # error checks for factor_interweaving
    if (is.numeric(factor_interweaving) && length(factor_interweaving) == 1) {
      factor_list$factor_interweaving <- as.integer(factor_interweaving)
    } else {
      stop("Argument 'factor_interweaving' must contain a single numeric value.")
    }

    if (factor_interweaving != 0 & factor_interweaving != 1 & factor_interweaving != 2 & factor_interweaving != 3 & factor_interweaving != 4 & factor_interweaving != 5 & factor_interweaving != 6 & factor_interweaving != 7) {
      stop("Argument 'factor_interweaving' must be one of: 0, 1, 2, 3, 4.")
    }

    # error checks 'factor_restrict'
    if(length(factor_restrict)<=2L & is.character(factor_restrict)){
      if(length(factor_restrict)==2L) factor_restrict <- factor_restrict[1]
      if(!(factor_restrict %in% c("none", "upper"))){
        stop("factor_restrict must be either 'none' or 'upper'.")
      }
    }else stop("factor_restrict must be either 'none' or 'upper'.")
    restr <- matrix(FALSE, nrow = M, ncol = factor_factors)
    if (factor_restrict == "upper") restr[upper.tri(restr)] <- TRUE

    if (factor_interweaving %in% c(1, 2) && any(diag(restr) == TRUE)) {
      stop("Setting 'factor_interweaving' to either 1 or 2 and restricting the diagonal elements of the factor loading matrix are not allowed at the same time.")
    }
    # factorstochvol sampler interpretes 0 as restricted and 1 as unrestricted
    factor_list$factor_restrinv <- matrix(as.integer(!restr), nrow = nrow(restr), ncol = ncol(restr))


    # error checks for 'factor_priorfacloadtype'
    if(length(factor_priorfacloadtype)<=3L & is.character(factor_priorfacloadtype)){
      if(length(factor_priorfacloadtype)>1L) factor_priorfacloadtype <- factor_priorfacloadtype[1]
      if(!(factor_priorfacloadtype %in% c("rowwiseng", "colwiseng", "normal"))){
        stop("factor_priorfacloadtype must be either 'rowwiseng' or 'colwiseng' or 'normal'.")
      }
    }else stop("factor_priorfacloadtype must be either 'rowwiseng' or 'colwiseng' or 'normal'.")
    if (factor_priorfacloadtype == "normal") {
      #factor_pfl <- 1L
      factor_list$factor_ngprior <- FALSE
    } else if (factor_priorfacloadtype == "rowwiseng") {
      factor_list$factor_ngprior <- TRUE
      factor_list$factor_columnwise <- FALSE
    } else if (factor_priorfacloadtype == "colwiseng") {
      factor_list$factor_ngprior <- TRUE
      factor_list$factor_columnwise <- TRUE
    }

    # error checks for 'factor_priorng'
    if (!is.numeric(factor_priorng) | length(factor_priorng) != 2 | any(factor_priorng <= 0)) {
      stop("Argument 'factor_priorng' (prior hyperhyperparameters for Normal-Gamma prior) must be numeric and of length 2.")
    }
    cShrink <- factor_priorng[1]
    dShrink <- factor_priorng[2]

    # error checks for 'factor_priorfacload'
    if(!is.numeric(factor_priorfacload) | any(factor_priorfacload <= 0)) {
      stop("Argument 'priorfacload' must be numeric and positive.")
    }

    if(is.matrix(factor_priorfacload)) {
      if(nrow(factor_priorfacload) != M || ncol(factor_priorfacload) != factor_factors) {
        stop("If argument 'priorfacload' is a matrix, it must be of appropriate dimensions.")
      }
      if (factor_priorfacloadtype == "normal") {
        factor_starttau2 <- factor_priorfacload^2
        aShrink <- as.numeric(NA)
        cShrink <- as.numeric(NA)
        dShrink <- as.numeric(NA)
      } else if (factor_priorfacloadtype == "rowwiseng") {
        factor_starttau2 <- matrix(1, nrow = M, ncol = factor_factors)
        aShrink <- factor_priorfacload[,1]
        warning("Only first column of 'priorfacload' is used.'")
        cShrink <- rep(cShrink, M)
        dShrink <- rep(dShrink, M)
      } else if (factor_priorfacloadtype == "colwiseng") {
        factor_starttau2 <- matrix(1, nrow = M, ncol = factor_factors)
        aShrink <- factor_priorfacload[1,]
        warning("Only first row of 'priorfacload' is used.'")
        cShrink <- rep(cShrink, factor_factors)
        dShrink <- rep(dShrink, factor_factors)
      } else if (factor_priorfacloadtype == "dl") {
        stop("'dl'prior for factorloading is not supported by bayesianVARs!")
        # factor_pfl <- 4L
        # factor_starttau2 <- matrix(1, nrow = M, ncol = factor_factors)
        # aShrink <- factor_priorfacload[1,1]
        # warning("Only first element of 'priorfacload' is used.'")
        # cShrink <- NA
        # dShrink <- NA
      }
    } else {
      if (length(factor_priorfacload) != 1) {
        stop("If argument 'priorfacload' isn't a matrix, it must be a single number.")
      }
      if (factor_priorfacloadtype == "normal") {
        factor_starttau2 <- matrix(factor_priorfacload^2, nrow = M, ncol = factor_factors)
        aShrink <- as.numeric(NA)
        cShrink <- as.numeric(NA)
        dShrink <- as.numeric(NA)
      } else if (factor_priorfacloadtype == "rowwiseng") {
        factor_starttau2 <- matrix(1, nrow = M, ncol = factor_factors)
        aShrink <- rep(factor_priorfacload, M)
        cShrink <- rep(cShrink, M)
        dShrink <- rep(dShrink, M)
      } else if (factor_priorfacloadtype == "colwiseng") {
        factor_starttau2 <- matrix(1, nrow = M, ncol = factor_factors)
        aShrink <- rep(factor_priorfacload, factor_factors)
        cShrink <- rep(cShrink, factor_factors)
        dShrink <- rep(dShrink, factor_factors)
      } else if (factor_priorfacloadtype == "dl") {
        stop("'dl' prior for factorloading is not supported by bayesianVARs!")
        # factor_pfl <- 4L
        # factor_starttau2 <- matrix(1, nrow = M, ncol = factor_factors)
        # aShrink <- factor_priorfacload
        # cShrink <- NA
        # dShrink <- NA
      }
    }
    factor_list$factor_shrinkagepriors <- list(a = aShrink,
                                               c = cShrink,
                                               d = dShrink)
    factor_list$factor_starttau2 <- factor_starttau2

    # error checks for 'factor_facloadtol'
    if(factor_facloadtol < 0){
      stop("Argument 'factor_facloadtol' (tolerance for the factor loadings) must be >=0.")
    }
    factor_list$factor_facloadtol <- factor_facloadtol

    # error checks for 'factor_priormu'
    if (!is.numeric(factor_priormu) | length(factor_priormu) != 2) {
      stop("Argument 'factor_priormu' (mean and sd for the Gaussian prior for mu) must be numeric and of length 2.")
    }
    if(any(factor_heteroskedastic==TRUE)){
      sv_priormu <- factor_priormu
    }


    # error checks for 'factor_priorphiidi'
    if (!is.numeric(factor_priorphiidi) | length(factor_priorphiidi) != 2) {
      stop("Argument 'factor_priorphiidi' (shape1 and shape2 parameters for the Beta prior for (phi+1)/2) must be numeric and of length 2.")
    }

    # error checks for 'factor_priorphifac'
    if (!is.numeric(factor_priorphifac) | length(factor_priorphifac) != 2) {
      stop("Argument 'factor_priorphifac' (shape1 and shape2 parameters for the Beta prior for (phi+1)/2) must be numeric and of length 2.")
    }
    if(any(factor_heteroskedastic==TRUE)){
      sv_priorphi <- c(factor_priorphiidi, factor_priorphifac)
    }

    # error checks for 'factor_priorsigmaidi'
    if (!is.numeric(factor_priorsigmaidi) | any(factor_priorsigmaidi <= 0)) {
      stop("Argument 'factor_priorsigmaidi' (scaling of the chi-squared(df = 1) prior for sigma^2) must be numeric and > 0.")
    }
    if(!(length(factor_priorsigmaidi)==1 | length(factor_priorsigmaidi) == M)){
      stop("Argument 'factor_priorsigmaidi' must be either of length 1 or M.")
    }
    factor_priorsigmaidi <- rep_len(factor_priorsigmaidi, M)

    # error checks for 'factor_priorsigmafac'
    if (!is.numeric(factor_priorsigmafac) | any(factor_priorsigmafac <= 0)) {
      stop("Argument 'factor_priorsigmafac' (scaling of the chi-squared(df = 1) prior for sigma^2) must be numeric and > 0.")
    }

    if (length(factor_priorsigmafac) == 1) {
      factor_priorsigmafac <- rep(factor_priorsigmafac, factor_factors)
    } else if (length(factor_priorsigmafac) == factor_factors) {
      factor_priorsigmafac <- factor_priorsigmafac
    } else {
      stop("Argument 'factor_priorsigmafac' (scaling of the chi-squared(df = 1) prior for sigma^2) must of length 1 or factor_factors")
    }

    factor_priorsigma <- c(factor_priorsigmaidi, factor_priorsigmafac)
    # factorstochvol specifies chi-squared prior, stochvol however is parametrized in gamma:
    if(any(factor_heteroskedastic==TRUE)){
      sv_priorsigma2 <- cbind(0.5,0.5/factor_priorsigma)
    }

    # error checks for factor_priorh0idi
    factor_priorh0idi[remember <- factor_priorh0idi == "stationary"] <- -1
    factor_priorh0idi[!remember] <- as.numeric(factor_priorh0idi[!remember])^2
    factor_priorh0idi <- as.numeric(factor_priorh0idi)
    if (any(factor_priorh0idi[!remember] < 0)) stop("Argument 'priorh0idi' must not contain negative values.")
    if(!(length(factor_priorh0idi) == 1 | length(factor_priorh0idi) == M)){
      stop("Argument 'factor_priorh0idi' must be either of length 1 or M.")
    }
    factor_priorh0idi <- rep_len(factor_priorh0idi,M)

    # error checks for factor_priorh0fac
    if (length(factor_priorh0fac) == 1) factor_priorh0fac <- rep(factor_priorh0fac, factor_factors)
    if (length(factor_priorh0fac) != factor_factors) stop("Argument 'factor_priorh0fac' must be of length 1 or factor_factors.")
    factor_priorh0fac[remember <- factor_priorh0fac == "stationary"] <- -1
    factor_priorh0fac[!remember] <- as.numeric(factor_priorh0fac[!remember])^2
    factor_priorh0fac <- as.numeric(factor_priorh0fac)
    if (any(factor_priorh0fac[!remember] < 0)) stop("Argument 'factor_priorh0fac' must not contain negative values.")

    sv_priorh0 <- c(factor_priorh0idi, factor_priorh0fac)

    # Some error checking for factor_heteroskedastic
    if (length(factor_heteroskedastic) == 1) factor_heteroskedastic <- rep(factor_heteroskedastic, M + factor_factors)
    if (length(factor_heteroskedastic) == 2) factor_heteroskedastic <- c(rep(factor_heteroskedastic[1], M), rep(factor_heteroskedastic[2], factor_factors))
    if (length(factor_heteroskedastic) != M + factor_factors) stop("Argument 'factor_heteroskedastic' must be of length 1, 2, or (ncol(y) + factor_factors).")
    if (!is.logical(factor_heteroskedastic)) stop("Argument 'factor_heteroskedastic' must be a vector containing only logical values.")
    if (is.null(factor_heteroskedastic)) factor_heteroskedastic <- rep(TRUE, M + factor_factors)
    if (!all(factor_heteroskedastic[M+seq_len(factor_factors)])) {
      if (factor_interweaving == 2L || factor_interweaving == 4L) {
        if(!quiet){
          cat("\nCannot do deep factor_interweaving if (some) factor_factors are homoskedastic. Setting 'factor_interweaving' to 3.\n")
        }
        factor_list$factor_interweaving <- 3L
      }
    }

    if (!all(factor_heteroskedastic)) {
      if (any(is.na(factor_priorhomoskedastic))) {
        factor_priorhomoskedastic <- matrix(c(1.1, 0.055), byrow = TRUE, nrow = M, ncol = 2)
        if(!quiet){
          cat(paste0("\nArgument 'factor_priorhomoskedastic' must be a matrix with dimension c(M, 2) if some of the
		  elements of 'factor_heteroskedastic' are FALSE. Setting factor_priorhomoskedastic to a matrix with
		  all rows equal to c(", factor_priorhomoskedastic[1], ", ", factor_priorhomoskedastic[2], ").\n"))
        }

      }
      if (!is.matrix(factor_priorhomoskedastic) || nrow(factor_priorhomoskedastic) != M ||
          ncol(factor_priorhomoskedastic) != 2 || any(factor_priorhomoskedastic <= 0)) {
        stop("Argument 'factor_priorhomoskedastic' must be a matrix with positive entries and dimension c(M, 2).")
      }
    }
    sv_heteroscedastic <- factor_heteroskedastic
    factor_list$factor_priorhomoskedastic <- as.matrix(factor_priorhomoskedastic)

  }else if(type == "cholesky"){

    if(!quiet){
      cat("\nSince argument 'type' is specified with 'cholesky', all arguments starting with 'factor_' are being ignored.\n")
    }

    n_U <- (M^2-M)/2 # number of free off diagonal elements in U

    # error checks for cholesky_heteroscedastic
    if(!is.logical(cholesky_heteroscedastic) | length(cholesky_heteroscedastic) > 1L){
      stop("Argument 'cholesky_heteroscedastic' must be a single logical.")
    }
    sv_heteroscedastic <- rep_len(cholesky_heteroscedastic, M)

    # error checks for cholesky_priormu
    if(!is.numeric(cholesky_priormu) | length(cholesky_priormu) != 2){
      stop("Argument 'choleksy_priormu' must be a numeric vector of length 2, where the
                      second element must be posistive.")
    }
    if(cholesky_priormu[2]<0){
      stop("Argument 'choleksy_priormu' must be a numeric vector of length 2, where the
                      second element must be posistive.")
    }
    if(cholesky_heteroscedastic){
      sv_priormu <- cholesky_priormu
    }


    #error checks for cholesky_priorphi
    if(!is.numeric(cholesky_priorphi) | any(cholesky_priorphi<0) | length(cholesky_priorphi) != 2){
      stop("Argument 'cholesky_priorphi' must be a  strictly positive numeric vector of length 2.")
    }
    if(cholesky_heteroscedastic){
      sv_priorphi <- cholesky_priorphi
    }

    # error checks for cholesky_priorsigma2
    if(is.matrix(cholesky_priorsigma2)){
      if(ncol(cholesky_priorsigma2)!=2 | any(cholesky_priorsigma2<0) | nrow(cholesky_priorsigma2) != M){
        stop("Argument 'cholesky_priorsigma2' must be either a  strictly positive numeric vector of length 2 or a matrix of dimension 'c(M,2)'.")
      }
    }else if(length(cholesky_priorsigma2) == 2){
      cholesky_priorsigma2 <- matrix(cholesky_priorsigma2, M, 2, byrow = TRUE)
    }else{
      stop("Argument 'cholesky_priorsigma2' must be either a  strictly positive numeric vector of length 2 or a matrix of dimension 'c(M,2)'.")
    }
    if(cholesky_heteroscedastic){
      sv_priorsigma2 <- cholesky_priorsigma2
    }

    # error checks for cholesky_priorh0
    cholesky_priorh0[remember <- cholesky_priorh0 == "stationary"] <- -1
    cholesky_priorh0[!remember] <- as.numeric(cholesky_priorh0[!remember])^2
    cholesky_priorh0 <- as.numeric(cholesky_priorh0)
    if (any(cholesky_priorh0[!remember] < 0)) stop("Argument 'priorh0idi' must not contain negative values.")
    if(!(length(cholesky_priorh0) == 1 | length(cholesky_priorh0) == M)){
      stop("Argument 'cholesky_priorh0' must be either of length 1 or M.")
    }
    sv_priorh0 <- rep_len(cholesky_priorh0,M)

    # error checks for expert_sv_offset
    if(length(expert_sv_offset)>1L){
      if(!is.null(data) & length(expert_sv_offset) !=M){
        stop("Argument 'expert_sv_offset' must be either a single non-negative number, or a vector of length 'ncol(data)' with non-negative entries.")
      }
    }
    if(any(expert_sv_offset<0) | !is.numeric(expert_sv_offset)){
      stop("Argument 'expert_sv_offset' must be greater than zero.")
    }

    cholesky_list$cholesky_sv_offset <- rep_len(expert_sv_offset, M)

    # error checks cholesky_priorhomoscedastic
    if(!cholesky_heteroscedastic){
      ph_error <- paste0("cholesky_priorHomoscedastic must be either a numeric matrix of dimension c(M,2),
             or a numeric vector of length 2, where all entries are greater than 0. \n")
      if(length(cholesky_priorhomoscedastic)==1L){
        if(is.na(cholesky_priorhomoscedastic)){
          cholesky_priorhomoscedastic <- c(.01,.01)
          if(!quiet){
            cat("\nArgument 'cholesky_priorhomoscedastic' not specified. Setting both shape and rate of inverse gamma prior equal to 0.01.\n")
          }
        }else stop(ph_error)
      }
      if(!identical(dim(cholesky_priorhomoscedastic), as.integer(c(M,2)))){
        if(length(cholesky_priorhomoscedastic) == 2 & is.numeric(cholesky_priorhomoscedastic) &
           all(cholesky_priorhomoscedastic>0)){
          cholesky_priorhomoscedastic <- matrix(cholesky_priorhomoscedastic, M, 2, byrow = TRUE)
        }else {
          stop(ph_error)
        }
      }else if(identical(dim(cholesky_priorhomoscedastic), as.integer(c(M,2))) &
               is.numeric(cholesky_priorhomoscedastic)){
        if(!all(cholesky_priorhomoscedastic>0)){
          stop(ph_error)
        }
      }else{
        stop(ph_error)
      }
    }
    cholesky_list$cholesky_priorhomoscedastic <- as.matrix(cholesky_priorhomoscedastic)

    # error checks for cholesky_U_prior
    cholesky_U_prior <- cholesky_U_prior[1]
    if(!(cholesky_U_prior %in% c("DL", "HMP", "SSVS", "normal", "R2D2", "NG", "HS"))){
      stop("Argument 'cholesky_U_prior' must be one of 'DL', 'SSVS', 'HMP' or 'normal'. \n")
    }
    cholesky_list$cholesky_U_prior <- cholesky_U_prior

    # error checks for 'cholesky_U_tol'
    if(cholesky_U_tol < 0){
      stop("Argument 'cholesky_U_tol' (tolerance for the factor loadings) must be >=0.")
    }
    cholesky_list$cholesky_U_tol <- cholesky_U_tol

    if(cholesky_U_prior == "DL"){
      text <- c("Argument 'cholesky_DL_a' must be either a single positive numeric or '1/n'. \n ")
      if(is.numeric(cholesky_DL_a) & any(cholesky_DL_a<=0)) stop(text)
      if(is.character(cholesky_DL_a)){#"hyperprior",
        if(length(cholesky_DL_a)>1 | !(cholesky_DL_a %in% c("1/n")) ){
          stop(text)
        }
      }

      cholesky_list$cholesky_U_prior <- cholesky_U_prior
      cholesky_list$cholesky_a <- cholesky_DL_a
      cholesky_list$cholesky_GL_tol <- cholesky_DL_tol

      if(is.numeric(cholesky_list$cholesky_a) & length(cholesky_list$cholesky_a) == 1L){
        cholesky_list$cholesky_DL_hyper <- FALSE
      }else if(is.character(cholesky_list$cholesky_a)) {
        if(cholesky_list$cholesky_a == "1/n"){
          cholesky_list$cholesky_a <- 1/n_U
          cholesky_list$cholesky_DL_hyper <- FALSE
        }else stop(text)

      }else if(is.matrix(cholesky_list$cholesky_a)){

        if(ncol(cholesky_list$cholesky_a)!=2){
          stop("If you specify 'DL_a' as a matrix, the first column represents
             the support points and the second column the weights of a discrete
             hyperprior on 'DL_a' !")
        }

        cholesky_list$cholesky_DL_hyper <- TRUE
        cholesky_list$cholesky_a_vec <- cholesky_list$cholesky_a[,1]
        cholesky_list$cholesky_a_weight <- cholesky_list$cholesky_a[,2]
        # precompute log normalizing constants of hyperprior
        cholesky_list$cholesky_norm_consts <- 0.5^cholesky_list$cholesky_a_vec -
          lgamma(cholesky_list$cholesky_a_vec)
        cholesky_list$cholesky_a <- cholesky_list$cholesky_a_vec[1] #initial value

      }
      # else if(cholesky_list$cholesky_a == "hyperprior"){
      #   stop("'cholesky_a' cannot be 'hyperperior' anymore.")
      # priorSigma_in$cholesky_DL_hyper <- TRUE
      #
      # grid_L <- 1000
      # priorSigma_in$cholesky_a_vec <- seq(1/(n_U),1/2,length.out = grid_L)
      # #priorSigma_in$cholesky_prep1 <- priorSigma_in$cholesky_b_vec - 1
      # #priorSigma_in$cholesky_prep2 <- lgamma(n_U*priorSigma_in$cholesky_b_vec) - n_U*lgamma(priorSigma_in$cholesky_b_vec)
      # priorSigma_in$cholesky_norm_consts <- 0.5^priorSigma_in$cholesky_a_vec -
      #   lgamma(priorSigma_in$cholesky_a_vec)
      # priorSigma_in$cholesky_a_weight <- rep(1,grid_L)
      # priorSigma_in$cholesky_a <- 1/2 # initial value
      # }

      if(!exists("cholesky_DL_plus", optionals) || base::isFALSE(optionals$cholesky_DL_plus)){
        cholesky_list$cholesky_DL_plus <- FALSE
      }else if(optionals$cholesky_DL_plus){
        cholesky_list$cholesky_DL_plus <- TRUE
        if(!exists("DL_b", optionals)){
          cholesky_list$cholesky_b <- 0.5
        }else{
          cholesky_list$cholesky_b <- cholesky_list$cholesky_DL_b
        }
        if(!exists("DL_c", optionals)){
          cholesky_list$cholesky_c <- 0.5*cholesky_list$cholesky_a
        }else{
          cholesky_list$cholesky_c <- cholesky_list$cholesky_DL_c
        }
      }else{
        stop("Never heard of DL_plus?")
      }

    }else if(cholesky_U_prior == "R2D2"){

      cholesky_list$cholesky_U_prior <- "GT"
      cholesky_list$cholesky_b <- cholesky_R2D2_b
      cholesky_list$cholesky_a <- cholesky_R2D2_a
      cholesky_list$cholesky_c = "0.5*a"
      cholesky_list$cholesky_GT_vs <- 1/2
      cholesky_list$cholesky_GT_priorkernel <- "exponential"
      cholesky_list$cholesky_GL_tol <- cholesky_R2D2_tol

    }else if(cholesky_U_prior == "NG"){

      cholesky_list$cholesky_U_prior <- "GT"
      cholesky_list$cholesky_a <- cholesky_NG_a
      cholesky_list$cholesky_b <- cholesky_NG_b
      cholesky_list$cholesky_c <- cholesky_NG_c
      cholesky_list$cholesky_GT_vs <- 1
      cholesky_list$cholesky_GT_priorkernel <- "normal"
      cholesky_list$cholesky_GL_tol <- cholesky_NG_tol

    }else if(cholesky_U_prior == "SSVS"){
      if(!(cholesky_SSVS_c0>0 & cholesky_SSVS_c1>0)){
        stop("'cholesky_SSVS_c0' and 'cholesky_SSVS_c1' must be positive numeric values.")
      }
      if(length(cholesky_SSVS_p)==2L){
        cholesky_SSVS_sa <- cholesky_SSVS_p[1]
        cholesky_SSVS_sb <- cholesky_SSVS_p[2]
        cholesky_SSVS_p <- 0.5 # initial value
        cholesky_SSVS_hyper <- TRUE
      }else if(length(cholesky_SSVS_p)==1L){
        cholesky_SSVS_p <- cholesky_SSVS_p
        cholesky_SSVS_sa <- cholesky_SSVS_sb <- NA
        cholesky_SSVS_hyper <- FALSE
      }else{
        stop("cholesky_SSVS_p must be either numeric vector of length 1L or 2L!")
      }
      cholesky_list$cholesky_U_prior <- cholesky_U_prior
      # cholesky_list$cholesky_SSVS_c0 <- cholesky_SSVS_c0
      # cholesky_list$cholesky_SSVS_c1 <- cholesky_SSVS_c1
      cholesky_list$cholesky_SSVS_tau0 <- rep(cholesky_SSVS_c0, n_U)
      cholesky_list$cholesky_SSVS_tau1 <- rep(cholesky_SSVS_c1, n_U)
      cholesky_list$cholesky_SSVS_s_a <- cholesky_SSVS_sa
      cholesky_list$cholesky_SSVS_s_b <- cholesky_SSVS_sb
      cholesky_list$cholesky_SSVS_p <- rep_len(cholesky_SSVS_p, n_U)
      cholesky_list$cholesky_SSVS_hyper <- cholesky_SSVS_hyper


    }else if(cholesky_U_prior == "normal"){
      if(!(all(cholesky_normal_sds>0))){
        stop("'cholesky_normal_sds' must be positive. \n")
      }
      cholesky_list$cholesky_U_prior <- cholesky_U_prior
      cholesky_list$cholesky_V_i <- rep_len(cholesky_normal_sds^2, n_U) # bvar expects variances!

    }else if(cholesky_U_prior == "HMP"){

      cholesky_list$cholesky_U_prior <- cholesky_U_prior
      cholesky_list$cholesky_lambda_3 <- cholesky_HMP_lambda3

    }else if(cholesky_U_prior == "HS"){

      cholesky_list$cholesky_U_prior <- cholesky_U_prior

    }

    if(cholesky_list$cholesky_U_prior == "GT"){

      if(is.matrix(cholesky_list$cholesky_a)){
        if(ncol(cholesky_list$cholesky_a)==2){
          cholesky_list$cholesky_GT_hyper <- TRUE
          cholesky_list$cholesky_a_vec <- cholesky_list$cholesky_a[,1]
          cholesky_list$cholesky_a_weight <- cholesky_list$cholesky_a[,2]
          cholesky_list$cholesky_norm_consts <- lgamma(cholesky_list$cholesky_a_vec)
          cholesky_list$cholesky_a <- sample(cholesky_list$cholesky_a_vec, 1, replace = TRUE, prob = cholesky_list$cholesky_a_weight) # initialize a
        }else if(ncol(cholesky_list$cholesky_a)>2){
          stop("The easiest way to specify 'R2D2_a', 'NG_a' or 'GT_a' is a single postive number!")
        }else{
          cholesky_list$cholesky_a <- as.vector(cholesky_list$cholesky_a)
        }
      }else if(is.null(dim(cholesky_list$cholesky_a))){
        cholesky_list$cholesky_GT_hyper <- FALSE
        cholesky_list$cholesky_a <- cholesky_list$cholesky_a
      }

      if(is.character(cholesky_list$cholesky_c)){
        cholesky_list$cholesky_c_rel_a <- TRUE # then c is always proportion of a (e.g. for R2D2 c=0.5a)
        mya <- cholesky_list$cholesky_a
        myc <- gsub("a","mya", cholesky_list$cholesky_c)
        cholesky_c <- eval(str2lang(myc)) # initial value
        if(base::isTRUE(cholesky_list$cholesky_GT_hyper)){
          myc2 <- gsub("a","cholesky_list$cholesky_a_vec", cholesky_list$cholesky_c) #cholesky_list$cholesky_c is still character
          cholesky_list$cholesky_c_vec <- eval(str2lang(myc2))
        }
        cholesky_list$cholesky_c <- cholesky_c # initial value to cholesky_list
      }else if(is.numeric(cholesky_list$cholesky_c)){
        cholesky_list$cholesky_c_rel_a <- FALSE
      }

    }
  }
  sv_list <- list(M=M, sv_priormu = sv_priormu, sv_priorphi = sv_priorphi, sv_priorsigma2 = sv_priorsigma2, sv_priorh0 = sv_priorh0, sv_heteroscedastic = sv_heteroscedastic)
  out <- c("type"=type, sv_list, factor_list, cholesky_list)
  class(out) <- "bayesianVARs_prior_sigma"
  return(out)
}

#'Specify prior on PHI
#'
#'Configures prior on PHI, the matrix of reduced-form VAR coefficients.
#'
#'For details concerning prior-elicitation for VARs please see Gruber & Kastner
#'(2023).
#'
#'Currently one can choose between six hierarchical shrinkage priors and a
#'normal prior: `prior="HS"` stands for the Horseshoe-prior, `prior="R2D2` for
#'the R\eqn{^2}-induced-Dirichlet-decompostion-prior, `prior="NG"` for the
#'normal-gamma-prior, `prior="DL"` for the Dirichlet-Laplace-prior,
#'`prior="SSVS"` for the stochastic-search-variable-selection-prior,
#'`prior="HMP"` for the semi-hierarchical Minnesota prior and `prior=normal` for
#'the normal-prior.
#'
#'Semi-global shrinkage, i.e. group-specific shrinkage for pre-specified
#'subgroups of the coefficients, can be achieved through the argument
#'`global_grouping`.
#'
#'@param data Optional. Data matrix (can be a time series object). Each of
#'  \eqn{M} columns is assumed to contain a single time-series of length
#'  \eqn{T}.
#'@param M positive integer indicating the number of time-series of the VAR.
#'@param lags positive integer indicating the order of the VAR, i.e. the number
#'  of lags of the dependent variables included as predictors.
#'@param prior character, one of \code{"HS"}, \code{"R2D2"}, `"NG"`,
#'  \code{"DL"}, \code{"SSVS"}, \code{"HMP"} or \code{"normal"}.
#'@param priormean real numbers indicating the prior means of the VAR
#'  coefficients. One single number means that the prior mean of all own-lag
#'  coefficients w.r.t. the first lag equals \code{priormean} and \code{0} else.
#'  A vector of length M means that the prior mean of the own-lag coefficients
#'  w.r.t. the first lag equals \code{priormean} and \code{0} else. If
#'  \code{priormean} is a matrix of dimension \code{c(lags*M,M)}, then each of
#'  the \eqn{M} columns is assumed to contain \code{lags*M} prior means for the
#'  VAR coefficients of the respective VAR equations.
#'@param PHI_tol Minimum number that the absolute value of a VAR coefficient
#'  draw can take. Prevents numerical issues that can appear when strong
#'  shrinkage is enforced if chosen to be greater than zero.
#'@param DL_a (Single) positive real number. The value is interpreted as the
#'  concentration parameter for the local scales. Smaller values enforce heavier
#'  shrinkage. If the argument \code{global_grouping} specifies e.g. \code{k}
#'  groups, then \code{DL_a} can be a numeric vector of length \code{k} and the
#'  elements indicate the shrinkage in each group.  A matrix of dimension
#'  \code{c(s,2)} specifies a discrete hyperprior, where the first column
#'  contains \code{s} support points and the second column contains the
#'  associated prior probabilities. \code{DL_a} has only to be specified if
#'  \code{prior="DL"}.
#'@param DL_tol Minimum number that a parameter draw of one of the shrinking
#'  parameters of the Dirichlet Laplace prior can take. Prevents numerical
#'  issues that can appear when strong shrinkage is enforced if chosen to be
#'  greater than zero. \code{DL_tol} has only to be specified if
#'  \code{prior="DL"}.
#'@param R2D2_a (Single) positive real number. The value is interpreted as the
#'  concentration parameter for the local scales. Smaller values enforce heavier
#'  shrinkage. If the argument \code{global_grouping} specifies e.g. \code{k}
#'  groups, then \code{R2D2_a} can be a numeric vector of length \code{k} and
#'  the elements indicate the shrinkage in each group.  A matrix of dimension
#'  \code{c(s,2)} specifies a discrete hyperprior, where the first column
#'  contains \code{s} support points and the second column contains the
#'  associated prior probabilities. \code{R2D2_a} has only to be specified if
#'  \code{prior="R2D2"}.
#'@param R2D2_b (Single) positive real number. The value indicates the shape
#'  parameter of the inverse gamma prior on the (semi-)global scales. If the
#'  argument \code{global_grouping} specifies e.g. \code{k} groups, then
#'  \code{NG_b} can be a numeric vector of length \code{k} and the elements
#'  determine the shape parameter in each group. \code{R2D2_b} has only to be
#'  specified if \code{prior="R2D2"}.
#'@param R2D2_tol Minimum number that a parameter draw of one of the shrinking
#'  parameters of the R2D2 prior can take. Prevents numerical issues that can
#'  appear when strong shrinkage is enforced if chosen to be greater than zero.
#'  \code{R2D2_tol} has only to be specified if \code{prior="R2D2"}.
#'@param NG_a (Single) positive real number. The value is interpreted as the
#'  concentration parameter for the local scales. Smaller values enforce heavier
#'  shrinkage. If the argument \code{global_grouping} specifies e.g. \code{k}
#'  groups, then \code{NG_a} can be a numeric vector of length \code{k} and the
#'  elements indicate the shrinkage in each group.  A matrix of dimension
#'  \code{c(s,2)} specifies a discrete hyperprior, where the first column
#'  contains \code{s} support points and the second column contains the
#'  associated prior probabilities. \code{NG_a} has only to be specified if
#'  \code{prior="NG"}.
#'@param NG_b (Single) positive real number. The value indicates the shape
#'  parameter of the inverse gamma prior on the (semi-)global scales. If the
#'  argument \code{global_grouping} specifies e.g. \code{k} groups, then
#'  \code{NG_b} can be a numeric vector of length \code{k} and the elements
#'  determine the shape parameter in each group. \code{NG_b} has only to be
#'  specified if \code{prior="NG"}.
#'@param NG_c (Single) positive real number. The value indicates the scale
#'  parameter of the inverse gamma prior on the (semi-)global scales. If the
#'  argument \code{global_grouping} specifies e.g. \code{k} groups, then
#'  \code{NG_c} can be a numeric vector of length \code{k} and the elements
#'  determine the scale parameter in each group. Expert option would be to set
#'  the scale parameter proportional to \code{NG_a}. E.g. in the case where a
#'  discrete hyperprior for \code{NG_a} is chosen, a desired proportion of let's
#'  say \code{0.2} is achieved by setting \code{NG_c="0.2a"} (character input!).
#'  \code{NG_c} has only to be specified if \code{prior="NG"}.
#'@param NG_tol  Minimum number that a parameter draw of one of the shrinking
#'  parameters of the normal-gamma prior can take. Prevents numerical issues
#'  that can appear when strong shrinkage is enforced if chosen to be greater
#'  than zero. \code{NG_tol} has only to be specified if \code{prior="NG"}.
#'@param SSVS_c0 single positive number indicating the (unscaled) standard
#'  deviation of the spike component. \code{SSVS_c0} has only to be specified if
#'  \code{prior="SSVS"}. It should be that \eqn{SSVS_{c0} \ll SSVS_{c1}}!
#'  \code{SSVS_c0} has only to be specified if \code{prior="SSVS"}.
#'@param SSVS_c1 single positive number indicating the (unscaled) standard
#'  deviation of the slab component. \code{SSVS_c0} has only to be specified if
#'  \code{prior="SSVS"}. It should be that \eqn{SSVS_{c0} \ll SSVS_{c1}}!
#'@param SSVS_p Either a single positive number in the range \code{(0,1)}
#'  indicating the (fixed) prior inclusion probability of each coefficient. Or
#'  numeric vector of length 2 with positive entries indicating the shape
#'  parameters of the Beta distribution. In that case a Beta hyperprior is
#'  placed on the prior inclusion probability. \code{SSVS_p} has only to be
#'  specified if \code{prior="SSVS"}.
#'@param SSVS_semiautomatic logical. If \code{SSVS_semiautomatic=TRUE} both
#'  \code{SSVS_c0} and \code{SSVS_c1} will be scaled by the variances of the
#'  posterior of PHI under a FLAT conjugate (dependent Normal-Wishart prior).
#'  \code{SSVS_semiautomatic} has only to be specified if \code{prior="SSVS"}.
#'@param HMP_lambda1 numeric vector of length 2. Both entries must be positive.
#'  The first indicates the shape and the second the rate of the Gamma
#'  hyperprior on own-lag coefficients. \code{HMP_lambda1} has only to be
#'  specified if \code{prior="HMP"}.
#'@param HMP_lambda2 numeric vector of length 2. Both entries must be positive.
#'  The first indicates the shape and the second the rate of the Gamma
#'  hyperprior on cross-lag coefficients. \code{HMP_lambda2} has only to be
#'  specified if \code{prior="HMP"}.
#'@param normal_sds numeric vector of length \eqn{n}, where \eqn{n = lags M^2} is
#'  the number of all VAR coefficients (excluding the intercept), indicating the
#'  prior variances. A single number will be recycled accordingly! Must be
#'  positive. \code{normal_sds} has only to be specified if
#'  \code{prior="normal"}.
#'@param global_grouping One of \code{"global"}, \code{"equation-wise"},
#'  \code{"covariate-wise"}, \code{"olcl-lagwise"} \code{"fol"} indicating the
#'  sub-groups of the semi-global(-local) modifications to HS, R2D2, NG, DL and
#'  SSVS prior. Works also with user-specified indicator matrix of dimension
#'  \code{c(lags*M,M)}.  Only relevant if \code{prior="HS"}, \code{prior="DL"},
#'  \code{prior="R2D2"}, `prior="NG"` or \code{prior="SSVS"}.
#'@param ... Do not use!
#'
#'@references  Gruber, L. and Kastner, G. (2023). Forecasting macroeconomic data
#'  with Bayesian VARs: Sparse or dense? It depends!
#'  \href{https://arxiv.org/abs/2206.04902}{arXiv:2206.04902}.
#'
#'@return A `baysianVARs_prior_phi`-object.
#'@export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Horseshoe prior for a VAR(2)
#' phi_hs <- specify_prior_phi(data = data, lags = 2L ,prior = "HS")
#'
#' # Semi-global-local Horseshoe prior for a VAR(2) with semi-global shrinkage parameters for
#' # cross-lag and own-lag coefficients in each lag
#' phi_hs_sg <- specify_prior_phi(data = data, lags = 2L, prior = "HS",
#' global_grouping = "olcl-lagwise")
#'
#' # Semi-global-local Horseshoe prior for a VAR(2) with equation-wise shrinkage
#' # construct indicator matrix for equation-wise shrinkage
#' semi_global_mat <- matrix(1:ncol(data), 2*ncol(data), ncol(data),
#' byrow = TRUE)
#' phi_hs_ew <- specify_prior_phi(data = data, lags = 2L, prior = "HS",
#' global_grouping = semi_global_mat)
#' # (for equation-wise shrinkage one can also use 'global_grouping = "equation-wise"')
#'
#'\donttest{
#' # Estimate model with your prior configuration of choice
#' mod <- bvar(data, lags = 2L, prior_phi = phi_hs_sg, quiet = TRUE)
#'}
specify_prior_phi <- function(data = NULL,
                              M = ncol(data),
                              lags = 1L,
                              prior = "HS",
                              priormean = 0,
                              PHI_tol = 1e-18,
                              DL_a = "1/K", DL_tol = 0,
                              R2D2_a =0.1, R2D2_b = 0.5, R2D2_tol = 0,
                              NG_a = 0.1, NG_b = 1, NG_c = 1, NG_tol = 0,
                              SSVS_c0 = 0.01, SSVS_c1 = 100,
                              SSVS_semiautomatic = TRUE, SSVS_p=0.5,
                              HMP_lambda1 = c(0.01,0.01), HMP_lambda2 = c(0.01,0.01),
                              normal_sds = 10,
                              global_grouping="global",
                              ...){

  if(!(prior %in% c("DL", "HMP", "SSVS", "normal", "R2D2", "SL", "HS", "NG"))){
    stop("Argument 'prior' must be one of 'DL', 'HS', 'NG', 'SSVS', 'HMP' or 'normal'. \n")
  }

  if(is.null(data) & is.null(M)){
    stop("\nBoth 'data' and 'M' are missing! (One of them has to be specified.)\n")
  }
  if(!is.null(data) & (!identical(ncol(data), as.integer(M)))){
    M <- ncol(data)
    warning(paste0("\nArgument 'M' does not align with 'ncol(data)'. Setting 'M'=", M,"!\n"))
  }

  K <- lags*M

  if(prior %in% c("DL", "SSVS",  "R2D2", "HS", "NG")){
    if(is.character(global_grouping)){
      if(!(global_grouping %in% c("global", "equation-wise", "covariate-wise", "fol", "olcl-lagwise"))){
        stop("Argument 'global_grouping' must be one of 'global',
           'equation-wise', 'covariate-wise', 'olcl-lagwise' or 'fol'. \n")
      }
    }

    if(all(is.numeric(global_grouping))){
      i_mat <- global_grouping
      if(!(identical(dim(i_mat), as.integer(c(K,M))) & all(is.numeric(i_mat))) )
        stop("\nArgument 'global_grouping' is not of dimension 'c(lags*M,M)'!\n")
    }else if(global_grouping=="global"){
      i_mat <- matrix(1, K, M)
    }else if(global_grouping=="fol"){
      i_mat <- matrix(1, K, M)
      diag(i_mat[1:M,1:M]) <- 2
    }else if(global_grouping=="olcl-lagwise"){
      i_mat <- matrix(1, K, M)
      diag(i_mat[1:M,1:M]) <- 2
      if(lags>1){
        for (j in 1:(lags-1)) {
          i_mat[(j*M+1):((j+1)*M),] <- i_mat[((j-1)*M+1):(j*M),] + 2
        }
      }
    }else if(global_grouping == "equation-wise"){
      i_mat <- matrix(
        rep(1:M, K),
        K,M,
        byrow = TRUE
      )
    }else if(global_grouping == "covariate-wise"){
      i_mat <- matrix(
        rep(1:K, M),
        K,M
      )
    }else {
      stop("Something went wrong specifying 'global_grouping'.")
    }

  }else{
    i_mat_1 <- diag(M)
    i_mat_1[upper.tri(i_mat_1)] <-
      i_mat_1[lower.tri(i_mat_1)] <- -1
    i_mat <- NULL
    for (j in seq_len(lags)) {
      i_mat <- rbind(i_mat, i_mat_1 * j)
    }
  }


  if(any(priormean<0)){
    stop("\nAt least one element in 'priormean' is negative!\n")
  }
  if(length(priormean)>1L){
    if(is.vector(priormean)){
      if(length(priormean)!=M) stop("\n'length(priormean)' does not equal 'M'.\n")
      priormean <- as.vector(priormean)
      PHI0 <- matrix(0, K, M)
      PHI0[1:M,1:M] <- diag(priormean)
    }
    if(is.matrix(priormean)){
      if(ncol(priormean)!=M) stop("\n'ncol(priormean)' does not equal 'M'!\n")
      if(nrow(priormean)!=K) stop("\n'nrow(priormean)' does not equal 'lags*M'!\n")
      PHI0 <- priormean
    }
  }else if(length(priormean) == 1L){
    priormean <- rep(priormean, M)
    PHI0 <- matrix(0, K, M)
    PHI0[1:M,1:M] <- diag(M)*priormean
  }

  if(prior == "DL"){
    text <- c("Argument 'DL_a' must be either a single positive numeric or one of 'hyperprior',
           '1/K' or '1/n'. \n ")

    if(any(DL_a <= 0)) {
      stop(text)
    }else if(all(is.character(DL_a))){
      if(!(any(is.character(DL_a)) &
           any(DL_a %in% c("hyperprior", "1/K", "1/n")) &
           length(DL_a)==1)){
        stop(text)
      }
    }

    if(is.character(global_grouping)){
      if(!(global_grouping %in% c("global", "equation-wise", "covariate-wise", "fol", "olcl-lagwise"))){
        stop("Argument 'global_grouping' must be one of 'global',
           'equation-wise', 'covariate-wise', 'olcl-lagwise' or 'fol'. \n")
      }
    }
    out <- list(prior = prior, PHI_tol = PHI_tol, a = DL_a, global_grouping = global_grouping,
                GL_tol = DL_tol, ...)

  }else if(prior == "R2D2"){
    if(is.character(global_grouping)){
      if(!(global_grouping %in% c("global", "equation-wise", "covariate-wise", "fol", "olcl-lagwise"))){
        stop("Argument 'global_grouping' must be one of 'global',
           'equation-wise', 'covariate-wise', 'olcl-lagwise' or 'fol'. \n")
      }
    }

    if(any(R2D2_a <= 0)) {
      stop("R2D2_a must be strictly greater than 0!")
    }

    out <- list(prior = "GT", PHI_tol = PHI_tol, b = R2D2_b, a = R2D2_a,
                global_grouping = global_grouping, c = "0.5*a", GT_vs = 1/2,
                GT_priorkernel = "exponential", GL_tol = R2D2_tol,...)

  }else if(prior == "SSVS"){
    if(!(all(SSVS_c0>0) & all(SSVS_c1>0))){
      stop("'SSVS_c0' and 'SSVS_c1' must be positive numeric values. \n")
    }
    if(length(SSVS_p)==2L){
      SSVS_sa <- SSVS_p[1]
      SSVS_sb <- SSVS_p[2]
      SSVS_p <- 0.5 # initial value
      SSVS_hyper <- TRUE
    }else if(length(SSVS_p)==1L){
      SSVS_p <- SSVS_p
      SSVS_sa <- SSVS_sb <- NA
      SSVS_hyper <- FALSE
    }else{
      stop("SSVS_p must be either numeric vector of length 1L or 2L!")
    }
    out <- list(prior = prior, PHI_tol = PHI_tol, SSVS_c0=SSVS_c0, SSVS_c1=SSVS_c1,
                semiautomatic=SSVS_semiautomatic, SSVS_s_a=SSVS_sa,
                SSVS_s_b=SSVS_sb, SSVS_p = SSVS_p, SSVS_hyper = SSVS_hyper,
                global_grouping = global_grouping)

  }else if(prior == "normal"){
    if(!(all(normal_sds>0))){
      stop("'normal_sds' must be positive. \n")
    }
    out <- list(prior = prior, PHI_tol = PHI_tol, V_i=normal_sds^2) # note to myself: bvar expects variances!
  }else if(prior == "HMP"){
    out <- list(prior = prior, PHI_tol = PHI_tol, lambda_1 = HMP_lambda1, lambda_2 = HMP_lambda2)
  }else if(prior == "SL"){
    out <- list(prior = prior, PHI_tol = PHI_tol, ...)
  }else if(prior == "HS"){
    if(is.character(global_grouping)){
      if(!(global_grouping %in% c("global", "equation-wise", "covariate-wise", "fol", "olcl-lagwise"))){
        stop("Argument 'global_grouping' must be one of 'global',
           'equation-wise', 'covariate-wise', 'olcl-lagwise' or 'fol'. \n")
      }
    }
    out <- list(prior = prior, PHI_tol = PHI_tol, global_grouping = global_grouping)
  }else if(prior == "NG"){
    if(is.character(global_grouping)){
      if(!(global_grouping %in% c("global", "equation-wise", "covariate-wise", "fol", "olcl-lagwise"))){
        stop("Argument 'global_grouping' must be one of 'global',
           'equation-wise', 'covariate-wise', 'olcl-lagwise' or 'fol'. \n")
      }
    }
    out <- list(prior = "GT", PHI_tol = PHI_tol, a = NG_a, b = NG_b, c = NG_c, GT_vs = 1,
                GT_priorkernel = "normal",
                GL_tol = NG_tol, global_grouping = global_grouping)
  }
  out$lags <- as.integer(lags)
  out$M <- as.integer(M)
  out$PHI0 <- PHI0
  out$i_mat <- i_mat
  class(out) <- "bayesianVARs_prior_phi"
  out
}


# Predict methods ---------------------------------------------------------

#' Predict method for Bayesian VARs
#'
#' Simulates from (out-of-sample) predictive density for Bayesian VARs estimated
#' via [`bvar()`] and computes log predictive likelhoods if ex-post
#' observed data is supplied.
#'
#' @param object A `bayesianVARs_bvar` object, obtained from [`bvar()`].
#' @param ahead Integer vector (or coercible to such), indicating the number of
#'   steps ahead at which to predict.
#' @param each Single integer (or coercible to such) indicating how often should
#'   be drawn from the posterior predictive distribution for each draw that has
#'   been stored during MCMC sampling.
#' @param stable logical indicating whether to consider only those draws from
#'   the posterior that fulfill the 'stable' criterion. Default is \code{TRUE}.
#' @param simulate_predictive logical, indicating whether the posterior
#'   predictive distribution should be simulated.
#' @param LPL logical indicating whether `ahead`-step-ahead log predictive
#'   likelihoods should be computed. If \code{LPL=TRUE}, \code{Y_obs} has to be
#'   specified.
#' @param Y_obs Data matrix of observed values for computation of log predictive
#'   likelihood. Each of `ncol(object$Yraw)` columns is assumed to contain a
#'   single time-series of length \code{length(ahead)}.
#' @param LPL_VoI either integer vector or character vector of column-names
#'   indicating for which subgroup of time-series in `object$Yraw` a joint log
#'   predictive likelihood shall be computed.
#' @param ... Currently ignored!
#'
#' @seealso [`stable_bvar()`], [plot.bayesianVARs_predict()], [pairs.bayesianVARs_predict()].
#'
#' @return Object of class `bayesianVARs_predict`, a list that may contain the
#'   following elements:
#' * `predictions` array of dimensions
#'   \code{c(length(ahead), ncol(object$Yraw), each * dim(object$PHI)[3])}
#'   containing the simulations from the predictive density (if
#'   `simulate_predictive=TRUE`).
#' * `LPL` vector of length `length(ahead)` containing the
#'   log-predictive-likelihoods (taking into account the joint distribution of
#'   all variables) (if `LPL=TRUE`).
#' * `LPL_univariate` matrix of dimension `c(length(ahead), ncol(object$Yraw)`
#'   containing the marginalized univariate log-predictive-likelihoods of each
#'   series (if `LPL=TRUE`).
#' * `LPL_VoI` vector of length `length(ahead)` containing the
#'   log-predictive-likelihoods for a subset of variables (if `LPL=TRUE` and
#'   `LPL_VoI != NA`).
#' * `Yraw` matrix containing the data used for the estimation of the VAR.
#' * `LPL_draws` matrix containing the simulations of the
#'   log-predictive-likelihood (if `LPL=TRUE`).
#' * `PL_univariate_draws` array containing the simulations of the univariate
#'   predictive-likelihoods (if `LPL=TRUE`).
#' * `LPL_sub_draws` matrix containing the simulations of the
#'   log-predictive-likelihood for a subset of variables (if `LPL=TRUE` and
#'   `LPL_VoI != NA`).
#' @export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Split data in train and test
#' train <- data[1:(nrow(data)-4),]
#' test <- data[-c(1:(nrow(data)-4)),]
#'
#' # Estimate model using train data only
#' mod <- bvar(train, quiet = TRUE)
#'
#' # Simulate from 1-step to 4-steps ahead posterior predictive and compute
#' # log-predictive-likelihoods
#' predictions <- predict(mod, ahead = 1:4, LPL = TRUE, Y_obs = test)
#'
#' # Summary
#' summary(predictions)
#'
#' # Visualize via fan-charts
#' plot(predictions)
#'
#' \donttest{
#' # In order to evaluate the joint predictive density of a subset of the
#' # variables (variables of interest), consider specifying 'LPL_VoI':
#' predictions <- predict(mod, ahead = 1:4, LPL = TRUE, Y_obs = test, LPL_VoI = c("GDPC1","FEDFUNDS"))
#' predictions$LPL_VoI
#' }
predict.bayesianVARs_bvar <- function(object, ahead = 1L, each = 1L, stable = TRUE,
                        simulate_predictive = TRUE, LPL = FALSE, Y_obs = NA,
                        LPL_VoI = NA, ...){

  factor <- object$sigma_type == "factor"
  ahead <- sort(as.integer(ahead))

  if(!simulate_predictive & !LPL){
    stop("At least one of 'simulate_predictive' or 'LPL' must be set 'TRUE'!")
  }

  if(stable){
    if(!inherits(object, "bayesianVARs_bvar_stable")){
      cat("'stable=TRUE': Calling 'stable_bvar()' to discard those posterior
          draws, that do not fulfill the stable criterion.\n")
      object <- stable_bvar(object, quiet = TRUE)
      cat("\n",dim(object$PHI)[3],"stable posterior draws remaining for prediction!\n")
    }
  }

  # relevant mod settings
  sv_indicator <- which(object$heteroscedastic==TRUE)
  intercept <- object$intercept

  # data preparation
  variables <- colnames(object$Y)
  # draws <- dim(object$PHI)[1]
  # lags <- object$lags # number of lags
  M <- ncol(object$Y)
  #if(any(SV)) {
  # extract sv parameters
  sv_mu <- matrix(object$sv_para[1,sv_indicator,], length(sv_indicator), dim(object$sv_para)[3])
  sv_phi <- matrix(object$sv_para[2,sv_indicator,], length(sv_indicator), dim(object$sv_para)[3])
  sv_sigma <- matrix(object$sv_para[3,sv_indicator,], length(sv_indicator), dim(object$sv_para)[3])
  #}
  # extract current state of logvariance (in case of homoscedasticity one could pick any, since constant...)
  sv_h_T <- object$logvar[dim(object$logvar)[1],,]
  X_T_plus_1 <- as.numeric(object$datamat[nrow(object$datamat), 1:(object$lags*M)])
  if(object$intercept) X_T_plus_1 <- c(X_T_plus_1, 1)


  LPL_subset <- logical(1L)
  VoI <- 1L
  if(!LPL){
    Y_obs <- matrix(nrow = 0, ncol = 0)
  }else{
    if(!any(is.na(LPL_VoI))){
      LPL_subset <- TRUE
      if(all(is.character(LPL_VoI))){
        VoI <- which(variables %in% LPL_VoI)
        if(length(LPL_VoI) != length(VoI)){
          stop("Cannot find variables of interest specified in 'LPL_VoI' in the data! \n")
        }
      }else if(all(is.numeric(LPL_VoI))){
        VoI <- as.integer(LPL_VoI)
      }
    }else{
      LPL_subset <- FALSE
    }
    Y_obs <- matrix(Y_obs, length(ahead), M)
    colnames(Y_obs) <- variables
  }

  ret <- out_of_sample(each,
                       X_T_plus_1,
                       object$PHI,
                       object$U,
                       object$facload,
                       sv_h_T,
                       ahead,
                       sv_mu,
                       sv_phi,
                       sv_sigma,
                       as.integer(sv_indicator-1),
                       factor,
                       LPL,
                       Y_obs,
                       LPL_subset,
                       as.integer(VoI-1),
                       simulate_predictive)


  if(LPL){
    numericalnormalizer <- apply(ret$LPL_draws,1,max) - 700
    LPL <- log(rowMeans(exp( ret$LPL_draws - numericalnormalizer))) + numericalnormalizer
    names(LPL) <- paste0("t+", ahead)
    ret$LPL <- LPL

    ret$LPL_univariate <- log(apply(ret$PL_univariate_draws, 1:2, mean))
    rownames(ret$LPL_univariate) <- paste0("t+", ahead)
    colnames(ret$LPL_univariate) <- variables
    if(LPL_subset){
      if(any(is.na(ret$LPL_sub_draws))){
        nrnas <- apply(ret$LPL_sub_draws,1,function(x) length(which(is.na(x))))
        warning(paste0(paste0("For t+", ahead," ", nrnas, " draws had to be discarded for the computation of LPL_Voi due to numerical issues!\n")))
        LPL_VoI <- rep_len(NA, length(ahead))
        for(e in seq_along(ahead)){
          numericalnormalizer2 <- max(ret$LPL_sub_draws[e,], na.rm = TRUE) - 700
          LPL_VoI[e] <- log(mean(exp(ret$LPL_sub_draws[e,] - numericalnormalizer2), na.rm=TRUE)) + numericalnormalizer2
        }
      }else{
        numericalnormalizer2 <- apply(ret$LPL_sub_draws,1,max) - 700
        LPL_VoI <- log(rowMeans(exp( ret$LPL_sub_draws - numericalnormalizer2))) +
          numericalnormalizer2
      }
      names(LPL_VoI) <- paste0("t+", ahead)
      ret$LPL_VoI <- LPL_VoI

      #ret$LPL_VoI_draws <- ret$LPL_sub_draws
      ret$VoI <- variables[VoI]
    }
  }else if(!LPL){
    ret$LPL_draws <- NULL
    ret$PL_univariate_draws <- NULL
    ret$LPL_sub_draws <- NULL
  }
  if(simulate_predictive){
    dimnames(ret$predictions) <- list(paste0("t+", ahead), variables, NULL)
  }
  ret$Y_obs <- Y_obs
  ret$Yraw <- object$Yraw

  class(ret) <- "bayesianVARs_predict"
  ret
}

predict_old <- function(object, n.ahead, stable = TRUE, LPL = FALSE, Y_obs = NA, LPL_VoI = NA, ...){

  if(stable){
    cat("\nArgument 'stable' is 'TRUE'. Therfore, only stable draws are considered for prediction.\n")
    object <- stable_bvar(object, quiet = TRUE)
    cat("\n",dim(object$PHI)[1],"stable posterior draws remaining for prediction!\n")
  }

  # relevant mod settings
  SV <- object$heteroscedastic
  sv_indicator <- which(SV==TRUE)
  intercept <- object$intercept

  # data preparation
  variables <- colnames(object$Y)
  draws <- dim(object$PHI)[3]
  lags <- object$lags # number of lags
  M <- ncol(object$Y) # dimensionality, number of time series
  # if(object$sigma_type == "factor"){
  #   r <- dim(object$facload)[3] # number of factors
  # }
  if(LPL){
    if(!any(is.na(LPL_VoI))){
      LPL_subset <- TRUE
      if(all(is.character(LPL_VoI))){
        VoI <- which(variables %in% LPL_VoI)
        if(length(LPL_VoI) != length(VoI)){
          stop("Cannot find variables of interest specified in 'LPL_VoI' in the data! \n")
        }
      }else if(all(is.numeric(LPL_VoI))){
        VoI <- LPL_VoI
      }
    }else{
      LPL_subset <- FALSE
    }
    Y_obs <- matrix(Y_obs, n.ahead, M)
    colnames(Y_obs) <- variables
  }
  if(any(SV)) {
    # extract sv parameters
    sv_mu <- (object$sv_para[1,,])
    sv_phi <- (object$sv_para[2,,])
    sv_sigma <- (object$sv_para[3,,])

  }
  # extract current state of logvariance (in case of homoscedasticity one could pick any, since constant...)
  sv_h_T <- object$logvar[dim(object$logvar)[1],,]

  # storage
  predictions <- array(as.numeric(NA), c(draws, n.ahead, M),
                       dimnames = list(NULL, paste0("t+", 1:n.ahead), variables))
  if(LPL){
    LPL_draws <- matrix(as.numeric(NA), draws, n.ahead)
    colnames(LPL_draws) <- paste0("t+", 1:n.ahead)
    PL_univariate_draws <- array(as.numeric(NA), c(draws, n.ahead, M),
                                 dimnames = list(NULL, paste0("t+", 1:n.ahead), variables))
    if(LPL_subset){
      LPL_sub_draws <- matrix(as.numeric(NA), draws, n.ahead)
      colnames(LPL_sub_draws) <- paste0("t+", 1:n.ahead)
    }
  }

  ## X_fore1: predictors for one-step ahead forecasts
  X_fore1 <- as.numeric(object$datamat[nrow(object$datamat), 1:(lags*M)])

  if(intercept) X_fore1 <- c(X_fore1, 1)

  U <- diag(1, nrow = M)
  upperind <- which(upper.tri(U))

  for (i in seq.int(draws)) {

    X_fore_k <- X_fore1

    # initialize latent logvariance at current state
    h_fore <- sv_h_T[, i]

    if(object$sigma_type == "cholesky"){
      U[upperind] <- object$U[,i]
      U_inv <- backsolve(U, diag(M))
    }

    for(k in seq.int(n.ahead)){

      mean_fore <- as.vector(X_fore_k%*%object$PHI[,,i])

      # compute prediction of variance-covariance matrix
      if(any(SV)){ # in case of SV, predict logvariances
        # compute k-step ahead forecast of latent log-volas
        h_fore[sv_indicator] <- sv_mu[sv_indicator, i] + sv_phi[sv_indicator, i]*(h_fore[sv_indicator] - sv_mu[sv_indicator, i]) +
          sv_sigma[sv_indicator, i]*stats::rnorm(length(sv_indicator), mean = 0, sd = 1)
      }

      if(object$sigma_type == "factor"){
        Sigma_fore <- crossprod(t(object$facload[,,i])*exp(h_fore[-c(1:M)]/2)) # facload %*% diag(exp(h)) %*% t(facload)
        diag(Sigma_fore) <- diag(Sigma_fore) +  exp(h_fore[1:M])
      }else if(object$sigma_type == "cholesky"){
        Sigma_chol_fore <- diag(exp(h_fore/2)) %*% U_inv
        Sigma_fore <- crossprod(Sigma_chol_fore)
      }

      predictions[i,k,] <-  if(object$sigma_type == "factor"){
        MASS::mvrnorm(1, mean_fore, Sigma_fore)
      }else if(object$sigma_type == "cholesky"){
        tryCatch(
          mean_fore + stats::rnorm(M) %*% Sigma_chol_fore,
          error = function(e) MASS::mvrnorm(1, mean_fore, Sigma_fore)
        )
      }

      if(LPL){
        LPL_draws[i,k] <- if(object$sigma_type == "factor"){
          mvtnorm::dmvnorm(as.vector(Y_obs[k,]),mean_fore,Sigma_fore, log = TRUE)
        }else if(object$sigma_type == "cholesky"){
          mydmvnorm(Y_obs[k,], mean_fore,Sigma_chol_fore, log = TRUE)
        }
        PL_univariate_draws[i, k,] <-  stats::dnorm(as.vector(Y_obs[k, ]), mean_fore, sqrt(diag(Sigma_fore)))
        if(LPL_subset){
          LPL_sub_draws[i, k] <- mvtnorm::dmvnorm(Y_obs[k, VoI], mean_fore[VoI], (Sigma_fore[VoI,VoI, drop = FALSE]), log = TRUE)
        }
      }

      if(k<n.ahead){
        if(lags==1){
          X_fore_k <- predictions[i,k,]
        }else{
          X_fore_k <- c(predictions[i,k,], X_fore_k[1:((lags-1)*M)])
        }

        if(intercept){
          X_fore_k <- c(X_fore_k,1)
        }
      }

    }# end 1:k
  }# end 1:draws

  out <- list(predictions = predictions)

  if(LPL){
    numericalnormalizer <- apply(LPL_draws,2,max) - 700
    LPL <- log(colMeans(exp( t(t(LPL_draws) - numericalnormalizer)))) + numericalnormalizer
    names(LPL) <- paste0("t+", 1:n.ahead)
    out$LPL <- LPL
    out$LPL_draws <- LPL_draws

    out$LPL_univariate <- log(apply(PL_univariate_draws, 2:3, mean))
    rownames(out$LPL_univariate) <- paste0("t+", 1:n.ahead)
    out$PL_univariate_draws <- PL_univariate_draws

    if(LPL_subset){
      if(any(is.na(LPL_sub_draws))){
        nrnas <- apply(LPL_sub_draws,2,function(x) length(which(is.na(x))))
        warning(paste0(paste0("For t+", 1:n.ahead," ", nrnas, " draws had to be discarded for the computation of LPL_Voi due to numerical issues!\n")))
        LPL_VoI <- rep(NA, n.ahead)
        for(e in seq.int(n.ahead)){
          numericalnormalizer2 <- max(LPL_sub_draws[,e], na.rm = TRUE) - 700
          LPL_VoI[e] <- log(mean(exp(LPL_sub_draws[,e] - numericalnormalizer2), na.rm=TRUE)) + numericalnormalizer2
        }
      }else{
        numericalnormalizer2 <- apply(LPL_sub_draws,2,max) - 700
        LPL_VoI <- log(colMeans(exp( t(t(LPL_sub_draws) - numericalnormalizer2)))) +
          numericalnormalizer2
      }
      names(LPL_VoI) <- paste0("t+", 1:n.ahead)
      out$LPL_VoI <- LPL_VoI
      out$LPL_VoI_draws <- LPL_sub_draws
      out$VoI <- variables[VoI]
    }
  }
  #class(out) <- "bayesianVARs_predict"
  return(out)

}

#' Print method for bayesianVARs_predict objects
#'
#' Print method for bayesianVARs_predict objects.
#'
#' @param x A `bayesianVARs_predict` object obtained via
#'   [`predict.bayesianVARs_bvar()`].
#' @param ... Currently ignored!
#'
#' @return Returns `x` invisibly.
#' @export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Split data in train and test
#' train <- data[1:(nrow(data)-4),]
#' test <- data[-c(1:(nrow(data)-4)),]
#'
#' # Estimate model using train data only
#' mod <- bvar(train, quiet = TRUE)
#'
#' # Simulate from 1-step ahead posterior predictive
#' predictions <- predict(mod, ahead = 1L)
#' print(predictions)
print.bayesianVARs_predict <- function(x, ...){
  cat(paste("\nGeneric functions for bayesianVARs_predict objects:\n",
            " - summary.bayesianVARs_predict(),\n",
            " - pairs.bayesianVARs_predict(),\n",
            " - plot.bayesianVARs_predict() (alias for pairs.bayesianVARs_predict()).\n"))
  invisible(x)
}

#' Summary method for bayesianVARs_predict objects
#'
#' Summary method for `bayesianVARs_predict` objects.
#'
#' @param object A `bayesianVARs_predict` object obtained via
#'   [`predict.bayesianVARs_bvar()`].
#' @param ... Currently ignored!
#'
#' @return A `summary.bayesianVARs_predict` object.
#' @export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Split data in train and test
#' train <- data[1:(nrow(data)-4),]
#' test <- data[-c(1:(nrow(data)-4)),]
#'
#' # Estimate model using train data only
#' mod <- bvar(train, quiet = TRUE)
#'
#' # Simulate from 1-step ahead posterior predictive
#' predictions <- predict(mod, ahead = 1L)
#' summary(predictions)
summary.bayesianVARs_predict <- function(object, ...){
  out <- list()
  if(!is.null(object$LPL)){
    out$LPL <- object$LPL
    out$LPL_univariate <- object$LPL_univariate
  }
  if(!is.null(object$LPL_VoI)){
    out$LPL_VoI <- object$LPL_VoI
    out$VoI <- object$VoI
  }

  if(!is.null(object$predictions)){
    out$prediction_quantiles <- apply(object$predictions, 1:2, stats::quantile, c(.05,.5,.95))
  }

  class(out) <- "summary.bayesianVARs_predict"
  out
}


#' Print method for summary.bayesianVARs_predict objects
#'
#' Print method for `summary.bayesianVARs_predict` objects.
#'
#' @param x A `summary.bayesianVARs_predict` object obtained via [`summary.bayesianVARs_predict()`].
#' @param ... Currently ignored.
#'
#' @return Returns `x` invisibly.
#' @export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Split data in train and test
#' train <- data[1:(nrow(data)-4),]
#' test <- data[-c(1:(nrow(data)-4)),]
#'
#' # Estimate model using train data only
#' mod <- bvar(train, quiet = TRUE)
#'
#' # Simulate from 1-step ahead posterior predictive
#' predictions <- predict(mod, ahead = 1L)
#' sum <- summary(predictions)
#' print(sum)
#'
print.summary.bayesianVARs_predict <- function(x, ...){
  digits <- max(3, getOption("digits") - 3)
  if(!is.null(x$LPL)){
    cat("\nLPL:\n")
    print(x$LPL, digits = digits)

    if(!is.null(x$LPL_VoI)){
      n <- length(x$VoI)
      cat("\nMarginal joint LPL of ")
      cat(x$VoI[1:(n-1)] ,sep = ", ")
      cat(" & ", x$VoI[n], ":\n", sep = "")
      print(x$LPL_VoI, digits = digits)
    }

    cat("\nMarginal univariate LPLs:\n")
    print(x$LPL_univariate, digits = digits)

  }

  if(!is.null(x$prediction_quantiles)){
    cat("\nPrediction quantiles:\n")
    print(x$prediction_quantiles, digits = digits)
  }

  invisible(x)
}
