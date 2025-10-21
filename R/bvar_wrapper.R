#'Markov Chain Monte Carlo Sampling for Bayesian Vectorautoregressions
#'
#'\code{bvar} simulates from the joint posterior distribution of the parameters
#'and latent variables and returns the posterior draws.
#'
#'The VAR(p) model is of the following form: \eqn{ \boldsymbol{y}^\prime_t = \boldsymbol{\iota}^\prime +
#'\boldsymbol{x}^\prime_t\boldsymbol{\Phi} + \boldsymbol{\epsilon}^\prime_t}, where
#'\eqn{\boldsymbol{y}_t} is a \eqn{M}-dimensional vector of dependent variables
#'and \eqn{\boldsymbol{\epsilon}_t} is the error term of the same dimension.
#'\eqn{\boldsymbol{x}_t} is a \eqn{K=pM}-dimensional vector containing
#'lagged/past values of the dependent variables \eqn{\boldsymbol{y}_{t-l}} for
#'\eqn{l=1,\dots,p} and \eqn{\boldsymbol{\iota}} is a constant term (intercept)
#'of dimension \eqn{M\times 1}. The reduced-form coefficient matrix
#'\eqn{\boldsymbol{\Phi}} is of dimension \eqn{K \times M}.
#'
#'\code{bvar} offers two different specifications for the errors: The user can
#'choose between a factor stochastic volatility structure or a cholesky
#'stochastic volatility structure. In both cases the disturbances
#'\eqn{\boldsymbol{\epsilon}_t} are assumed to follow a \eqn{M}-dimensional
#'multivariate normal distribution with zero mean and variance-covariance matrix
#'\eqn{\boldsymbol{\Sigma}_t}. In case of the
#'cholesky specification \eqn{\boldsymbol{\Sigma}_t = \boldsymbol{U}^{\prime -1} \boldsymbol{D}_t
#'\boldsymbol{U}^{-1}}, where \eqn{\boldsymbol{U}^{-1}} is upper unitriangular (with ones on
#'the diagonal). The diagonal matrix \eqn{\boldsymbol{D}_t} depends upon latent
#'log-variances, i.e. \eqn{\boldsymbol{D}_t=diag(exp(h_{1t}),\dots,
#'exp(h_{Mt})}. The log-variances follow a priori independent autoregressive
#'processes \eqn{h_{it}\sim N(\mu_i + \phi_i(h_{i,t-1}-\mu_i),\sigma_i^2)} for
#'\eqn{i=1,\dots,M}. In case of the factor structure,
#'\eqn{\boldsymbol{\Sigma}_t = \boldsymbol{\Lambda} \boldsymbol{V}_t \boldsymbol{\Lambda}^\prime +
#'\boldsymbol{G}_t}. The diagonal matrices \eqn{\boldsymbol{V}_t} and
#'\eqn{\boldsymbol{G}_t} depend upon latent log-variances, i.e.
#'\eqn{\boldsymbol{G}_t=diag(exp(h_{1t}),\dots, exp(h_{Mt})} and
#'\eqn{\boldsymbol{V}_t=diag(exp(h_{M+1,t}),\dots, exp(h_{M+r,t})}. The
#'log-variances
#'follow a priori independent autoregressive processes \eqn{h_{it}\sim N(\mu_i +
#'\phi_i(h_{i,t-1}-\mu_i),\sigma_i^2)} for \eqn{i=1,\dots,M} and
#'\eqn{h_{M+j,t}\sim N(\phi_ih_{M+j,t-1},\sigma_{M+j}^2)} for \eqn{j=1,\dots,r}.
#'
#'@param data Data matrix (can be a time series object). Each of \eqn{M} columns
#'  is assumed to contain a single time-series of length \eqn{T}.
#'@param lags Integer indicating the order of the VAR, i.e. the number of lags
#'  of the dependent variables included as predictors.
#'@param draws single integer indicating the number of draws after the burnin
#'@param burnin single integer indicating the number of draws discarded as
#'  burnin
#'@param thin single integer. Every \eqn{thin}th draw will be stored. Default is
#'  \code{thin=1L}.
#'@param prior_intercept Either \code{prior_intercept=FALSE} and no constant
#'  term (intercept) will be included. Or a numeric vector of length \eqn{M}
#'  indicating the (fixed) prior standard deviations on the constant term. A
#'  single number will be recycled accordingly. Default is
#'  \code{prior_intercept=10}.
#'@param prior_phi \code{bayesianVARs_prior_phi} object specifying prior for the
#'  reduced form VAR coefficients. Best use constructor
#'  \code{\link{specify_prior_phi}}.
#'@param prior_sigma \code{bayesianVARs_prior_sigma} object specifying prior for
#'  the variance-covariance matrix of the VAR. Best use constructor
#'  \code{\link{specify_prior_sigma}}.
#'@param sv_keep String equal to \code{"all"} or \code{"last"}. In case of
#'  \code{sv_keep = "last"}, the default, only draws for the very last
#'  log-variance \eqn{h_T} are stored.
#'@param quiet logical value indicating whether information about the progress
#'  during sampling should be displayed during sampling (default is
#'  \code{TRUE}).
#'@param startvals optional list with starting values.
#'@param expert optional list with expert settings.
#'
#'@section MCMC algorithm: To sample efficiently the reduced-form VAR
#'  coefficients assuming a **factor structure for the errors**, the equation
#'  per equation algorithm in Kastner & Huber (2020) is implemented. All
#'  parameters and latent variables associated with the factor-structure are
#'  sampled using package \code{\link[factorstochvol]{factorstochvol-package}}'s
#'  function `update_fsv` callable on the C-level only.
#'
#'  To sample efficiently the reduced-form VAR coefficients, assuming a
#'  **cholesky-structure for the errors**, the corrected triangular algorithm in
#'  Carriero et al. (2021) is implemented. The SV parameters and latent
#'  variables are sampled using package \code{\link[stochvol]{stochvol}}'s
#'  \code{\link[stochvol]{update_fast_sv}} function. The precision parameters,
#'  i.e. the free off-diagonal elements in \eqn{\boldsymbol{U}}, are sampled as
#'  in Cogley and Sargent (2005).
#'
#'@references  Gruber, L. and Kastner, G. (2025). Forecasting macroeconomic data
#'  with Bayesian VARs: Sparse or dense? It depends! \emph{International Journal
#'  of Forecasting}. \doi{10.1016/j.ijforecast.2025.02.001}.
#'
#'@references  Kastner, G. and Huber, F. Sparse (2020). Bayesian vector
#'  autoregressions in huge dimensions. \emph{Journal of Forecasting}.
#'  \bold{39}, 1142--1165, \doi{10.1002/for.2680}.
#'
#'@references Kastner, G. (2019). Sparse Bayesian Time-Varying Covariance
#'  Estimation in Many Dimensions \emph{Journal of Econometrics}, \bold{210}(1),
#'  98--115, \doi{10.1016/j.jeconom.2018.11.007}.
#'
#'@references Carriero, A. and Chan, J. and  Clark, T. E. and Marcellino, M.
#'  (2021). Corrigendum to “Large Bayesian vector autoregressions with
#'  stochastic volatility and non-conjugate priors” \[J. Econometrics 212 (1)
#'  (2019) 137–154\]. \emph{Journal of Econometrics},
#'  \doi{10.1016/j.jeconom.2021.11.010}.
#'
#'@references Cogley, S. and Sargent, T. (2005). Drifts and volatilities:
#'  monetary policies and outcomes in the post WWII US. \emph{Review of Economic
#'  Dynamics}, \bold{8}, 262--302, \doi{10.1016/j.red.2004.10.009}.
#'
#'@references Hosszejni, D. and Kastner, G. (2021). Modeling Univariate and
#'  Multivariate Stochastic Volatility in R with stochvol and factorstochvol.
#'  \emph{Journal of Statistical Software}, \emph{100}, 1–-34.
#'  \doi{10.18637/jss.v100.i12}.
#'
#'@seealso
#'* Helpers for prior configuration: [specify_prior_phi()], [specify_prior_sigma()].
#'* Plotting: [`plot.bayesianVARs_bvar()`].
#'* Extractors: [`coef.bayesianVARs_bvar()`], [`vcov.bayesianVARs_bvar()`].
#'* 'stable' bvar: [stable_bvar()].
#'* summary method: [`summary.bayesianVARs_bvar()`].
#'* predict method: [`predict.bayesianVARs_bvar()`].
#'* fitted method: [`fitted.bayesianVARs_bvar()`].
#'
#'@return An object of type `bayesianVARs_bvar`, a list containing the following
#'  objects:
#'* `PHI`: A `bayesianVARs_coef` object, an array, containing the posterior draws
#'  of the VAR coefficients (including the intercept).
#'* `U`: A `bayesianVARs_draws` object, a matrix, containing the posterior draws
#'  of the contemporaneous coefficients (if cholesky decomposition for sigma is
#'  specified).
#'* `logvar`: A `bayesianVARs_draws` object containing the log-variance draws.
#'* `sv_para`: A `baysesianVARs_draws` object containing the posterior draws of
#'  the stochastic volatility related parameters.
#'* `phi_hyperparameter`: A matrix containing the posterior draws of the
#'  hyperparameters of the conditional normal prior on the VAR coefficients.
#'* `u_hyperparameter`: A matrix containing the posterior draws of the
#'  hyperparameters of the conditional normal prior on U (if cholesky
#'  decomposition for sigma is specified).
#'* `bench`: Numerical indicating the average time it took to generate one
#'  single draw of the joint posterior distribution of all parameters.
#'* `V_prior`: An array containing the posterior draws of the variances of the
#'  conditional normal prior on the VAR coefficients.
#'* `facload`: A `bayesianVARs_draws` object, an array, containing draws from the
#'  posterior distribution of the factor loadings matrix (if factor
#'  decomposition for sigma is specified).
#'* `fac`: A `bayesianVARs_draws` object, an array, containing factor draws from
#'  the posterior distribution (if factor decomposition for sigma is specified).
#'* `Y`: Matrix containing the dependent variables used for estimation.
#'* `X` matrix containing the lagged values of the dependent variables, i.e.
#'  the covariates.
#'* `lags`: Integer indicating the lag order of the VAR.
#'* `intercept`: Logical indicating whether a constant term is included.
#'* `heteroscedastic` logical indicating whether heteroscedasticity is assumed.
#'* `Yraw`: Matrix containing the dependent variables, including the initial
#'  'lags' observations.
#'* `Traw`: Integer indicating the total number of observations.
#'* `sigma_type`: Character specifying the decomposition of the
#'  variance-covariance matrix.
#'* `datamat`: Matrix containing both 'Y' and 'X'.
#'* `config`: List containing information on configuration parameters.
#'@export
#'
#'@examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Estimate a model
#' mod <- bvar(data, sv_keep = "all", quiet = TRUE)
#'
#' # Plot
#' plot(mod)
#'
#' # Summary
#' summary(mod)
#'
bvar <- function(data,
                 lags=1L,
                 draws=1000L,
                 burnin=1000L,
                 thin = 1L,
                 prior_intercept = 10,
                 prior_phi = specify_prior_phi(data = data, lags = lags, prior = "HS"),
                 prior_sigma = specify_prior_sigma(data = data, type = "factor", quiet = TRUE),
                 sv_keep = "last",
                 quiet=FALSE,
                 startvals = list(),
                 expert = list()
                      ){


# Data preliminaries ------------------------------------------------------

  # M: number of variables,
  # Tobs: number of observations used for estimation,
  # K: number of covariates per equation (without intercepts!!!),
  # n: number of VAR coefficients (without intercepts!!!),
  # Y: data without first lags observations,
  # X: lagged values of data

  thin <- as.integer(thin)
  if(thin<1) stop("Argument 'thin' must be an integer greater than 0!")
  draws <- as.integer(draws)
  burnin <- as.integer(burnin)
  if(!sv_keep%in%c("last", "all")){
    stop("Argument 'sv_keep' must be one of 'last' or 'all'.")
  }

  lags <- as.integer(lags)
  M <- ncol(data)
  K <- lags*M
  Traw <- nrow(data)
  Y_tmp <- as.matrix(data)
  if (any(is.na(Y_tmp))){
    stop("\nNAs in data.\n")
  }
  if (ncol(Y_tmp) < 2) {
    stop("The matrix 'data' should contain at least two variables. \n")
  }
  if (is.null(colnames(Y_tmp))) {
    colnames(Y_tmp) <- paste("y", 1:ncol(Y_tmp), sep = "")
    warning(paste("No column names supplied in data, using:",
                  paste(colnames(Y_tmp), collapse = ", "), ", instead.\n"))
  }
  colnames(Y_tmp) <- make.names(colnames(Y_tmp))
  # embed: get lagged values
  X <- embed(Y_tmp, dimension = lags + 1)[, -(1:M)]
  colnames(X) <- paste0(colnames(Y_tmp), ".l", sort(rep(1:lags,M)))
  if(is.numeric(prior_intercept)){
    X <- cbind(X,1)
    colnames(X)[ncol(X)] <- c("intercept")
  }
  Y <- Y_tmp[-c(1:lags), ]

  Tobs <- Traw - lags
  n <- K*M

# Input checks ------------------------------------------------------------

  if(!inherits(prior_phi, "bayesianVARs_prior_phi")){
    stop("\nArgument 'prior_phi' must be of class 'bayesianVARs_prior_phi'. Please use helper function 'specify_prior_phi()'!\n")
  }
  if(prior_phi[["general_settings"]][["M"]] != M) stop(paste0("\nObject 'prior_phi' is specified for ", prior_phi[["general_settings"]][["M"]], "time series. 'data', however, consists of ", M, "time series!\n"))
  if(prior_phi[["general_settings"]][["lags"]] != lags) stop("\nLag-length 'lags' does not coincide with lag-length specified in object 'prior_phi'!\n")

  if(!inherits(prior_sigma, "bayesianVARs_prior_sigma")){
    stop("\nArgument 'prior_sigma' must be of class 'bayesianVARs_prior_sigma'. Please use helper function 'specify_prior_sigma()'!\n")
  }
  if(prior_sigma[["general_settings"]][["M"]] != M){
    stop(paste0("\n'prior_sigma' is specified for ", prior_sigma[["general_settings"]][["M"]], " timeseries. 'data', however, consits of ",M," timeseries!\n"))
  }

# Indicator matrix --------------------------------------------------------

  if(is.numeric(prior_intercept)){
    i_intercept <- rep(0,M)
  }else i_intercept <- NULL

  # indicator for for coefficients of lagged ys.
  i_mat <- prior_phi[["general_settings"]][["i_mat"]]
  # add indicator for intercept
  i_mat <- rbind(i_mat, i_intercept)
  mode(i_mat) <- "integer"
  i_vec <- as.vector(i_mat)


# PHI settings ------------------------------------------------------------

  if(is.logical(prior_intercept)){
    if(any(prior_intercept == TRUE)){
      stop(paste0("\nArgument 'prior_intercept' must be one of 'FALSE', a single numeric greater than 0, or a vector of length M.\n"))
    }
    intercept <- 0L
    priorIntercept <- vector("numeric")
  }else{
    prior_intercept <- prior_intercept^2 # transform to variances
    if(length(prior_intercept)==1L){
      priorIntercept <- rep(prior_intercept, M)
    }else if(length(prior_intercept)==M){
      priorIntercept <- prior_intercept
    }else{
      stop(paste0("\nArgument 'prior_intercept' must be one of 'FALSE', a single numeric greater than 0, or a vector of length ", M, ".\n"))
    }

    intercept <- 1L
  }
  PHI0 <- matrix(0, K+intercept, M) # prior mean of intercept is 0
  PHI0[1:K,] <- prior_phi[["general_settings"]][["PHI0"]]

  # Initialize PHI --------

  if(prior_phi[["prior_phi_cpp"]][["prior"]] == "SSVS" | !exists("PHI", startvals) | (!exists("U", startvals) & prior_sigma[["prior_sigma_cpp"]][["type"]] == "cholesky")){
    # Note: if startvals are provided and if(!prior_phi[["general_settings"]][["SSVS_semiautomatic"]]) then
    # the following computations are redundant

    # Posterior mean of a flat conjugate Normal inverse Wishart prior
    # exists even when OLS estimate does not exist (in situations where Tobs < K)
    # N(0, 10^3) on PHI, and invWish(I, M+2) on Sigma
    XX <- crossprod(X)
    V_post_flat <- chol2inv(chol(diag(1/rep(10^3, (K+intercept))) + XX))
    PHI_flat <- V_post_flat %*% (diag(1/rep(10^3, (K+intercept)))%*%PHI0 + t(X)%*%Y)
    S_post <- diag(M) + crossprod(Y - X%*%PHI_flat) + t(PHI_flat - PHI0) %*%
      diag(1/rep(10^3, (K+intercept))) %*% (PHI_flat - PHI0)
    Sigma_flat <- (S_post)/(M +2 + Tobs - M - 1)
    U <- chol(Sigma_flat)
    D <- diag(U)^2
    L_inv <- U/sqrt(D)
    L_flat <- backsolve(L_inv, diag(M))
  }

  # "empirical" prior settings
  if(prior_phi[["prior_phi_cpp"]][["prior"]] == "SSVS"){

    if(prior_phi[["general_settings"]][["SSVS_semiautomatic"]]) {

      # standard errors of flat phi posterior estimates
      sigma_phi <- sqrt(diag(Sigma_flat %x% V_post_flat))
      # select all, except those of intercept
      sigma_phi <- sigma_phi[which(i_vec!=0)]
      if(length(sigma_phi)!=lags*M^2){
        stop(paste0(length(sigma_phi), " length(sigma_phi)!= lags*M^2"))
      }

      # scale tau with variances
      prior_phi[["prior_phi_cpp"]][["SSVS_tau0"]] <- prior_phi[["prior_phi_cpp"]][["SSVS_tau0"]]*sigma_phi
      prior_phi[["prior_phi_cpp"]][["SSVS_tau1"]] <- prior_phi[["prior_phi_cpp"]][["SSVS_tau1"]]*sigma_phi
    }

  }else if(prior_phi[["prior_phi_cpp"]][["prior"]] == "HMP"){
    sigma_sq <- MP_sigma_sq(data, 6)
    # prepare prior variances down to lambdas
    prior_phi[["prior_phi_cpp"]][["V_i_prep"]] <- MP_V_prior_prep(sigma_sq, (K+intercept), M, intercept>0)
  }

# Sigma -------------------------------------------------------------------

  # everything is done within specify_prior_sigma()

# startvals ---------------------------------------------------------------

  # Note to myself: for prior_sigma$prior_sigma_cpp$type=="factor", result does not depend on starting value of PHI
  startvals[["PHI"]] <- if(exists("PHI", startvals)) startvals[["PHI"]] else PHI_flat
  startvals[["U"]] <- if(!prior_sigma[["prior_sigma_cpp"]][["type"]]=="cholesky") matrix(0,1,1) else if(exists("U", startvals)) startvals[["U"]] else L_flat
  startvals[["sv_para"]] <- matrix(0,1,1)
  startvals[["sv_logvar"]] <- matrix(0,1,1)
  startvals[["sv_logvar0"]] <- numeric(1L)
  startvals[["factor_startval"]] <- list(facload = matrix(0,1,1), fac = matrix(0,1,1), tau2 = matrix(0,1,1))

  if(prior_sigma[["prior_sigma_cpp"]][["type"]] == "cholesky"){
    startvals[["sv_para"]] <- matrix(data= c(rep(-10,M), rep(0.9,M), rep(0.2,M)), #, rep(-10,M)
                                         nrow = 3, ncol = M, byrow = TRUE)
    startvals[["sv_logvar0"]] <- rep(-10,M)
    startvals[["sv_logvar"]] <- matrix(rep(-10, Tobs*M), Tobs,M)
  }else if(prior_sigma[["prior_sigma_cpp"]][["type"]] == "factor"){
    startfacload <- matrix(rnorm(M*prior_sigma[["prior_sigma_cpp"]][["factor_factors"]], sd = .5)^2, nrow=M, ncol=prior_sigma[["prior_sigma_cpp"]][["factor_factors"]])
    startfac <- matrix(rnorm(prior_sigma[["prior_sigma_cpp"]][["factor_factors"]]*Tobs, 0, sd=.1), nrow=prior_sigma[["prior_sigma_cpp"]][["factor_factors"]])
    startpara <- rbind(mu = c(rep(-3, M) + rnorm(M), rep(0, prior_sigma[["prior_sigma_cpp"]][["factor_factors"]])),
                      phi = c(rep(.8, M), rep(.8, prior_sigma[["prior_sigma_cpp"]][["factor_factors"]])) + pmin(rnorm(M + prior_sigma[["prior_sigma_cpp"]][["factor_factors"]], sd=.06), .095),
                      sigma = rep(.1, M + prior_sigma[["prior_sigma_cpp"]][["factor_factors"]]) + rgamma(M + prior_sigma[["prior_sigma_cpp"]][["factor_factors"]], 1, 10))
    startlogvar <- matrix(startpara["mu",][1] + rnorm(Tobs*(M + prior_sigma[["prior_sigma_cpp"]][["factor_factors"]])), Tobs, M + prior_sigma[["prior_sigma_cpp"]][["factor_factors"]])
    startlogvar[,M+which(prior_sigma[["prior_sigma_cpp"]][["sv_heteroscedastic"]][-c(1:M)]==FALSE)] <- 0 # !!! important, needed for factorstochvol: if factor is assumed to be homoscedastic, the corresponding column in logvar has to be 0!!!
    startlogvar0 <- startpara["mu",][1] + rnorm(M + prior_sigma[["prior_sigma_cpp"]][["factor_factors"]])
    starttau2 <- if(!prior_sigma[["prior_sigma_cpp"]][["factor_ngprior"]]){ # if prior is 'normal'
      prior_sigma[["general_settings"]][["factor_starttau2"]]
    }else{
      matrix(1, nrow = M, ncol = prior_sigma[["prior_sigma_cpp"]][["factor_factors"]])
    }
    startvals[["factor_startval"]] <- list(facload = startfacload,
                                          fac = startfac,
                                          tau2 = starttau2)
    startvals[["sv_para"]] <- startpara
    startvals[["sv_logvar0"]] <- startlogvar0
    startvals[["sv_logvar"]] <- startlogvar
  }

  ######### expert settings
  if(!is.null(expert)){
    if(exists("huge", expert)){
      expert[["huge"]] <- expert[["huge"]]
    }else{
      expert[["huge"]] <- FALSE
    }
  }

  if(!quiet){
    cat("\nCalling MCMC sampler for a ", prior_sigma[["prior_sigma_cpp"]][["type"]], "-VAR(",lags,") for ",M, " series of length ", Traw, ".\n", sep = "")
  }

  res <- bvar_cpp(Y,
                  X,
                  M,
                  Tobs,
                  K,
                  draws,
                  burnin,
                  thin,
                  sv_keep,
                  intercept,
                  priorIntercept,
                  PHI0,
                  prior_phi[["prior_phi_cpp"]],
                  prior_sigma[["prior_sigma_cpp"]],
                  startvals,
                  i_mat,
                  i_vec,
                  !quiet,
                  prior_phi[["general_settings"]][["PHI_tol"]],
                  prior_sigma[["general_settings"]][["cholesky_U_tol"]],
                  expert[["huge"]]
                  )

  #Rcpp timer is in nanoseconds
  #conversion to secs per iteration
  res[["bench"]] <- diff(res[["bench"]])/(10^(9)*(draws+burnin))
  attributes(res[["bench"]]) <- list("names" = "secs/itr")
  dimnames(res[["PHI"]])[1] <- list(colnames(X))
  dimnames(res[["PHI"]])[2] <- list(colnames(Y))
  phinames <- as.vector((vapply(seq_len(M), function(i) paste0(colnames(Y)[i], "~", colnames(X[,1:(ncol(X)-intercept)])), character(K))))
  if(prior_phi[["prior_phi_cpp"]][["prior"]] %in% c("DL","DL_h", "GT")){

    rownames(res[["phi_hyperparameter"]]) <- c(paste0("a",1:prior_phi[["prior_phi_cpp"]][["n_groups"]]),
                                          paste0("xi",1:prior_phi[["prior_phi_cpp"]][["n_groups"]]),
                                          paste0("lambda: ", phinames),
                                          paste0("psi: ", phinames))

  }else if(prior_phi[["prior_phi_cpp"]][["prior"]] == "HS"){

    rownames(res[["phi_hyperparameter"]]) <- c(paste0("zeta", 1:prior_phi[["prior_phi_cpp"]][["n_groups"]]),
                                          paste0("varpi", 1:prior_phi[["prior_phi_cpp"]][["n_groups"]]),
                                          paste0("theta: ", phinames),
                                          paste0("nu: ", phinames))
  }else if(prior_phi[["prior_phi_cpp"]][["prior"]] == "SSVS"){

    rownames(res[["phi_hyperparameter"]]) <- c(paste0("gamma: ", phinames), paste0("p_i: ", phinames))#

  }else if(prior_phi[["prior_phi_cpp"]][["prior"]] == "HMP"){

    rownames(res[["phi_hyperparameter"]]) <- c("lambda_1", "lambda_2")

  }

  if(prior_sigma[["prior_sigma_cpp"]][["type"]]=="cholesky"){
    lnames <- rep(as.character(NA), (M^2-M)/2)
    ii <- 1
    for(j in 2:M){
      #lnames <- c(lnames, paste0(colnames(Y)[j],"~", colnames(Y)[1:(j-1)]))
      lnames[ii:(ii+j-2)] <- paste0(colnames(Y)[j],"~", colnames(Y)[1:(j-1)])
      ii <- ii+j-1
    }
    rownames(res[["U"]]) <- lnames
    if(prior_sigma[["prior_sigma_cpp"]][["cholesky_U_prior"]] %in% c("DL", "GT")){

      rownames(res[["u_hyperparameter"]]) <- c("a",
                                          "xi",
                                          paste0("lambda: ", lnames),
                                          paste0("psi: ", lnames))

    }else if(prior_sigma[["prior_sigma_cpp"]][["cholesky_U_prior"]] == "HS"){
      rownames(res[["u_hyperparameter"]]) <- c("zeta",
                                          "varpi",
                                          paste0("theta: ", lnames),
                                          paste0("nu: ", lnames))
    }else if(prior_sigma[["prior_sigma_cpp"]][["cholesky_U_prior"]] == "SSVS"){
      rownames(res[["u_hyperparameter"]]) <- c(paste0("gamma: ", lnames), paste0("p_i: ", lnames))
    }else if(prior_sigma[["prior_sigma_cpp"]][["cholesky_U_prior"]] == "HMP"){
      rownames(res[["u_hyperparameter"]]) <- c("lambda_3")
    }

  }
  res[["fac"]] <- aperm(res[["fac"]], c(2,1,3))
  res[["Y"]] <- Y
  res[["X"]] <- X
  res[["lags"]] <- lags
  res[["intercept"]] <- ifelse(intercept>0, TRUE, FALSE)
  res[["heteroscedastic"]] <- prior_sigma[["prior_sigma_cpp"]][["sv_heteroscedastic"]]
  res[["Yraw"]] <- Y_tmp
  res[["Traw"]] <- Traw
  res[["sigma_type"]] <- prior_sigma[["prior_sigma_cpp"]][["type"]]
  res[["datamat"]] <- data.frame(cbind(Y, X))
  res[["config"]] <- list(draws = draws, burnin = burnin, thin = thin,
                     sv_keep = sv_keep)
  class(res[["PHI"]]) <- c("bayesianVARs_coef", "bayesianVARs_draws")
  class(res[["U"]]) <- "bayesianVARs_draws"
  class(res[["facload"]]) <- "bayesianVARs_draws"
  class(res[["fac"]]) <- "bayesianVARs_draws"
  class(res[["sv_para"]]) <- "bayesianVARs_draws"
  class(res[["logvar"]]) <- "bayesianVARs_draws"
  class(res) <- "bayesianVARs_bvar"
  res
}

#' Extract or Replace Parts of a bayesianVARs_draws object
#'
#' Extract or replace parts of a `bayesianVARs_draws` object.
#'
#' @param x An object of type `bayesianVARs_draws`.
#' @param i indices
#' @param j indices
#' @param ... further indices
#'
#' @return An object of type `bayesianVARs_draws`.
#' @export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Estimate a model
#' mod <- bvar(data, sv_keep = "all", quiet = TRUE)
#'
#' # Extract coefficients, which are of class bayesianVARs_draws
#' phi <- coef(mod)
#' phi[1,1,1]
`[.bayesianVARs_draws` <- function(x, i, j, ...){
  structure(NextMethod(), class = "bayesianVARs_draws")
}

#' Extract or Replace Parts of a bayesianVARs_coef object
#'
#' Extract or replace parts of a `bayesianVARs_coef` object.
#'
#' @param x An object of type `bayesianVARs_coef`.
#' @param i indices
#' @param j indices
#' @param ... further indices
#'
#' @return An object of type `bayesianVARs_coef`.
#' @export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Estimate a model
#' mod <- bvar(data, sv_keep = "all", quiet = TRUE)
#'
#' # Extract coefficients, which are of class bayesianVARs_coef
#' phi <- coef(mod)
#' phi[1,1,1]
`[.bayesianVARs_coef` <- function(x, i, j, ...){
  structure(NextMethod(), class = c("bayesianVARs_coef", "bayesianVARs_draws"))
}

#' Stable posterior draws
#'
#' `stable_bvar()` detects and discards all posterior draws of an
#' \code{bayesianVARs_bvar} object that do not fulfill the stability condition:
#' A VAR(p) model is considered as stable only if the eigenvalues of the
#' companion form matrix lie inside the unit circle.
#'
#' @param object A \code{bayesianVARs_bvar} object obtained via [`bvar()`].
#'
#' @param quiet logical indicating whether informative output should be omitted.
#'
#' @return An object of type `bayesianVARs_bvar`.
#' @export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Estimate a model
#' mod <- bvar(data, sv_keep = "all", quiet = TRUE)
#'
#' # Discard "unstable" draws
#' stable_mod <- stable_bvar(mod)
#'
stable_bvar <- function(object, quiet = FALSE){

  draws <- dim(object[["PHI"]])[3]

  if(!quiet){
    cat("\nOriginal 'bayesianVARs_bvar' object consists of",draws, "posterior draws.\n")
  }

  stable_indicator <- vapply(
    seq.int(draws), FUN = function(i) {
      comp <- get_companion(object[["PHI"]][,,i])
      bool <- if(max(Mod(eigen(comp, only.values = TRUE)[["values"]]))<1) TRUE else FALSE
      bool
    },
    FUN.VALUE = logical(1L)
  )
  if(!quiet){
    cat("\nDetected",sum(!stable_indicator), "unstable draws.\n")
  }

  object[["PHI"]] <- object[["PHI"]][,,stable_indicator, drop = FALSE]

  if(ncol(object[["U"]]) == draws){
    object[["U"]] <- object[["U"]][,stable_indicator, drop = FALSE]
  }
  if(dim(object[["logvar"]])[3] == draws){
    object[["logvar"]] <- object[["logvar"]][,,stable_indicator, drop = FALSE]
  }
  if(dim(object[["sv_para"]])[3] == draws){
    object[["sv_para"]] <- object[["sv_para"]][,,stable_indicator, drop = FALSE]
  }
  if(ncol(object[["phi_hyperparameter"]]) == draws){
    object[["phi_hyperparameter"]] <- object[["phi_hyperparameter"]][,stable_indicator, drop = FALSE]
  }
  if(ncol(object[["u_hyperparameter"]]) == draws){
    object[["u_hyperparameter"]] <- object[["u_hyperparameter"]][,stable_indicator, drop = FALSE]
  }
  if(dim(object[["V_prior"]])[3] == draws){
    object[["V_prior"]] <- object[["V_prior"]][,,stable_indicator, drop = FALSE]
  }
  if(dim(object[["facload"]])[3] == draws){
    object[["facload"]] <- object[["facload"]][,,stable_indicator, drop = FALSE]
  }
  if(dim(object[["fac"]])[3] == draws){
    object[["fac"]] <- object[["fac"]][,,stable_indicator, drop = FALSE]
  }

  if(!quiet){
    cat("\nRemaining draws:",dim(object[["PHI"]])[3], "!\n")
  }

  class(object) <- c(oldClass(object), "bayesianVARs_bvar_stable")
  object
}
