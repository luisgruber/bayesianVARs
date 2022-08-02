#'Markov Chain Monte Carlo Sampling for Bayesian Vectorautoregressions
#'
#'\code{bvar} simulates from the joint posterior distribution of the parameters
#'and latent variables and returns the posterior draws.
#'
#'The VAR(p) model is of the following form:  \eqn{ \bold{y}^\prime_t =
#'\bold{x}^\prime_t\bold{\Phi} + \bold{\epsilon}^\prime_t}, where
#'\eqn{\bold{y}_t} is a \eqn{M \times 1} vector of dependent variables and
#'\eqn{\bold{\epsilon}_t} is the error term of the same dimension.
#'\eqn{\bold{x}_t} is a \eqn{K \times 1} vector containing lagged/past values of
#'the independent variables \eqn{\bold{y}_{t-l}} for \eqn{l=1,\dots,p} and if
#'specified the constant (intercept) term. If the constant term is included
#'\eqn{K=pM+1}, else \eqn{K=pM}. The reduced-form coefficient matrix
#'\eqn{\bold{\Phi}} is of dimension \eqn{K \times M}.
#'
#'The disturbances \eqn{\bold{\epsilon}_t} are assumed to follow a
#'\eqn{M}-dimensional multivariate normal distribution with zero mean and
#'variance-covariance matrix \eqn{\bold{\Sigma}_t}. This matrix can be
#'decomposed as \eqn{\bold{\Sigma}_t = \bold{L}^{\prime -1} \bold{D}_t
#'\bold{L}^{-1}}, where \eqn{\bold{L}^{-1}} is upper unitriangular (with ones on
#'the diagonal) and \eqn{\bold{D}_t} diagonal. If \code{SV=TRUE} the matrix
#'\eqn{\bold{D}_t} is assumed to be time-varying and the orthogonalized
#'disturbances \eqn{\bold{\xi}^\prime_t = \bold{\epsilon}_t^\prime\bold{L}} will
#'be modeled with stochastic volatility specification. By setting
#'\code{SV=FALSE} a constant matrix \eqn{\bold{D}_t=\bold{D}} will be estimated.
#'
#'@section MCMC algorithm: To sample efficiently the reduced form VAR
#'  coefficients, the corrected triangular algorithm in Carriero et al. (2021)
#'  is implemented. The SV parameters and latent variables are sampled using
#'  package \code{\link{stochvol}}'s \code{\link[stochvol]{update_fast_sv}}
#'  function. The precision parameters, i.e. the free off-diagonal elements in
#'  \eqn{\bold{L}}, are sampled as in Cogley and Sargent (2005).
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
#'@param Yraw Data matrix (can be a time series object). Each of \eqn{M} columns
#'  is assumed to contain a single time-series of length \eqn{T}.
#'@param p Integer indicating the order of the VAR, i.e. the number of lags of
#'  the dependent variables included as predictors.
#'@param intercept Either \code{intercept=FALSE} and no constant term
#'  (intercept) will be included. Or a numeric vector of length \eqn{M}
#'  indicating the (fixed) prior variances on the constant term. A single number
#'  will be recycled accordingly. Default is \code{intercept=100}, which imposes
#'  hardly any regularization w.r.t. the intercept.
#'@param draws single number indicating the number of draws after the burnin
#'@param burnin single number indicating the number of draws discarded as burnin
#'@param persistence single number greater or equal to 0 indicating the prior
#'  mean of the first own-lag coefficients (default is 0).
#'@param priorPHI List object from \code{\link{specify_priorPHI}} with prior
#'  settings. Used to select and customize the prior on the reduced-form VAR
#'  coefficients.
#'@param priorL List object from \code{\link{specify_priorL}} with prior
#'  settings. Used to select and customize the prior on the constant covariance
#'  (precision) parameters.
#'@param SV logical indicating whether time-varying (\code{TRUE}) or constant
#'  (\code{FALSE}) variance should be estimated.
#'@param sv_spec list with stochastic volatility specification if
#'  \code{SV=TRUE}. Must contain
#'@param priorHomoscedastic Only used if \code{SV=FALSE}. In that case
#'  \code{priorHomoscedastic} must be a \eqn{M \times 2} matrix, where the
#'  first column entries indicate the shape and the second column entries
#'  indicate the scale parameters of the inverse Gamma prior distributions of
#'  the orthogonalized variances. All entries must be greater than 0. A vector
#'  of length 2 will be recycled accordingly.
#'@param progressbar logical value indicating whether a progressbar should be
#'  displayed during sampling (default is \code{TRUE}).
#'@param PHI_in numeric matrix used as initial draw for the reduced-form VAR
#'  coefficients. Must be of dimension c(K,M), where K is the number of
#'  covariables per equation and M the number of endogenous variables. If not
#'  specified the posterior mean of a flat conjugate prior will be used.
#'@param L_in numeric upper unitriangular matrix (with ones on the diagonal) of
#'  dimension c(M,M). Used as initial draw for covariance paramters. If not
#'  specified the normalized cholesky factor of the posterior mean of the
#'  variance-covariance matrix of a flat conjugate prior will be used.
#'
#'@export
bvar <- function(Yraw,
                      p,
                      intercept = 100,
                      draws=1000,
                      burnin=1000,
                      persistence = 0,
                      priorPHI,
                      priorL,
                      SV=TRUE,
                      sv_spec = list(priormu = c(0,100),
                                     priorphi = c(20, 1.5),
                                     priorsigma2 = c(0.5,0.5),
                                     sv_offset = 0),
                      priorHomoscedastic = NA,
                      progressbar=TRUE,
                      PHI_in=NULL,
                      L_in=NULL
                      ){

# Data preliminaries ------------------------------------------------------

  # M: number of variables,
  # T: number of observations used for estimation,
  # K: number of covariates per equation (without intercepts!!!),
  # n: number of VAR coefficients (without intercepts!!!),
  # n_L: number of free off-diagonal elements in L
  # Y: Yraw without first p observations,
  # X: lagged values of Yraw

  M <- ncol(Yraw)
  K <- p*M
  Traw <- nrow(Yraw)
  Y_tmp <- as.matrix(Yraw)
  if (any(is.na(Y_tmp))){
    stop("\nNAs in Yraw.\n")
  }
  if (ncol(Y_tmp) < 2) {
    stop("The matrix 'Yraw' should contain at least two variables. \n")
  }
  if (is.null(colnames(Y_tmp))) {
    colnames(Y_tmp) <- paste("y", 1:ncol(Y_tmp), sep = "")
    warning(paste("No column names supplied in Yraw, using:",
                  paste(colnames(Y_tmp), collapse = ", "), ", instead.\n"))
  }
  colnames(Y_tmp) <- make.names(colnames(Y_tmp))
  # embed: get lagged values
  X <- embed(Y_tmp, dimension = p + 1)[, -(1:M)]
  colnames(X) <- paste0(colnames(Y_tmp), ".l", sort(rep(1:p,M)))
  if(is.numeric(intercept)){
    X <- cbind(X,1)
    colnames(X)[ncol(X)] <- c("intercept")
  }
  Y <- Y_tmp[-c(1:p), ]

  T <- Traw - p
#  if(T!=nrow(Y) | T!=nrow(X)){
#    stop("Something went wrong: T != nrow(Y). \n")
#  }
  n <- K*M
  n_L <- (M^2 - M)/2

# Indicator matrix --------------------------------------------------------

  if(is.numeric(intercept)){
    i_intercept <- rep(0,M)
  }else i_intercept <- NULL

    if(priorPHI$prior == "DL" | priorPHI$prior == "R2D2" |
       priorPHI$prior == "DL_h" | priorPHI$prior == "SSVS"){

      if(all(is.numeric(priorPHI$global_grouping))){
        i_mat <- priorPHI$global_grouping
        if(!(identical(dim(i_mat), as.integer(c(K,M))) & all(is.numeric(i_mat))) )
          stop("Something went wrong specifying 'global_grouping'.")
      }else if(priorPHI$global_grouping=="global"){
        i_mat <- matrix(1, K, M)
      }else if(priorPHI$global_grouping=="fol"){
        i_mat <- matrix(1, K, M)
        diag(i_mat[1:M,1:M]) <- 2
      }else if(priorPHI$global_grouping=="olcl-lagwise"){
        i_mat <- matrix(1, K, M)
        diag(i_mat[1:M,1:M]) <- 2
        if(p>1){
          for (j in 1:(p-1)) {
            i_mat[(j*M+1):((j+1)*M),] <- i_mat[((j-1)*M+1):(j*M),] + 2
          }
        }
      }else if(priorPHI$global_grouping == "equation-wise"){
        i_mat <- matrix(
          rep(1:M, K),
          K,M,
          byrow = TRUE
        )
      }else if(priorPHI$global_grouping == "covariate-wise"){
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
      for (j in seq_len(p)) {
        i_mat <- rbind(i_mat, i_mat_1 * j)
      }
    }

    i_mat <- rbind(i_mat, i_intercept)
    i_vec <- as.vector(i_mat)

# Hyperparameter settigns -------------------------------------------------

  if(!intercept){
    intercept <- 0
    priorIntercept <- vector("numeric")
  }else{
    ##some checks???
    priorIntercept <- rep(intercept, M)
    intercept <- 1
  }

  if(persistence == 0) {

    PHI0 <- matrix(0, K+intercept, M)

  }else {

    PHI0 <- matrix(0, K+intercept, M)
    PHI0[1:M, 1:M] <- diag(M)*persistence

  }

  if(!(priorPHI$prior %in% c("DL", "DL_h", "HMP", "SSVS", "normal", "R2D2", "SL"))){
    stop("Argument 'priorPHI$prior' must be one of
           'DL', 'SSVS', 'HMP' or 'normal'. \n")
  }

# Initialize --------------------------------------------------------------

  if(priorPHI$prior == "SSVS" | is.null(PHI_in) | is.null(L_in)){
    # Posterior mean of a flat conjugate Normal inverse Wishart prior
    # exists even when OLS estimate does not exist (in situations where T < K)
    # N(0, 10^3) on PHI, and invWish(I, M+2) on Sigma
    XX <- crossprod(X)
    V_post_flat <- chol2inv(chol(diag(1/rep(10^3, (K+intercept))) + XX))
    PHI_flat <- V_post_flat %*% (diag(1/rep(10^3, (K+intercept)))%*%PHI0 + t(X)%*%Y)
    S_post <- diag(M) + crossprod(Y - X%*%PHI_flat) + t(PHI_flat - PHI0) %*%
      diag(1/rep(10^3, (K+intercept))) %*% (PHI_flat - PHI0)
    Sigma_flat <- (S_post)/(M +2 + T - M - 1)
    U <- chol(Sigma_flat)
    D <- diag(U)^2
    L_inv <- U/sqrt(D)
    L_flat <- backsolve(L_inv, diag(M))
  }

    #some proper checks missing!!!
  if(is.null(PHI_in)){
    PHI_in <- PHI_flat
  }
  if(is.null(L_in)){
    L_in <- L_flat
  }


  sv_para_init <- matrix(data= c(rep(-10,M), rep(0.9,M), rep(0.2,M), rep(-10,M)),
                                nrow = 4, ncol = M, byrow = TRUE)
  rownames(sv_para_init) <- c("mu", "phi", "sigma", "h0")

  h_init <- matrix(rep(-10, T*M), T,M)

  # creating placeholders (for cpp, maybe move to cpp code)
  priorPHI_in <- list()
  priorPHI_in$prior <- priorPHI$prior

  #GL priors
  priorPHI_in$n_groups <- 1
  priorPHI_in$groups <- 1

  #DL
  priorPHI_in$DL_a <- double(1)
  priorPHI_in$a_vec <- double(1)
  priorPHI_in$prep1 <- double(1)
  priorPHI_in$prep2 <- matrix(0)
  priorPHI_in$DL_hyper <- logical(1)
  #R2D2
  priorPHI_in$R2D2_hyper <- logical(1)
  priorPHI_in$R2D2_api <- double(1)
  priorPHI_in$R2D2_b <- double(1)
  priorPHI_in$api_vec <- double(1)
  priorPHI_in$b_vec <- double(1)
  #SSVS
  priorPHI_in$SSVS_tau0 <- double(1)
  priorPHI_in$SSVS_tau1 <- double(1)
  priorPHI_in$SSVS_s_a <- double(1)
  priorPHI_in$SSVS_s_b <- double(1)
  priorPHI_in$SSVS_p <- double(1)
  priorPHI_in$SSVS_hyper <- logical(1)
  #HM
  priorPHI_in$lambda_1 <- double(1)
  priorPHI_in$lambda_2 <- double(1)
  priorPHI_in$V_i_prep <- double(1)

  # prior specification for PHI

  if(priorPHI$prior == "R2D2" | priorPHI$prior == "DL" |
     priorPHI$prior == "DL_h" | priorPHI$prior == "SSVS"){

    groups <- unique(i_vec[i_vec!=0])
    n_groups <- length(groups)
    priorPHI_in$n_groups <- n_groups
    priorPHI_in$groups <- groups

  }

  if(priorPHI$prior == "R2D2"){

    if(all(is.numeric(priorPHI$R2D2_b))){
      R2D2_hyper <- FALSE
    }else if("hyperprior" %in% priorPHI$R2D2_b & length(priorPHI$R2D2_b) == 1){

      R2D2_hyper <- TRUE

    }else{
      stop("Something went wrong specifying R2D2_b!")
    }
    priorPHI_in$R2D2_hyper <- R2D2_hyper

    if(!R2D2_hyper){

      priorPHI_in$R2D2_b <- rep_len(priorPHI$R2D2_b, n_groups)
      priorPHI_in$R2D2_api <- rep(as.numeric(NA), n_groups)
      for(j in seq.int(n_groups)){
        priorPHI_in$R2D2_api[j] <- 1/(n^(priorPHI$R2D2_b[j]/2) *
                                     (T^(priorPHI$R2D2_b[j]/2)) *log(T))
      }

    }else if(R2D2_hyper){

      grid <- 100
      priorPHI_in$b_vec <- seq(.01,1,length.out=grid)
      priorPHI_in$api_vec <- rep(as.numeric(NA), grid)
      for(j in seq.int(grid)){
        priorPHI_in$api_vec[j] <- 1/(n^(priorPHI_in$b_vec[j]/2) *
                                    (T^(priorPHI_in$b_vec[j]/2)) *log(T))
      }
      priorPHI_in$R2D2_b <- rep(0.5, n_groups)
      priorPHI_in$R2D2_api <- rep(priorPHI_in$api_vec[which(priorPHI_in$b_vec==0.5)],
                               n_groups)

    }


  }else if(priorPHI$prior == "DL"){

    if(all(is.numeric(priorPHI$DL_a))){
      DL_hyper <- FALSE
    }else if("hyperprior" %in% priorPHI$DL_a & length(priorPHI$DL_a) == 1){

      DL_hyper <- TRUE
    }else if(any(c("1/K", "1/n") %in% priorPHI$DL_a) &
             length(priorPHI$DL_a) == 1){
      DL_hyper <- FALSE
      if(priorPHI$DL_a == "1/K") {
        DL_a <- rep(1/K, n_groups)
      }else if(priorPHI$DL_a == "1/n") {
        DL_a <- rep(1/n, n_groups)
      }
    }else{
      stop("Something went wrong specifying DL_a!")
    }
    priorPHI_in$DL_hyper <- DL_hyper
    if(DL_hyper){#priorPHI$
      grid <- 1000
      priorPHI_in$a_vec <- seq(1/(n),1/2,length.out = grid)
      # prep1 &prep2: some preparations for the evaluation of the
      # Dirichlet-density
      priorPHI_in$prep1 <- priorPHI_in$a_vec-1
      prep2 <- matrix(NA, n_groups, grid)
      ii <- 1
      for(j in groups){
        n_tmp <- length(which(i_vec == j))
        # normalizing constant of symmetric Dirichlet density in logs
        prep2[ii,] <- lgamma(n_tmp*priorPHI_in$a_vec) - n_tmp*lgamma(priorPHI_in$a_vec)
        ii <- ii + 1
      }
      priorPHI_in$prep2 <- prep2
      priorPHI_in$DL_a <- rep_len(1/2, n_groups) #initial value
    }else if(!DL_hyper){
      if(all(is.numeric(priorPHI$DL_a))){
        priorPHI_in$DL_a <- rep_len(priorPHI$DL_a, n_groups)
      }
    }

  }else if(priorPHI$prior == "SSVS"){
    priorPHI_in$SSVS_hyper <- priorPHI$SSVS_hyper
    priorPHI_in$SSVS_p <- rep_len(priorPHI$SSVS_p,n)
    if(priorPHI$semiautomatic) {
      # shape parameters of beta hyperprior on prior-inclusion-probability
      priorPHI_in$SSVS_s_a <- priorPHI$SSVS_s_a
      priorPHI_in$SSVS_s_b <- priorPHI$SSVS_s_b

      # standard errors of flat phi posterior estimates
      sigma_phi <- sqrt(diag(Sigma_flat %x% V_post_flat))
      # select all, except those of intercept
      sigma_phi <- sigma_phi[which(i_vec!=0)]
      if(length(sigma_phi)!=p*M^2){
        stop(length(sigma_phi))
      }

      # scale tau with variances
      priorPHI_in$SSVS_tau0 <- priorPHI$SSVS_c0*sigma_phi
      priorPHI_in$SSVS_tau1 <- priorPHI$SSVS_c1*sigma_phi
    }else{
      priorPHI_in$SSVS_tau0 <- rep_len(priorPHI$SSVS_c0, n)
      priorPHI_in$SSVS_tau1 <- rep_len(priorPHI$SSVS_c1, n)
    }

  }else if(priorPHI$prior == "normal"){

    priorPHI_in$V_i <- rep(priorPHI$V_i, length = n)

  }else if(priorPHI$prior == "HMP"){

    priorPHI_in$lambda_1 <- priorPHI$lambda_1
    priorPHI_in$lamda_2 <- priorPHI$lambda_2
    sigma_sq <- MP_sigma_sq(Yraw, 6)

    # prepare prior variances down to lambdas
    priorPHI_in$V_i_prep <- MP_V_prior_prep(sigma_sq, (K+intercept), M, intercept>0)
  }

  # creating placeholders (for cpp, maybe move to cpp code)
  priorL_in <- list()
  priorL_in$prior <- priorL$prior

  #DL
  priorL_in$DL_b <- double(1)
  priorL_in$b_vec <- double(1)
  priorL_in$prep1 <- double(1)
  priorL_in$prep2 <- double(1)
  priorL_in$DL_hyper <- logical(1)
  #R2D2
  priorL_in$R2D2_hyper <- logical(1)
  priorL_in$R2D2_api <- double(1)
  priorL_in$R2D2_b <- double(1)
  priorL_in$api_vec <- double(1)
  priorL_in$b_vec <- double(1)
  #SSVS
  priorL_in$SSVS_tau0 <- double(1)
  priorL_in$SSVS_tau1 <- double(1)
  priorL_in$SSVS_s_a <- double(1)
  priorL_in$SSVS_s_b <- double(1)
  priorL_in$SSVS_hyper <- logical(1)
  priorL_in$SSVS_p <- double(1L)
  #HM
  priorL_in$lambda_3 <- double(1)

  # prior specification for L

  if(priorL$prior == "DL"){
    if(is.numeric(priorL$DL_b)){
      priorL_in$DL_hyper <- FALSE
      priorL_in$DL_b <- priorL$DL_b
    }else if(priorL$DL_b == "1/n") {
      priorL_in$DL_b <- 1/n_L
      priorL_in$DL_hyper <- FALSE
    }else if(priorL$DL_b == "hyperprior"){

      priorL_in$DL_hyper <- TRUE

      grid_L <- 1000
      priorL_in$b_vec <- seq(1/(n_L),1/2,length.out = grid_L)
      priorL_in$prep1 <- priorL_in$b_vec - 1
      priorL_in$prep2 <- lgamma(n_L*priorL_in$b_vec) - n_L*lgamma(priorL_in$b_vec)
      priorL_in$DL_b <- 1/2
    }

  }else if(priorL$prior == "SSVS"){
    priorL_in$SSVS_hyper <- priorL$SSVS_hyper
    priorL_in$SSVS_p <- rep_len(priorL$SSVS_p, n_L)
    priorL_in$SSVS_tau0 <- rep(priorL$SSVS_c0, n_L)
    priorL_in$SSVS_tau1 <- rep(priorL$SSVS_c1, n_L)
    priorL_in$SSVS_s_a <- priorL$SSVS_s_a
    priorL_in$SSVS_s_b <- priorL$SSVS_s_b

  }else if(priorL$prior == "normal"){

    priorL_in$V_i <- rep(priorL$V_i, length = n_L)

  }else if(priorL$prior == "R2D2"){

    if(is.numeric(priorL$R2D2_b)){
      priorL_in$R2D2_hyper <- FALSE
    }else if(priorL$R2D2_b == "hyperprior"){

      priorL_in$R2D2_hyper <- TRUE

    }else{
      stop("Something went wrong specifying R2D2_b (specify_L...)!")
    }

    if(!priorL_in$R2D2_hyper){

      priorL_in$R2D2_b <- priorL$R2D2_b
      priorL_in$R2D2_api <- 1/(n_L^(priorL_in$R2D2_b/2) *
                                     (T^(priorL_in$R2D2_b/2)) *log(T))


    }else if(priorL_in$R2D2_hyper){

      grid <- 100
      priorL_in$b_vec <- seq(.01,1,length.out=grid)
      priorL_in$api_vec <- rep(as.numeric(NA), grid)
      for(j in seq.int(grid)){
        priorL_in$api_vec[j] <- 1/(n_L^(priorL_in$b_vec[j]/2) *
                                    (T^(priorL_in$b_vec[j]/2)) *log(T))
      }
      priorL_in$R2D2_b <- 0.5
      priorL_in$R2D2_api <- priorL_in$api_vec[which(priorL_in$b_vec==0.5)]

    }
  }else if(priorL$prior == "HMP"){
    priorL_in$lambda_3 <- priorL$lambda_3
  }


#  if(is.null(sv_spec)){
#    sv_spec <- list(priormu = c(0,100),
#                    priorphi = c(20, 1.5),
#                    priorsigma2 = c(0.5,0.5)#,
#                    #priorh0 = -1 #h0 from stationary distribution
#    )
#  }
  if(SV == TRUE){
    sv_spec_error <- "sv_spec is a list that must supply the following
                      components: \n
                      priormu: a numeric vector of length 2, where the
                      second element must be posistive \n
                      priorphi: a strictly positive numeric vector of length 2 \n
                      priorsigma: a strictly positive numeric vector of length 2 \n
                      sv_offset: either a single non-negative number, or a vector of length M with non-negative entries. \n"
    sv_specs <- c("priormu", "priorphi", "priorsigma2", "sv_offset")
    if(any(!(sv_specs %in% names(sv_spec)))){
      stop(sv_spec_error)
    }
    if(!(length(sv_spec[["priormu"]]) == 2 & is.numeric(sv_spec[["priormu"]]) &
       sv_spec[["priormu"]][2]>0)){
      stop(sv_spec_error)
    }
    if(!(length(sv_spec[["priorphi"]]) == 2 & is.numeric(sv_spec[["priorphi"]]) &
       all(sv_spec[["priorphi"]]>0))){
      stop(sv_spec_error)
    }
    if(!(length(sv_spec[["priorsigma2"]]) == 2 & is.numeric(sv_spec[["priorsigma2"]]) &
         all(sv_spec[["priorsigma2"]]>0))){
      stop(sv_spec_error)
    }
    if(length(sv_spec[["sv_offset"]]) == 1 | length(sv_spec[["sv_offset"]]) == M){
      sv_spec$sv_offset <- rep_len(sv_spec$sv_offset, M)
    }else {
      stop(sv_spec_error)
      }


    if(!(length(priorHomoscedastic) == 1 )){
      warning("'priorHomoscedastic' setting will be ignored, because 'SV=TRUE'! \n",
              immediate. = TRUE)
    }else if(length(priorHomoscedastic) == 1){
      if(!is.na(priorHomoscedastic)){
        warning("'priorHomoscedastic' setting will be ignored, because 'SV=TRUE'! \n",
                immediate. = TRUE)
      }
    }
    priorHomoscedastic <- matrix(as.numeric(NA), M, 2) #0.01
  }else if(SV == FALSE){
    ph_error <- paste0("priorHomoscedastic must be either a numeric matrix of dimension c(M,2),
             or a numeric vector of length 2, where all entries are greater than 0. \n")
    if(!identical(dim(priorHomoscedastic), as.integer(c(M,2)))){
      if(length(priorHomoscedastic) == 2 & is.numeric(priorHomoscedastic) &
         all(priorHomoscedastic>0)){
        priorHomoscedastic <- matrix(priorHomoscedastic, M, 2, byrow = TRUE)
      }else {
        stop(ph_error)
      }
    }else if(identical(dim(priorHomoscedastic), as.integer(c(M,2))) &
             is.numeric(priorHomoscedastic)){
      if(!all(priorHomoscedastic>0)){
        stop(ph_error)
      }
    }else{
      stop(ph_error)
    }
  }

  res <- bvar_cpp(Y,
                  X,
                  M,
                  T,
                  K,
                  draws,
                  burnin,
                  intercept,
                  priorIntercept,
                  PHI_in,
                  PHI0,
                  priorPHI_in,
                  priorL_in,
                  L_in,
                  SV,
                  priorHomoscedastic,
                  sv_spec,
                  h_init,
                  sv_para_init,
                  i_mat,
                  i_vec,
                  progressbar
                  )

  #Rcpp timer is in nanoseconds
  #conversion to secs per iteration
  res$bench <- diff(res$bench)/(10^(9)*(draws+burnin))
  attributes(res$bench) <- list("names" = "secs/itr")
  dimnames(res$PHI)[2] <- list(colnames(X))
  dimnames(res$PHI)[3] <- list(colnames(Y))
  dimnames(res$L)[2] <- dimnames(res$L)[3] <- list(colnames(Y))
  phinames <- as.vector((vapply(seq_len(M), function(i) paste0(colnames(Y)[i], "~", colnames(X[,1:(ncol(X)-intercept)])), character(K))))
  if(priorPHI$prior %in% c("DL","DL_h")){
    colnames(res$phi_hyperparameter) <- c(paste0("zeta",1:priorPHI_in$n_groups), paste0("psi: ", phinames), paste0("theta: ", phinames), paste0("a",1:priorPHI_in$n_groups))#

  }else if(priorPHI$prior == "R2D2"){
    colnames(res$phi_hyperparameter) <- c(paste0("zeta",1:priorPHI_in$n_groups),
                                          paste0("psi: ", phinames),
                                          paste0("theta: ", phinames),
                                          paste0("xi",1:priorPHI_in$n_groups),
                                          paste0("b",1:priorPHI_in$n_groups),
                                          paste0("api",1:priorPHI_in$n_groups))#

  }else if(priorPHI$prior == "SSVS"){
  colnames(res$phi_hyperparameter) <- c(paste0("gamma: ", phinames), paste0("p_i: ", phinames))#

  }else if(priorPHI$prior == "HMP"){
    colnames(res$phi_hyperparameter) <- c("lambda_1", "lambda_2")
  }
  lnames <- NULL
  for(j in 2:M){
    lnames <- c(lnames, paste0(colnames(Y)[j],"~", colnames(Y)[1:(j-1)]))
  }
  if(priorL$prior %in% c("DL","DL_h")){
    colnames(res$l_hyperparameter) <- c("zeta", paste0("psi: ", lnames), paste0("theta: ", lnames), "b")
  }else if(priorL$prior == "R2D2"){
    colnames(res$l_hyperparameter) <- c("zeta", paste0("psi: ", lnames), paste0("theta: ", lnames), "xi", "b", "api")#

  }else if(priorL$prior == "SSVS"){
    colnames(res$l_hyperparameter) <- c(paste0("gamma: ", lnames), paste0("p_i: ", lnames))
  }else if(priorL$prior == "HMP"){
    colnames(res$l_hyperparameter) <- c("lambda_3")
  }
  res$phi_hyperparameter <- as.data.frame(res$phi_hyperparameter)
  res$l_hyperparameter <- as.data.frame(res$l_hyperparameter)
  res$Y <- Y
  res$X <- X
  res$p <- p
  res$intercept <- ifelse(intercept>0, TRUE, FALSE)
  res$SV <- SV
  res$Yraw <- Y_tmp
  res$Traw <- Traw
  res$datamat <- data.frame(cbind(Y, X))
  class(res$PHI) <- "PHI"
  class(res$L) <- "L"
  class(res) <- "bvar"
  res
}
