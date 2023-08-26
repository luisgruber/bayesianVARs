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
                 draws=1000L,
                 burnin=1000L,
                 thin = 1L,
                 persistence = 0,
                 priorPHI,
                 priorSigma, #priorL,
###                 SV=TRUE,
###                 sv_spec = list(priormu = c(0,100),
###                                priorphi = c(20, 1.5),
###                                priorsigma2 = c(0.5,0.5),#
###                                sv_offset = 0),
###                 priorHomoscedastic = NA,
                 tvp_keep = "last",
                 progressbar=TRUE,
                 PHI_in=NULL,
                 L_in=NULL,
                 PHI_tol = 1e-100,
                 L_tol = 1e-100,
                 expert = list()
                      ){

# Data preliminaries ------------------------------------------------------

  # M: number of variables,
  # T: number of observations used for estimation,
  # K: number of covariates per equation (without intercepts!!!),
  # n: number of VAR coefficients (without intercepts!!!),
  # n_L: number of free off-diagonal elements in L
  # Y: Yraw without first p observations,
  # X: lagged values of Yraw

  thin <- as.integer(thin)
  if(thin<1) stop("Argument 'thin' must be an integer greater than 0!")
  draws <- as.integer(draws)
  burnin <- as.integer(burnin)
  if(!tvp_keep%in%c("last", "all")){
    stop("Argument 'tvp_keep' must be one of 'last' or 'all'.")
  }

  p <- as.integer(p)
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
       priorPHI$prior == "HS" | priorPHI$prior == "SSVS" |
       priorPHI$prior == "NG" | priorPHI$prior == "GT"){

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
    mode(i_mat) <- "integer"
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

  if(!(priorPHI$prior %in% c("DL", "HS","NG", "HMP", "SSVS", "normal", "R2D2", "GT"))){
    stop("Argument 'priorPHI$prior' must be one of
           'DL', 'R2D2', 'HS', 'NG', 'SSVS', 'HMP' or 'normal'. \n")
  }

# Initialize PHI --------

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


# Hyperparameter settings w.r.t PHI ---------------------------------------

  # creating placeholders (for cpp, maybe move to cpp code)
  priorPHI_in <- list()
  priorPHI_in$prior <- priorPHI$prior

  #GL priors
  priorPHI_in$n_groups <- 1
  priorPHI_in$groups <- 1
  priorPHI_in$GL_tol <- double(1L)
  priorPHI_in$a <- double(1L)
  priorPHI_in$b <- double(1L)
  priorPHI_in$c <- double(1L)
  priorPHI_in$GT_vs <- double(1L)
  priorPHI_in$GT_priorkernel <- character(1L)
  priorPHI_in$a_vec <- double(1)
  priorPHI_in$a_weight <- double(1)
  priorPHI_in$norm_consts <- double(1)
  priorPHI_in$c_vec <- double(1)
  priorPHI_in$c_rel_a <- logical(1L)

  #DL
##?  priorPHI_in$GL_tol <- double(1)
##?  priorPHI_in$a <- double(1)
##?  priorPHI_in$DL_b <- double(1)
##?  priorPHI_in$DL_c <- double(1)
  priorPHI_in$prep2 <- matrix(0)
  priorPHI_in$prep1 <- double(1)
  priorPHI_in$DL_hyper <- logical(1)
  priorPHI_in$DL_plus <- logical(1)
  #GT (Normal Gamma and R2D2 belong to GT)
  priorPHI_in$GT_hyper <- logical(1)
  #R2D2
##?  priorPHI_in$R2D2_method <- integer(1)
##?  priorPHI_in$R2D2_kernel <- character(1L)
##?  priorPHI_in$R2D2_hyper <- logical(1)
##?  priorPHI_in$R2D2_api <- double(1)
##?  priorPHI_in$R2D2_b <- double(1)
##?  priorPHI_in$api_vec <- double(1)
##?  priorPHI_in$b_vec <- double(1)
  #NG
##?  priorPHI_in$NG_a <- double(1)
##?  priorPHI_in$NG_varrho0 <- double(1)
##?  priorPHI_in$NG_varrho1 <- double(1)
##?  priorPHI_in$NG_hyper <- logical(1)
##?  priorPHI_in$NG_a_vec <- double(1)
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
     priorPHI$prior == "DL_h" | priorPHI$prior == "SSVS" |
     priorPHI$prior == "HS" | priorPHI$prior == "NG" | priorPHI$prior == "GT"){

    groups <- unique(i_vec[i_vec!=0])
    n_groups <- length(groups)
    priorPHI_in$n_groups <- n_groups
    priorPHI_in$groups <- groups

  }

  if(priorPHI$prior == "GT"){

    priorPHI_in$GT_priorkernel <- priorPHI$GT_priorkernel
    priorPHI_in$GL_tol <- priorPHI$GL_tol
    priorPHI_in$GT_vs <- priorPHI$GT_vs
    priorPHI_in$b <- rep_len(priorPHI$b, n_groups)

    if(is.matrix(priorPHI$a)){
      if(ncol(priorPHI$a)==2){
        priorPHI_in$GT_hyper <- TRUE
        priorPHI_in$a_vec <- priorPHI$a[,1]
        priorPHI_in$a_weight <- priorPHI$a[,2]
        priorPHI_in$norm_consts <- lgamma(priorPHI_in$a_vec)
        priorPHI_in$a <- sample(priorPHI_in$a_vec, n_groups, replace = TRUE, prob = priorPHI_in$a_weight) # initialize a
      }else if(ncol(priorPHI$a)>2){
        stop("The easiest way to specify 'R2D2_a', 'NG_a' or 'GT_a' is a single postive number!")
      }else{
        priorPHI$a <- as.vector(priorPHI$a)
      }
    }
    if(is.null(dim(priorPHI$a))){
      priorPHI_in$GT_hyper <- FALSE
      priorPHI_in$a <- rep_len(priorPHI$a, n_groups)
    }

    if(is.character(priorPHI$c)){
      priorPHI_in$c_rel_a <- TRUE # then c is always proportion of a (e.g. for R2D2 c=0.5a)
      mya <- priorPHI_in$a
      myc <- gsub("a","mya", priorPHI$c)
      priorPHI_in$c <- eval(str2lang(myc))
      if(base::isTRUE(priorPHI_in$GT_hyper)){
        myc2 <- gsub("a","priorPHI_in$a_vec", priorPHI$c)
        priorPHI_in$c_vec <- eval(str2lang(myc2))
      }
    }else if(is.numeric(priorPHI$c)){
      priorPHI_in$c_rel_a <- FALSE
      priorPHI_in$c <- rep_len(priorPHI$c, n_groups)
    }


  }else if(priorPHI$prior == "DL"){

    priorPHI_in$GL_tol <- priorPHI$GL_tol
    # if(!exists("DL_method", priorPHI)){
    #   priorPHI_in$DL_method <- 2L
    # }else{
    #   priorPHI_in$DL_method <- priorPHI$DL_method
    # }

    if(all(is.numeric(priorPHI$a))&is.null(dim(priorPHI$a))){
      DL_hyper <- FALSE
    }else if(("hyperprior" %in% priorPHI$a & length(priorPHI$a) == 1) ||
             (is.matrix(priorPHI$a))){

      DL_hyper <- TRUE
    }else if(any(c("1/K", "1/n") %in% priorPHI$a) &
             length(priorPHI$a) == 1){
      DL_hyper <- FALSE
      if(priorPHI$a == "1/K") {
        priorPHI_in$a <- rep(1/K, n_groups)
      }else if(priorPHI$a == "1/n") {
        priorPHI_in$a <- rep(1/n, n_groups)
      }
    }else{
      stop("Something went wrong specifying DL_a!")
    }
    priorPHI_in$DL_hyper <- DL_hyper
    if(DL_hyper){
      if(is.matrix(priorPHI$a)){
        priorPHI_in$a_vec <- priorPHI$a[,1]
        priorPHI_in$a_weight <- priorPHI$a[,2]
        priorPHI_in$norm_consts <- lgamma(priorPHI_in$a_vec)
        priorPHI_in$a <- rep_len(priorPHI_in$a_vec[1], n_groups) #initial value
        #priorPHI_in$DL_method <- 2L
      }else{
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
        priorPHI_in$a <- rep_len(1/2, n_groups) #initial value
        #if(priorPHI_in$DL_method==2L){
          priorPHI_in$a_weight <- rep(1,grid)
          priorPHI_in$norm_consts <- lgamma(priorPHI_in$a_vec)
        #}
      }

    }else if(!DL_hyper){
      if(all(is.numeric(priorPHI$a))){
        priorPHI_in$a <- rep_len(priorPHI$a, n_groups)
      }
    }

    if(!exists("DL_plus", priorPHI) || base::isFALSE(priorPHI$DL_plus)){
      priorPHI_in$DL_plus <- FALSE
    }else if(priorPHI$DL_plus){
      priorPHI_in$DL_plus <- TRUE
      if(!exists("DL_b", priorPHI)){
        priorPHI_in$b <- rep_len(0.5, n_groups)
      }else{
        priorPHI_in$b <- rep_len(priorPHI$DL_b, n_groups)
      }
      if(!exists("DL_c", priorPHI)){
        priorPHI_in$c <- 0.5*priorPHI_in$a
        if(DL_hyper){
          priorPHI_in$c_vec <- 0.5*priorPHI_in$a_vec
          priorPHI_in$c_rel_a <- TRUE
        }
      }else{
        if(!is.numeric(priorPHI$DL_c)){
          stop("If you specify DL_c for model DL_plus, then it must be numeric!")
        }
        priorPHI_in$c <- rep_len(priorPHI$DL_c, n_groups)
        priorPHI_in$c_rel_a <- FALSE

      }
      # priorPHI_in$DL_method <- 2L
    }else{
      stop("Never heard of DL_plus?")
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
  # else if(priorPHI$prior == "R2D2"){
  #
  #   if(!exists("R2D2_method", priorPHI)){
  #     priorPHI_in$R2D2_method <- 1L
  #   }else{
  #     priorPHI_in$R2D2_method <- priorPHI$R2D2_method
  #   }
  #
  #   if(!exists("R2D2_kernel", priorPHI)){
  #     priorPHI_in$R2D2_kernel <- "laplace"
  #   }else{
  #     priorPHI_in$R2D2_kernel <- priorPHI$R2D2_kernel
  #   }
  #
  #   if(all(is.numeric(priorPHI$R2D2_b))){
  #     R2D2_hyper <- FALSE
  #   }else if("hyperprior" %in% priorPHI$R2D2_b & length(priorPHI$R2D2_b) == 1){
  #
  #     R2D2_hyper <- TRUE
  #
  #   }else{
  #     stop("Something went wrong specifying R2D2_b!")
  #   }
  #   priorPHI_in$R2D2_hyper <- R2D2_hyper
  #
  #   if(!R2D2_hyper){
  #
  #     priorPHI_in$R2D2_b <- rep_len(priorPHI$R2D2_b, n_groups)
  #     if(priorPHI$R2D2_api=="default"){
  #       priorPHI_in$R2D2_api <- rep(as.numeric(NA), n_groups)
  #       for(j in seq.int(n_groups)){
  #         priorPHI_in$R2D2_api[j] <- 1/(n^(priorPHI_in$R2D2_b[j]/2) *
  #                                         (T^(priorPHI_in$R2D2_b[j]/2)) *log(T))
  #       }
  #     }else{
  #       priorPHI_in$R2D2_api <- priorPHI$R2D2_api
  #     }
  #
  #
  #   }else if(R2D2_hyper){
  #
  #     grid <- 100
  #     priorPHI_in$b_vec <- seq(.01,1,length.out=grid)
  #     if(priorPHI$R2D2_api=="default"){
  #       priorPHI_in$api_vec <- rep(as.numeric(NA), grid)
  #       for(j in seq.int(grid)){
  #         priorPHI_in$api_vec[j] <- 1/(n^(priorPHI_in$b_vec[j]/2) *
  #                                        (T^(priorPHI_in$b_vec[j]/2)) *log(T))
  #       }
  #     }else{
  #       priorPHI_in$api_vec <- rep_len(priorPHI$R2D2_api, grid)
  #     }
  #
  #     #initialize
  #     priorPHI_in$R2D2_b <- rep(0.5, n_groups)
  #     priorPHI_in$R2D2_api <- rep(priorPHI_in$api_vec[which(priorPHI_in$b_vec==0.5)],
  #                                 n_groups)
  #   }
  #
  #
  #
  # }else if(priorPHI$prior == "NG"){
  #   #priorPHI_in$NG_a <- priorPHI$NG_a
  #   #priorPHI_in$NG_varrho0 <- priorPHI$NG_varrho0
  #   #priorPHI_in$NG_varrho1 <- priorPHI$NG_varrho1
  #
  #   if(all(is.numeric(priorPHI$NG_a))){
  #     NG_hyper <- FALSE
  #   }else if(priorPHI$NG_a == "hyperprior"){
  #     NG_hyper <- TRUE
  #   }else{
  #     stop("Something went wrong specifying NG_a!")
  #   }
  #   priorPHI_in$NG_hyper <- NG_hyper
  #   if(NG_hyper){
  #     priorPHI_in$NG_a_vec <- seq(0.01,1,length.out=100)
  #     priorPHI_in$NG_a <- rep_len(0.1, n_groups)
  #   }else{
  #     priorPHI_in$NG_a <- rep_len(priorPHI$NG_a, n_groups)
  #   }
  #
  # }

# Sigma -------------------------------------------------------------------


  if(!isa(priorSigma, "specify_Sigma")){
    stop("\nArgument 'priorSigma' must be of class 'specify_Sigma'. Please use helper function 'specify_Sigma()'!\n")
  }
  if(priorSigma$M != M){
    stop(paste0("\n'priorSigma' is specified for ",priorSigma$M, " timeseries. 'Yraw', however, consits of ",M," timeseries!\n"))
  }
    # # creating placeholders (for cpp, maybe move to cpp code)
    # priorL_in <- list()
    # priorL_in$prior <- priorL$prior
    #
    # ## GL priors
    # priorL_in$GL_tol <- double(1L)
    # priorL_in$a <- double(1L)
    # priorL_in$b <- double(1L)
    # priorL_in$c <- double(1L)
    # priorL_in$GT_vs <- double(1L)
    # priorL_in$GT_priorkernel <- character(1L)
    # priorL_in$a_vec <- double(1L)
    # priorL_in$a_weight <- double(1L)
    # priorL_in$norm_consts <- double(1L)
    # priorL_in$c_vec <- double(1)
    # priorL_in$c_rel_a <- logical(1L)
    # priorL_in$GT_hyper <- logical(1)
    #
    # #DL
    # priorL_in$DL_hyper <- logical(1L)
    # priorL_in$DL_plus <- logical(1L)
    # #R2D2
    # priorL_in$R2D2_hyper <- logical(1L)
    # priorL_in$R2D2_api <- double(1L)
    # priorL_in$R2D2_b <- double(1L)
    # priorL_in$api_vec <- double(1L)
    # priorL_in$b_vec <- double(1L)
    # #SSVS
    # priorL_in$SSVS_tau0 <- double(1L)
    # priorL_in$SSVS_tau1 <- double(1L)
    # priorL_in$SSVS_s_a <- double(1L)
    # priorL_in$SSVS_s_b <- double(1L)
    # priorL_in$SSVS_hyper <- logical(1L)
    # priorL_in$SSVS_p <- double(1L)
    # #HM
    # priorL_in$lambda_3 <- double(1L)

  ##priorSigma_in <- priorSigma
  # prior specification for L

  # if(priorSigma$cholesky_priorU == "GT"){
  #
  #   priorSigma_in$cholesky_GT_priorkernel <- priorSigma$cholesky_GT_priorkernel
  #   priorSigma_in$cholesky_GL_tol <- priorSigma$cholesky_GL_tol
  #   priorSigma_in$cholesky_GT_vs <- priorSigma$cholesky_GT_vs
  #  # priorSigma_in$cholesky_a <- priorSigma$cholesky_a
  #   priorSigma_in$cholesky_b <- priorSigma$cholesky_b
  #
  #   if(is.matrix(priorSigma$cholesky_a)){
  #     if(ncol(priorSigma$cholesky_a)==2){
  #       priorSigma_in$cholesky_GT_hyper <- TRUE
  #       priorSigma_in$cholesky_a_vec <- priorSigma$cholesky_a[,1]
  #       priorSigma_in$cholesky_a_weight <- priorSigma$cholesky_a[,2]
  #       priorSigma_in$cholesky_norm_consts <- lgamma(priorSigma_in$cholesky_a_vec)
  #       priorSigma_in$cholesky_a <- sample(priorSigma_in$cholesky_a_vec, 1, replace = TRUE, prob = priorSigma_in$cholesky_a_weight) # initialize a
  #     }else if(ncol(priorSigma$cholesky_a)>2){
  #       stop("The easiest way to specify 'R2D2_a', 'NG_a' or 'GT_a' is a single postive number!")
  #     }else{
  #       priorSigma$cholesky_a <- as.vector(priorSigma$cholesky_a)
  #     }
  #   }
  #   if(is.null(dim(priorSigma$cholesky_a))){
  #     priorSigma_in$cholesky_GT_hyper <- FALSE
  #     priorSigma_in$cholesky_a <- priorSigma$cholesky_a
  #   }
  #
  #   if(is.character(priorSigma$cholesky_c)){
  #     priorSigma_in$cholesky_c_rel_a <- TRUE # then c is always proportion of a (e.g. for R2D2 c=0.5a)
  #     mya <- priorSigma_in$cholesky_a
  #     myc <- gsub("a","mya", priorSigma$cholesky_c)
  #     priorSigma_in$cholesky_c <- eval(str2lang(myc))
  #     if(base::isTRUE(priorSigma_in$cholesky_GT_hyper)){
  #       myc2 <- gsub("a","priorSigma_in$cholesky_a_vec", priorSigma$cholesky_c)
  #       priorSigma_in$cholesky_c_vec <- eval(str2lang(myc2))
  #     }
  #   }else if(is.numeric(priorSigma$cholesky_c)){
  #     priorSigma_in$cholesky_c_rel_a <- FALSE
  #     priorSigma_in$cholesky_c <- priorSigma$cholesky_c
  #   }
  #
  # }else if(priorSigma$cholesky_priorU == "DL"){
  #   priorSigma_in$cholesky_GL_tol <- priorSigma$cholesky_GL_tol
  #   if(is.numeric(priorSigma$cholesky_a) & length(priorSigma$cholesky_a) == 1L){
  #     priorSigma_in$cholesky_DL_hyper <- FALSE
  #     priorSigma_in$cholesky_a <- priorSigma$cholesky_a
  #   }else if(priorSigma$cholesky_a == "1/n") {
  #     priorSigma_in$cholesky_a <- 1/n_L
  #     priorSigma_in$cholesky_DL_hyper <- FALSE
  #   }else if(priorSigma$cholesky_a == "hyperprior"){
  #
  #     priorSigma_in$cholesky_DL_hyper <- TRUE
  #
  #     grid_L <- 1000
  #     priorSigma_in$cholesky_a_vec <- seq(1/(n_L),1/2,length.out = grid_L)
  #     #priorSigma_in$cholesky_prep1 <- priorSigma_in$cholesky_b_vec - 1
  #     #priorSigma_in$cholesky_prep2 <- lgamma(n_L*priorSigma_in$cholesky_b_vec) - n_L*lgamma(priorSigma_in$cholesky_b_vec)
  #     priorSigma_in$cholesky_norm_consts <- 0.5^priorSigma_in$cholesky_a_vec -
  #       lgamma(priorSigma_in$cholesky_a_vec)
  #     priorSigma_in$cholesky_a_weight <- rep(1,grid_L)
  #     priorSigma_in$cholesky_a <- 1/2 # initial value
  #   }else if(is.matrix(priorSigma$cholesky_a)){
  #
  #     if(ncol(priorSigma$cholesky_a)!=2){
  #       stop("If you specify 'DL_a' as a matrix, the first column represents
  #            the support points and the second column the weights of a discrete
  #            hyperprior on 'DL_a' !")
  #     }
  #
  #     priorSigma_in$cholesky_DL_hyper <- TRUE
  #     priorSigma_in$cholesky_a_vec <- priorSigma$cholesky_a[,1]
  #     priorSigma_in$cholesky_a_weight <- priorSigma$cholesky_a[,2]
  #     # precompute log normalizing constants of hyperprior
  #     priorSigma_in$cholesky_norm_consts <- 0.5^priorSigma_in$cholesky_a_vec -
  #       lgamma(priorSigma_in$cholesky_a_vec)
  #     priorSigma_in$cholesky_a <- priorSigma_in$cholesky_a_vec[1] #initial value
  #
  #   }
  #
  #   if(!exists("DL_plus", priorSigma) || base::isFALSE(priorSigma$cholesky_DL_plus)){
  #     priorSigma_in$cholesky_DL_plus <- FALSE
  #   }else if(priorSigma$cholesky_DL_plus){
  #     priorSigma_in$cholesky_DL_plus <- TRUE
  #     if(!exists("DL_b", priorSigma)){
  #       priorSigma_in$cholesky_b <- 0.5
  #     }else{
  #       priorSigma_in$cholesky_b <- priorSigma$cholesky_DL_b
  #     }
  #     if(!exists("DL_c", priorSigma)){
  #       priorSigma_in$cholesky_c <- 0.5*priorSigma_in$cholesky_a
  #     }else{
  #       priorSigma_in$cholesky_c <- priorSigma$cholesky_DL_c
  #     }
  #   }else{
  #     stop("Never heard of DL_plus?")
  #   }
  #
  # }else if(priorSigma$cholesky_priorU == "SSVS"){
  #   priorSigma_in$cholesky_SSVS_hyper <- priorSigma$cholesky_SSVS_hyper
  #   priorSigma_in$cholesky_SSVS_p <- rep_len(priorSigma$cholesky_SSVS_p, n_L)
  #   priorSigma_in$cholesky_SSVS_tau0 <- rep(priorSigma$cholesky_SSVS_c0, n_L)
  #   priorSigma_in$cholesky_SSVS_tau1 <- rep(priorSigma$cholesky_SSVS_c1, n_L)
  #   priorSigma_in$cholesky_SSVS_s_a <- priorSigma$cholesky_SSVS_s_a
  #   priorSigma_in$cholesky_SSVS_s_b <- priorSigma$cholesky_SSVS_s_b
  #
  # }else if(priorSigma$cholesky_priorU == "normal"){
  #
  #   priorSigma_in$cholesky_V_i <- rep(priorSigma$cholesky_V_i, length = n_L)
  #
  # }else if(priorSigma$cholesky_priorU == "HMP"){
  #   priorSigma_in$cholesky_lambda_3 <- priorSigma$cholesky_lambda_3
  # }
  # else if(priorSigma$cholesky_priorU == "R2D2"){
  #
  #   if(is.numeric(priorSigma$cholesky_R2D2_b)){
  #     priorSigma_in$cholesky_R2D2_hyper <- FALSE
  #   }else if(priorSigma$cholesky_R2D2_b == "hyperprior"){
  #
  #     priorSigma_in$cholesky_R2D2_hyper <- TRUE
  #
  #   }else{
  #     stop("Something went wrong specifying R2D2_b (specify_L...)!")
  #   }
  #
  #   if(!priorSigma_in$cholesky_R2D2_hyper){
  #
  #     priorSigma_in$cholesky_R2D2_b <- priorSigma$cholesky_R2D2_b
  #     priorSigma_in$cholesky_R2D2_api <- 1/(n_L^(priorSigma_in$cholesky_R2D2_b/2) *
  #                                             (T^(priorSigma_in$cholesky_R2D2_b/2)) *log(T))
  #
  #
  #   }else if(priorSigma_in$cholesky_R2D2_hyper){
  #
  #     grid <- 100
  #     priorSigma_in$cholesky_b_vec <- seq(.01,1,length.out=grid)
  #     priorSigma_in$cholesky_api_vec <- rep(as.numeric(NA), grid)
  #     for(j in seq.int(grid)){
  #       priorSigma_in$cholesky_api_vec[j] <- 1/(n_L^(priorSigma_in$cholesky_b_vec[j]/2) *
  #                                                 (T^(priorSigma_in$cholesky_b_vec[j]/2)) *log(T))
  #     }
  #     priorSigma_in$cholesky_R2D2_b <- 0.5
  #     priorSigma_in$cholesky_R2D2_api <- priorSigma_in$cholesky_api_vec[which(priorSigma_in$cholesky_b_vec==0.5)]
  #
  #   }
  # }


# startvals_in ---------------------------------------------------------------

  startvals_in <- list(
    PHI = matrix(0,1,1),
    L = matrix(0,1,1),
    sv_para = matrix(0,1,1),
    sv_logvar = matrix(0,1,1),
    sv_logvar0 = numeric(1L),
    factor_startval = list(facload = matrix(0,1,1), fac = matrix(0,1,1), tau2 = matrix(0,1,1))
  )

  startvals_in$PHI <- if(is.null(PHI_in)) PHI_flat else PHI_in
  startvals_in$L <- if(is.null(L_in)) L_flat else L_in

  if(priorSigma$type == "cholesky"){
    startvals_in$sv_para <- matrix(data= c(rep(-10,M), rep(0.9,M), rep(0.2,M)), #, rep(-10,M)
                                         nrow = 3, ncol = M, byrow = TRUE)
    startvals_in$sv_logvar0 <- rep(-10,M)
    startvals_in$sv_logvar <- matrix(rep(-10, T*M), T,M)
  }else if(priorSigma$type == "factor"){
    startfacload <- matrix(rnorm(M*priorSigma$factor_factors, sd = .5)^2, nrow=M, ncol=priorSigma$factor_factors)
    startfac <- matrix(rnorm(priorSigma$factor_factors*T, 0, sd=.1), nrow=priorSigma$factor_factors)
    startpara <- rbind(mu = c(rep(-3, M) + rnorm(M), rep(0, priorSigma$factor_factors)),
                      phi = c(rep(.8, M), rep(.8, priorSigma$factor_factors)) + pmin(rnorm(M + priorSigma$factor_factors, sd=.06), .095),
                      sigma = rep(.1, M + priorSigma$factor_factors) + rgamma(M + priorSigma$factor_factors, 1, 10))
    startlogvar <- matrix(startpara["mu",][1] + rnorm(T*(M + priorSigma$factor_factors)), T, M + priorSigma$factor_factors)
    startlogvar[,M+which(isFALSE(priorSigma$sv_heteroscedastic[-c(1:M)]))] <- 0 # !!! important, needed for factorstochvol: if factor is assumed to be homoscedastic, the corresponding column in logvar has to be 0!!!
    startlogvar0 <- startpara["mu",][1] + rnorm(M + priorSigma$factor_factors)
    starttau2 <- if(!priorSigma$factor_ngprior){ # if prior is 'normal'
      priorSigma$factor_starttau2
    }else{
      matrix(1, nrow = M, ncol = priorSigma$factor_factors)
    }
    startvals_in$factor_startval <- list(facload = startfacload,
                                          fac = startfac,
                                          tau2 = starttau2)
    startvals_in$sv_para <- startpara
    startvals_in$sv_logvar0 <- startlogvar0
    startvals_in$sv_logvar <- startlogvar
  }

  ######### expert settings
  if(!is.null(expert)){
    if(exists("huge", expert)){
      expert$huge <- expert$huge
    }else{
      expert$huge <- FALSE
    }
  }

  res <- bvar_cpp(Y,
                  X,
                  M,
                  T,
                  K,
                  draws,
                  burnin,
                  thin,
                  tvp_keep,
                  intercept,
                  priorIntercept,
                  #PHI_in, to starvals_in
                  PHI0,
                  priorPHI_in,
                  priorSigma,# priorL_in
                  startvals_in,
                  #L_in, to starvals_in
                  # cholesky_SV, # SV
                  # cholesky_priorHomoscedastic, #priorHomoscedastic
                  # cholesky_sv_spec, #sv_spec
                  # cholesky_h_init, #h_init to starvals_in
                  # cholesky_sv_para_init, #sv_para_init to starvals_in
                  i_mat,
                  i_vec,
                  progressbar,
                  PHI_tol,
                  L_tol,
                  expert$huge
                  )

  #Rcpp timer is in nanoseconds
  #conversion to secs per iteration
  res$bench <- diff(res$bench)/(10^(9)*(draws+burnin))
  attributes(res$bench) <- list("names" = "secs/itr")
  dimnames(res$PHI)[2] <- list(colnames(X))
  dimnames(res$PHI)[3] <- list(colnames(Y))
  phinames <- as.vector((vapply(seq_len(M), function(i) paste0(colnames(Y)[i], "~", colnames(X[,1:(ncol(X)-intercept)])), character(K))))
  if(priorPHI$prior %in% c("DL","DL_h")){
    # if(priorPHI_in$DL_method==1L){
    #   colnames(res$phi_hyperparameter) <- c(paste0("a",1:priorPHI_in$n_groups), paste0("psi: ", phinames), paste0("theta: ", phinames), paste0("zeta",1:priorPHI_in$n_groups))#
    # }else if(priorPHI_in$DL_method==2L){
      colnames(res$phi_hyperparameter) <- c(paste0("a",1:priorPHI_in$n_groups), paste0("psi: ", phinames), paste0("lambda: ", phinames), paste0("xi",1:priorPHI_in$n_groups))
    # }

  }else if(priorPHI$prior == "GT"){

    colnames(res$phi_hyperparameter) <- c(paste0("a",1:priorPHI_in$n_groups),
                                          paste0("xi",1:priorPHI_in$n_groups),
                                          paste0("psi: ", phinames),
                                          paste0("lambda: ", phinames))

  }else if(priorPHI$prior == "R2D2"){
    colnames(res$phi_hyperparameter) <- c(paste0("zeta",1:priorPHI_in$n_groups),
                                          paste0("psi: ", phinames),
                                          paste0("theta: ", phinames),
                                          paste0("xi",1:priorPHI_in$n_groups),
                                          paste0("b",1:priorPHI_in$n_groups),
                                          paste0("api",1:priorPHI_in$n_groups))#

  }else if(priorPHI$prior == "HS"){
    colnames(res$phi_hyperparameter) <- c(paste0("zeta", 1:priorPHI_in$n_groups),
                                          paste0("theta: ", phinames),
                                          paste0("varpi", 1:priorPHI_in$n_groups),
                                          paste0("nu: ", phinames))
  }else if(priorPHI$prior == "NG"){
    colnames(res$phi_hyperparameter) <- c(paste0("zeta", 1:priorPHI_in$n_groups),
                                          paste0("a",1:priorPHI_in$n_groups),
                                          paste0("theta_tilde: ", phinames))

  }else if(priorPHI$prior == "SSVS"){
  colnames(res$phi_hyperparameter) <- c(paste0("gamma: ", phinames), paste0("p_i: ", phinames))#

  }else if(priorPHI$prior == "HMP"){
    colnames(res$phi_hyperparameter) <- c("lambda_1", "lambda_2")
  }

  if(priorSigma$type=="cholesky"){
    dimnames(res$L)[2] <- dimnames(res$L)[3] <- list(colnames(Y))
    lnames <- NULL
    for(j in 2:M){
      lnames <- c(lnames, paste0(colnames(Y)[j],"~", colnames(Y)[1:(j-1)]))
    }
    if(priorSigma$cholesky_priorU %in% c("DL","DL_h")){
      colnames(res$l_hyperparameter) <- c("a", paste0("psi: ", lnames), paste0("lambda: ", lnames), "xi")
    }else if(priorSigma$cholesky_priorU == "GT"){

      colnames(res$l_hyperparameter) <- c("a","xi",
                                          paste0("psi: ", lnames),
                                          paste0("lambda: ", lnames))

    }else if(priorSigma$cholesky_priorU == "HS"){
      colnames(res$l_hyperparameter) <- c("zeta",
                                          "varpi",
                                          paste0("theta: ", lnames),
                                          paste0("nu: ", lnames))
    }else if(priorSigma$cholesky_priorU == "R2D2"){
      colnames(res$l_hyperparameter) <- c("zeta", paste0("psi: ", lnames), paste0("theta: ", lnames), "xi", "b", "api")#

    }else if(priorSigma$cholesky_priorU == "SSVS"){
      colnames(res$l_hyperparameter) <- c(paste0("gamma: ", lnames), paste0("p_i: ", lnames))
    }else if(priorSigma$cholesky_priorU == "HMP"){
      colnames(res$l_hyperparameter) <- c("lambda_3")
    }
    res$l_hyperparameter <- as.data.frame(res$l_hyperparameter)
    class(res$L) <- "L"
  }

  res$phi_hyperparameter <- as.data.frame(res$phi_hyperparameter)

  res$Y <- Y
  res$X <- X
  res$p <- p
  res$intercept <- ifelse(intercept>0, TRUE, FALSE)
  #res$SV <- if(all(priorSigma$sv_heteroscedastic==TRUE)) TRUE else FALSE
  res$heteroscedastic <- priorSigma$sv_heteroscedastic #if(all(priorSigma$cholesky_heteroscedastic==TRUE) || all(priorSigma$factor_heteroskedastic==TRUE)) TRUE else FALSE
  res$Yraw <- Y_tmp
  res$Traw <- Traw
  res$Sigma_type <- priorSigma$type
  res$datamat <- data.frame(cbind(Y, X))
  class(res$PHI) <- "PHI"
  class(res) <- "bvar"
  res
}

#
# bvar_factor <- function(Yraw,
#                      p,
#                      intercept = 100,
#                      draws=1000,
#                      burnin=1000,
#                      persistence = 0,
#                      priorPHI,
#                      priorSigma,
#                      progressbar=TRUE,
#                      PHI_in=NULL,
#                      L_in=NULL,
#                      PHI_tol = 1e-100,
#                      L_tol = 1e-100
# ){
#
#   # Data preliminaries ------------------------------------------------------
#
#   # M: number of variables,
#   # T: number of observations used for estimation,
#   # K: number of covariates per equation (without intercepts!!!),
#   # n: number of VAR coefficients (without intercepts!!!),
#   # n_L: number of free off-diagonal elements in L
#   # Y: Yraw without first p observations,
#   # X: lagged values of Yraw
#
#   p <- as.integer(p)
#   M <- ncol(Yraw)
#   K <- p*M
#   Traw <- nrow(Yraw)
#   Y_tmp <- as.matrix(Yraw)
#   if (any(is.na(Y_tmp))){
#     stop("\nNAs in Yraw.\n")
#   }
#   if (ncol(Y_tmp) < 2) {
#     stop("The matrix 'Yraw' should contain at least two variables. \n")
#   }
#   if (is.null(colnames(Y_tmp))) {
#     colnames(Y_tmp) <- paste("y", 1:ncol(Y_tmp), sep = "")
#     warning(paste("No column names supplied in Yraw, using:",
#                   paste(colnames(Y_tmp), collapse = ", "), ", instead.\n"))
#   }
#   colnames(Y_tmp) <- make.names(colnames(Y_tmp))
#   # embed: get lagged values
#   X <- embed(Y_tmp, dimension = p + 1)[, -(1:M)]
#   if(p>0){
#     colnames(X) <- paste0(colnames(Y_tmp), ".l", sort(rep(1:p,M)))
#   }
#   if(is.numeric(intercept)){
#     X <- cbind(X,1)
#     colnames(X)[ncol(X)] <- c("intercept")
#   }
#
#   Y <- if(p>0) Y_tmp[-c(1:p), ] else if(p==0) Y_tmp
#
#   T <- Traw - p
#   #  if(T!=nrow(Y) | T!=nrow(X)){
#   #    stop("Something went wrong: T != nrow(Y). \n")
#   #  }
#   n <- K*M
#   n_L <- (M^2 - M)/2
#
#   # Indicator matrix --------------------------------------------------------
#
#   if(is.numeric(intercept)){
#     i_intercept <- rep(0,M)
#   }else i_intercept <- NULL
#
#   if(priorPHI$prior == "DL" | priorPHI$prior == "R2D2" |
#      priorPHI$prior == "HS" | priorPHI$prior == "SSVS" |
#      priorPHI$prior == "NG" | priorPHI$prior == "GT"){
#
#     if(all(is.numeric(priorPHI$global_grouping))){
#       i_mat <- priorPHI$global_grouping
#       if(!(identical(dim(i_mat), as.integer(c(K,M))) & all(is.numeric(i_mat))) )
#         stop("Something went wrong specifying 'global_grouping'.")
#     }else if(priorPHI$global_grouping=="global"){
#       i_mat <- matrix(1, K, M)
#     }else if(priorPHI$global_grouping=="fol"){
#       i_mat <- matrix(1, K, M)
#       diag(i_mat[1:M,1:M]) <- 2
#     }else if(priorPHI$global_grouping=="olcl-lagwise"){
#       i_mat <- matrix(1, K, M)
#       diag(i_mat[1:M,1:M]) <- 2
#       if(p>1){
#         for (j in 1:(p-1)) {
#           i_mat[(j*M+1):((j+1)*M),] <- i_mat[((j-1)*M+1):(j*M),] + 2
#         }
#       }
#     }else if(priorPHI$global_grouping == "equation-wise"){
#       i_mat <- matrix(
#         rep(1:M, K),
#         K,M,
#         byrow = TRUE
#       )
#     }else if(priorPHI$global_grouping == "covariate-wise"){
#       i_mat <- matrix(
#         rep(1:K, M),
#         K,M
#       )
#     }else {
#       stop("Something went wrong specifying 'global_grouping'.")
#     }
#
#   }else{
#     i_mat_1 <- diag(M)
#     i_mat_1[upper.tri(i_mat_1)] <-
#       i_mat_1[lower.tri(i_mat_1)] <- -1
#     i_mat <- NULL
#     for (j in seq_len(p)) {
#       i_mat <- rbind(i_mat, i_mat_1 * j)
#     }
#   }
#
#   i_mat <- rbind(i_mat, i_intercept)
#   i_vec <- as.vector(i_mat)
#
#   # Hyperparameter settigns -------------------------------------------------
#
#   if(!intercept){
#     intercept <- 0
#     priorIntercept <- vector("numeric")
#   }else{
#     ##some checks???
#     priorIntercept <- rep(intercept, M)
#     intercept <- 1
#   }
#
#   if(persistence == 0) {
#
#     PHI0 <- matrix(0, K+intercept, M)
#
#   }else {
#
#     PHI0 <- matrix(0, K+intercept, M)
#     PHI0[1:M, 1:M] <- diag(M)*persistence
#
#   }
#
#   if(!(priorPHI$prior %in% c("DL", "HS","NG", "HMP", "SSVS", "normal", "R2D2", "GT"))){
#     stop("Argument 'priorPHI$prior' must be one of
#            'DL', 'R2D2', 'HS', 'NG', 'SSVS', 'HMP' or 'normal'. \n")
#   }
#
#   # Initialize --------------------------------------------------------------
#
#   if(priorPHI$prior == "SSVS" | is.null(PHI_in) | is.null(L_in)){
#     # Posterior mean of a flat conjugate Normal inverse Wishart prior
#     # exists even when OLS estimate does not exist (in situations where T < K)
#     # N(0, 10^3) on PHI, and invWish(I, M+2) on Sigma
#     XX <- crossprod(X)
#     V_post_flat <- chol2inv(chol(diag(1/rep(10^3, (K+intercept)), nrow = K + intercept) + XX))
#     PHI_flat <- V_post_flat %*% (diag(1/rep(10^3, (K+intercept)), nrow = K + intercept)%*%PHI0 + crossprod(X,Y))
#     S_post <- diag(M) + crossprod(Y - X%*%PHI_flat) + t(PHI_flat - PHI0) %*%
#       diag(1/rep(10^3, (K+intercept)), nrow = K + intercept) %*% (PHI_flat - PHI0)
#     Sigma_flat <- (S_post)/(M +2 + T - M - 1)
#     U <- chol(Sigma_flat)
#     D <- diag(U)^2
#     L_inv <- U/sqrt(D)
#     L_flat <- backsolve(L_inv, diag(M))
#   }
#
#   #some proper checks missing!!!
#   if(is.null(PHI_in)){
#     PHI_in <- PHI_flat
#   }
#   if(is.null(L_in)){
#     L_in <- L_flat
#   }
#
#   # Hyperparameter settings w.r.t PHI ---------------------------------------
#
#   # creating placeholders (for cpp, maybe move to cpp code)
#   priorPHI_in <- list()
#   priorPHI_in$prior <- priorPHI$prior
#
#   #GL priors
#   priorPHI_in$n_groups <- 1
#   priorPHI_in$groups <- 1
#   priorPHI_in$GL_tol <- double(1L)
#   priorPHI_in$a <- double(1L)
#   priorPHI_in$b <- double(1L)
#   priorPHI_in$c <- double(1L)
#   priorPHI_in$GT_vs <- double(1L)
#   priorPHI_in$GT_priorkernel <- character(1L)
#   priorPHI_in$a_vec <- double(1)
#   priorPHI_in$a_weight <- double(1)
#   priorPHI_in$norm_consts <- double(1)
#   priorPHI_in$c_vec <- double(1)
#   priorPHI_in$c_rel_a <- logical(1L)
#
#   #DL
#   priorPHI_in$DL_method <- integer(1)
#   priorPHI_in$prep2 <- matrix(0)
#   priorPHI_in$prep1 <- double(1)
#   priorPHI_in$DL_hyper <- logical(1)
#   priorPHI_in$DL_plus <- logical(1)
#   #GT (Normal Gamma and R2D2 belong to GT)
#   priorPHI_in$GT_hyper <- logical(1)
#   #SSVS
#   priorPHI_in$SSVS_tau0 <- double(1)
#   priorPHI_in$SSVS_tau1 <- double(1)
#   priorPHI_in$SSVS_s_a <- double(1)
#   priorPHI_in$SSVS_s_b <- double(1)
#   priorPHI_in$SSVS_p <- double(1)
#   priorPHI_in$SSVS_hyper <- logical(1)
#   #HM
#   priorPHI_in$lambda_1 <- double(1)
#   priorPHI_in$lambda_2 <- double(1)
#   priorPHI_in$V_i_prep <- double(1)
#
#   # prior specification for PHI
#
#   if(priorPHI$prior == "R2D2" | priorPHI$prior == "DL" |
#      priorPHI$prior == "DL_h" | priorPHI$prior == "SSVS" |
#      priorPHI$prior == "HS" | priorPHI$prior == "NG" | priorPHI$prior == "GT"){
#
#     groups <- unique(i_vec[i_vec!=0])
#     n_groups <- length(groups)
#     priorPHI_in$n_groups <- n_groups
#     priorPHI_in$groups <- groups
#
#   }
#
#   if(priorPHI$prior == "GT"){
#
#     priorPHI_in$GT_priorkernel <- priorPHI$GT_priorkernel
#     priorPHI_in$GL_tol <- priorPHI$GL_tol
#     priorPHI_in$GT_vs <- priorPHI$GT_vs
#     priorPHI_in$b <- rep_len(priorPHI$b, n_groups)
#
#     if(is.matrix(priorPHI$a)){
#       if(ncol(priorPHI$a)==2){
#         priorPHI_in$GT_hyper <- TRUE
#         priorPHI_in$a_vec <- priorPHI$a[,1]
#         priorPHI_in$a_weight <- priorPHI$a[,2]
#         priorPHI_in$norm_consts <- lgamma(priorPHI_in$a_vec)
#         priorPHI_in$a <- sample(priorPHI_in$a_vec, n_groups, replace = TRUE, prob = priorPHI_in$a_weight) # initialize a
#       }else if(ncol(priorPHI$a)>2){
#         stop("The easiest way to specify 'R2D2_a', 'NG_a' or 'GT_a' is a single postive number!")
#       }else{
#         priorPHI$a <- as.vector(priorPHI$a)
#       }
#     }
#     if(is.null(dim(priorPHI$a))){
#       priorPHI_in$GT_hyper <- FALSE
#       priorPHI_in$a <- rep_len(priorPHI$a, n_groups)
#     }
#
#     if(is.character(priorPHI$c)){
#       priorPHI_in$c_rel_a <- TRUE # then c is always proportion of a (e.g. for R2D2 c=0.5a)
#       mya <- priorPHI_in$a
#       myc <- gsub("a","mya", priorPHI$c)
#       priorPHI_in$c <- eval(str2lang(myc))
#       if(base::isTRUE(priorPHI_in$GT_hyper)){
#         myc2 <- gsub("a","priorPHI_in$a_vec", priorPHI$c)
#         priorPHI_in$c_vec <- eval(str2lang(myc2))
#       }
#     }else if(is.numeric(priorPHI$c)){
#       priorPHI_in$c_rel_a <- FALSE
#       priorPHI_in$c <- rep_len(priorPHI$c, n_groups)
#     }
#
#   }else if(priorPHI$prior == "R2D2"){
#
#     if(!exists("R2D2_method", priorPHI)){
#       priorPHI_in$R2D2_method <- 1L
#     }else{
#       priorPHI_in$R2D2_method <- priorPHI$R2D2_method
#     }
#
#     if(!exists("R2D2_kernel", priorPHI)){
#       priorPHI_in$R2D2_kernel <- "laplace"
#     }else{
#       priorPHI_in$R2D2_kernel <- priorPHI$R2D2_kernel
#     }
#
#     if(all(is.numeric(priorPHI$R2D2_b))){
#       R2D2_hyper <- FALSE
#     }else if("hyperprior" %in% priorPHI$R2D2_b & length(priorPHI$R2D2_b) == 1){
#
#       R2D2_hyper <- TRUE
#
#     }else{
#       stop("Something went wrong specifying R2D2_b!")
#     }
#     priorPHI_in$R2D2_hyper <- R2D2_hyper
#
#     if(!R2D2_hyper){
#
#       priorPHI_in$R2D2_b <- rep_len(priorPHI$R2D2_b, n_groups)
#       if(priorPHI$R2D2_api=="default"){
#         priorPHI_in$R2D2_api <- rep(as.numeric(NA), n_groups)
#         for(j in seq.int(n_groups)){
#           priorPHI_in$R2D2_api[j] <- 1/(n^(priorPHI_in$R2D2_b[j]/2) *
#                                           (T^(priorPHI_in$R2D2_b[j]/2)) *log(T))
#         }
#       }else{
#         priorPHI_in$R2D2_api <- priorPHI$R2D2_api
#       }
#
#
#     }else if(R2D2_hyper){
#
#       grid <- 100
#       priorPHI_in$b_vec <- seq(.01,1,length.out=grid)
#       if(priorPHI$R2D2_api=="default"){
#         priorPHI_in$api_vec <- rep(as.numeric(NA), grid)
#         for(j in seq.int(grid)){
#           priorPHI_in$api_vec[j] <- 1/(n^(priorPHI_in$b_vec[j]/2) *
#                                          (T^(priorPHI_in$b_vec[j]/2)) *log(T))
#         }
#       }else{
#         priorPHI_in$api_vec <- rep_len(priorPHI$R2D2_api, grid)
#       }
#
#       #initialize
#       priorPHI_in$R2D2_b <- rep(0.5, n_groups)
#       priorPHI_in$R2D2_api <- rep(priorPHI_in$api_vec[which(priorPHI_in$b_vec==0.5)],
#                                   n_groups)
#     }
#
#
#
#   }else if(priorPHI$prior == "DL"){
#
#     priorPHI_in$GL_tol <- priorPHI$GL_tol
#     if(!exists("DL_method", priorPHI)){
#       priorPHI_in$DL_method <- 2L
#     }else{
#       priorPHI_in$DL_method <- priorPHI$DL_method
#     }
#
#     if(all(is.numeric(priorPHI$a))&is.null(dim(priorPHI$a))){
#       DL_hyper <- FALSE
#     }else if(("hyperprior" %in% priorPHI$a & length(priorPHI$a) == 1) ||
#              (is.matrix(priorPHI$a))){
#
#       DL_hyper <- TRUE
#     }else if(any(c("1/K", "1/n") %in% priorPHI$a) &
#              length(priorPHI$a) == 1){
#       DL_hyper <- FALSE
#       if(priorPHI$a == "1/K") {
#         priorPHI_in$a <- rep(1/K, n_groups)
#       }else if(priorPHI$a == "1/n") {
#         priorPHI_in$a <- rep(1/n, n_groups)
#       }
#     }else{
#       stop("Something went wrong specifying DL_a!")
#     }
#     priorPHI_in$DL_hyper <- DL_hyper
#     if(DL_hyper){
#       if(is.matrix(priorPHI$a)){
#         priorPHI_in$a_vec <- priorPHI$a[,1]
#         priorPHI_in$a_weight <- priorPHI$a[,2]
#         priorPHI_in$norm_consts <- lgamma(priorPHI_in$a_vec)
#         priorPHI_in$a <- rep_len(priorPHI_in$a_vec[1], n_groups) #initial value
#         priorPHI_in$DL_method <- 2L
#       }else{
#         grid <- 1000
#         priorPHI_in$a_vec <- seq(1/(n),1/2,length.out = grid)
#         # prep1 &prep2: some preparations for the evaluation of the
#         # Dirichlet-density
#         priorPHI_in$prep1 <- priorPHI_in$a_vec-1
#         prep2 <- matrix(NA, n_groups, grid)
#         ii <- 1
#         for(j in groups){
#           n_tmp <- length(which(i_vec == j))
#           # normalizing constant of symmetric Dirichlet density in logs
#           prep2[ii,] <- lgamma(n_tmp*priorPHI_in$a_vec) - n_tmp*lgamma(priorPHI_in$a_vec)
#           ii <- ii + 1
#         }
#         priorPHI_in$prep2 <- prep2
#         priorPHI_in$a <- rep_len(1/2, n_groups) #initial value
#         if(priorPHI_in$DL_method==2L){
#           priorPHI_in$a_weight <- rep(1,grid)
#           priorPHI_in$norm_consts <- lgamma(priorPHI_in$a_vec)
#         }
#       }
#
#     }else if(!DL_hyper){
#       if(all(is.numeric(priorPHI$a))){
#         priorPHI_in$a <- rep_len(priorPHI$a, n_groups)
#       }
#     }
#
#     if(!exists("DL_plus", priorPHI) || base::isFALSE(priorPHI$DL_plus)){
#       priorPHI_in$DL_plus <- FALSE
#     }else if(priorPHI$DL_plus){
#       priorPHI_in$DL_plus <- TRUE
#       if(!exists("DL_b", priorPHI)){
#         priorPHI_in$b <- rep_len(0.5, n_groups)
#       }else{
#         priorPHI_in$b <- rep_len(priorPHI$DL_b, n_groups)
#       }
#       if(!exists("DL_c", priorPHI)){
#         priorPHI_in$c <- 0.5*priorPHI_in$a
#         if(DL_hyper){
#           priorPHI_in$c_vec <- 0.5*priorPHI_in$a_vec
#           priorPHI_in$c_rel_a <- TRUE
#         }
#       }else{
#         if(!is.numeric(priorPHI$DL_c)){
#           stop("If you specify DL_c for model DL_plus, then it must be numeric!")
#         }
#         priorPHI_in$c <- rep_len(priorPHI$DL_c, n_groups)
#         priorPHI_in$c_rel_a <- FALSE
#
#       }
#       priorPHI_in$DL_method <- 2L
#     }else{
#       stop("Never heard of DL_plus?")
#     }
#
#   }else if(priorPHI$prior == "NG"){
#     #priorPHI_in$NG_a <- priorPHI$NG_a
#     #priorPHI_in$NG_varrho0 <- priorPHI$NG_varrho0
#     #priorPHI_in$NG_varrho1 <- priorPHI$NG_varrho1
#
#     if(all(is.numeric(priorPHI$NG_a))){
#       NG_hyper <- FALSE
#     }else if(priorPHI$NG_a == "hyperprior"){
#       NG_hyper <- TRUE
#     }else{
#       stop("Something went wrong specifying NG_a!")
#     }
#     priorPHI_in$NG_hyper <- NG_hyper
#     if(NG_hyper){
#       priorPHI_in$NG_a_vec <- seq(0.01,1,length.out=100)
#       priorPHI_in$NG_a <- rep_len(0.1, n_groups)
#     }else{
#       priorPHI_in$NG_a <- rep_len(priorPHI$NG_a, n_groups)
#     }
#
#   }else if(priorPHI$prior == "SSVS"){
#     priorPHI_in$SSVS_hyper <- priorPHI$SSVS_hyper
#     priorPHI_in$SSVS_p <- rep_len(priorPHI$SSVS_p,n)
#     if(priorPHI$semiautomatic) {
#       # shape parameters of beta hyperprior on prior-inclusion-probability
#       priorPHI_in$SSVS_s_a <- priorPHI$SSVS_s_a
#       priorPHI_in$SSVS_s_b <- priorPHI$SSVS_s_b
#
#       # standard errors of flat phi posterior estimates
#       sigma_phi <- sqrt(diag(Sigma_flat %x% V_post_flat))
#       # select all, except those of intercept
#       sigma_phi <- sigma_phi[which(i_vec!=0)]
#       if(length(sigma_phi)!=p*M^2){
#         stop(length(sigma_phi))
#       }
#
#       # scale tau with variances
#       priorPHI_in$SSVS_tau0 <- priorPHI$SSVS_c0*sigma_phi
#       priorPHI_in$SSVS_tau1 <- priorPHI$SSVS_c1*sigma_phi
#     }else{
#       priorPHI_in$SSVS_tau0 <- rep_len(priorPHI$SSVS_c0, n)
#       priorPHI_in$SSVS_tau1 <- rep_len(priorPHI$SSVS_c1, n)
#     }
#
#   }else if(priorPHI$prior == "normal"){
#
#     priorPHI_in$V_i <- rep(priorPHI$V_i, length = n)
#
#   }else if(priorPHI$prior == "HMP"){
#
#     priorPHI_in$lambda_1 <- priorPHI$lambda_1
#     priorPHI_in$lamda_2 <- priorPHI$lambda_2
#     sigma_sq <- MP_sigma_sq(Yraw, 6)
#
#     # prepare prior variances down to lambdas
#     priorPHI_in$V_i_prep <- MP_V_prior_prep(sigma_sq, (K+intercept), M, intercept>0)
#   }
#
#   # Hyperparameter settings w.r.t Sigma -----------------------------------------
#
#
#   ### priorSigma
#   priorSigma_in <- priorSigma
#
#   ### SV settings
#   if(priorSigma_in$type == "factor"){
#
#     if(!(length(priorSigma_in$factor_priorh0idi) == 1 | length(priorSigma_in$factor_priorh0idi) == M)){
#       stop("Argument 'factor_priorh0idi' must be either of length 1 or 'ncol(data)'.")
#     }
#     priorSigma_in$factor_priorh0idi <- rep_len(priorSigma_in$factor_priorh0idi,M)
#     priorSigma_in$factor_priorh0 <- c(priorSigma_in$factor_priorh0idi, priorSigma_in$factor_priorh0fac)
#
#     if(!(length(priorSigma_in$factor_heteroskedastic)==1 | length(priorSigma_in$factor_heteroskedastic) == 2 | length(priorSigma_in$factor_heteroskedastic) == (M + priorSigma_in$factor_factors))){
#       stop("Argument 'factor_heteroskedastic' must be of length 1, 2, or 'ncol(data) + factor_factors.'")
#     }
#     if(length(priorSigma_in$factor_heteroskedastic)==1){
#       priorSigma_in$factor_heteroskedastic <- rep(priorSigma_in$factor_heteroskedastic, M + priorSigma_in$factor_factors)
#     }else if(length(priorSigma_in$factor_heteroskedastic)==2){
#       priorSigma_in$factor_heteroskedastic <- c(rep(priorSigma_in$factor_heteroskedastic[1], M), rep(priorSigma_in$factor_heteroskedastic[2], priorSigma_in$factor_factors))
#     }
#     if (!all(priorSigma_in$factor_heteroskedastic[M+seq_len(priorSigma_in$factor_factors)])) {
#       if (priorSigma_in$factor_interweaving == 2L || priorSigma_in$factor_interweaving == 4L) {
#         warning("Cannot do deep interweaving if (some) factors are homoskedastic. Setting 'interweaving' to 3.")
#         priorSigma_in$factor_interweaving <- 3L
#       }
#     }
#
#     if(!all(priorSigma_in$factor_priorheteroskedastic)){
#       if(nrow(priorSigma_in$factor_priorhomoskedastic)!=M & nrow(priorSigma_in$factor_priorhomoskedastic) > 1L){
#         stop("Argument 'priorhomoskedastic' must be either  of length two or a matrix with positive entries and dimension 'number of time series times two'.")
#       }
#       if(ncol(priorSigma_in$factor_priorhomoskedastic)!=2){
#         stop("Argument 'priorhomoskedastic' must be either  of length two or a matrix with positive entries and dimension 'number of time series times two'.")
#       }
#       if(nrow(priorSigma_in$factor_priorhomoskedastic) == 1L){
#         priorSigma_in$factor_priorhomoskedastic <- matrix(priorSigma_in$factor_priorhomoskedastic, M, 2, byrow = TRUE)
#       }
#     }
#     ### To do priorhomoskedastic
#
#
#
#     if(!(length(priorSigma_in$factor_priorsigmaidi)==1 | length(priorSigma_in$factor_priorsigmaidi) == M)){
#       stop("Argument 'factor_priorsigmaidi' must be either of length 1 or ncol(data).")
#     }
#     priorSigma_in$factor_priorsigmaidi <- rep_len(priorSigma_in$factor_priorsigmaidi, M)
#
#     priorSigma_in$factor_priorsigma <- c(priorSigma_in$factor_priorsigmaidi, priorSigma_in$factor_priorsigmafac)
#
#     cShrink <- priorSigma_in$factor_priorng[1]
#     dShrink <- priorSigma_in$factor_priorng[2]
#
#     if(is.matrix(priorSigma_in$factor_priorfacload)) {
#       if(nrow(priorSigma_in$factor_priorfacload) != M || ncol(priorSigma_in$factor_priorfacload) != priorSigma_in$factor_factors) {
#         stop("If argument 'priorfacload' is a matrix, it must be of appropriate dimensions.")
#       }
#       if (priorSigma_in$factor_priorfacloadtype == "normal") {
#         priorSigma_in$factor_pfl <- 1L
#         priorSigma_in$factor_starttau2 <- priorSigma_in$factor_priorfacload^2
#         aShrink <- NA
#         cShrink <- NA
#         dShrink <- NA
#       } else if (priorSigma_in$factor_priorfacloadtype == "rowwiseng") {
#         priorSigma_in$factor_pfl <- 2L
#         priorSigma_in$factor_starttau2 <- matrix(1, nrow = M, ncol = priorSigma_in$factor_factors)
#         aShrink <- priorSigma_in$factor_priorfacload[,1]
#         warning("Only first column of 'priorfacload' is used.'")
#         cShrink <- rep(cShrink, M)
#         dShrink <- rep(dShrink, M)
#       } else if (priorSigma_in$factor_priorfacloadtype == "colwiseng") {
#         priorSigma_in$factor_pfl <- 3L
#         priorSigma_in$factor_starttau2 <- matrix(1, nrow = M, ncol = priorSigma_in$factor_factors)
#         aShrink <- priorSigma_in$factor_priorfacload[1,]
#         warning("Only first row of 'priorfacload' is used.'")
#         cShrink <- rep(cShrink, priorSigma_in$factor_factors)
#         dShrink <- rep(dShrink, priorSigma_in$factor_factors)
#       } else if (priorSigma_in$factor_priorfacloadtype == "dl") {
#         stop("dlprior for factorloading is not supported by bayesianVARs!")
#         # priorSigma_in$factor_pfl <- 4L
#         # priorSigma_in$factor_starttau2 <- matrix(1, nrow = M, ncol = priorSigma_in$factor_factors)
#         # aShrink <- priorSigma_in$factor_priorfacload[1,1]
#         # warning("Only first element of 'priorfacload' is used.'")
#         # cShrink <- NA
#         # dShrink <- NA
#       }
#     } else {
#       if (length(priorSigma_in$factor_priorfacload) != 1) {
#         stop("If argument 'priorfacload' isn't a matrix, it must be a single number.")
#       }
#       if (priorSigma_in$factor_priorfacloadtype == "normal") {
#         priorSigma_in$factor_pfl <- 1L
#         priorSigma_in$factor_starttau2 <- matrix(priorSigma_in$factor_priorfacload^2, nrow = M, ncol = priorSigma_in$factor_factors)
#         aShrink <- NA
#         cShrink <- NA
#         dShrink <- NA
#       } else if (priorSigma_in$factor_priorfacloadtype == "rowwiseng") {
#         priorSigma_in$factor_pfl <- 2L
#         priorSigma_in$factor_starttau2 <- matrix(1, nrow = M, ncol = priorSigma_in$factor_factors)
#         aShrink <- rep(priorSigma_in$factor_priorfacload, M)
#         cShrink <- rep(cShrink, M)
#         dShrink <- rep(dShrink, M)
#       } else if (priorSigma_in$factor_priorfacloadtype == "colwiseng") {
#         priorSigma_in$factor_pfl <- 3L
#         priorSigma_in$factor_starttau2 <- matrix(1, nrow = M, ncol = priorSigma_in$factor_factors)
#         aShrink <- rep(priorSigma_in$factor_priorfacload, priorSigma_in$factor_factors)
#         cShrink <- rep(cShrink, priorSigma_in$factor_factors)
#         dShrink <- rep(dShrink, priorSigma_in$factor_factors)
#       } else if (priorSigma_in$factor_priorfacloadtype == "dl") {
#         stop("dlprior for factorloading is not supported by bayesianVARs!")
#         # priorSigma_in$factor_pfl <- 4L
#         # priorSigma_in$factor_starttau2 <- matrix(1, nrow = M, ncol = priorSigma_in$factor_factors)
#         # aShrink <- priorSigma_in$factor_priorfacload
#         # cShrink <- NA
#         # dShrink <- NA
#       }
#     }
#     priorSigma_in$factor_shrinkagepriors <- list(a = aShrink,
#                                                  c = cShrink,
#                                                  d = dShrink)
#
#     restr <- matrix(FALSE, nrow = M, ncol = priorSigma_in$factor_factors)
#     if (priorSigma_in$factor_restrict == "upper") restr[upper.tri(restr)] <- TRUE
#
#     if (priorSigma_in$factor_interweaving %in% c(1, 2) && any(diag(restr) == TRUE)) {
#       stop("Setting 'interweaving' to either 1 or 2 and restricting the diagonal elements of the factor loading matrix are not allowed at the same time.")
#     }
#
#     priorSigma_in$factor_restrinv <- matrix(as.integer(!restr), nrow = nrow(restr), ncol = ncol(restr))
#
#
#     startfacload <- matrix(rnorm(M*priorSigma_in$factor_factors, sd = .5)^2, nrow=M, ncol=priorSigma_in$factor_factors)
#     startfac <- matrix(rnorm(priorSigma_in$factor_factors*T, 0, sd=.1), nrow=priorSigma_in$factor_factors)
#     startpara <- list(mu = c(rep(-3, M) + rnorm(M), rep(0, priorSigma_in$factor_factors)),
#                       phi = c(rep(.8, M), rep(.8, priorSigma_in$factor_factors)) + pmin(rnorm(M + priorSigma_in$factor_factors, sd=.06), .095),
#                       sigma = rep(.1, M + priorSigma_in$factor_factors) + rgamma(M + priorSigma_in$factor_factors, 1, 10))
#     startlogvar <- matrix(startpara[["mu"]][1] + rnorm(T*(M + priorSigma_in$factor_factors)), T, M + priorSigma_in$factor_factors)
#     startlogvar0 <- startpara[["mu"]][1] + rnorm(M + priorSigma_in$factor_factors)
#     starttau2 <- if(priorSigma_in$factor_priorfacloadtype == "normal"){
#       priorSigma_in$factor_starttau2
#     }else{
#       matrix(1, nrow = M, ncol = priorSigma_in$factor_factors)
#     }
#     priorSigma_in$factor_startval <- list(facload = startfacload,
#                                           fac = startfac,
#                                           para = startpara,
#                                           latent = startlogvar,
#                                           latent0 = startlogvar0,
#                                           tau2 = starttau2)
#
#   }
#
#   ##maybe
#   ## adjust length of sv_offset (expert option...)
#   ##priorSigma_in$cholesky_sv_spec$sv_offset <- rep_len(priorSigma_in$cholesky_sv_spec$sv_offset, M)
#   ##
#   res <- bvar_fsv(Y,
#                   X,
#                   M,
#                   T,
#                   K,
#                   draws,
#                   burnin,
#                   intercept,
#                   priorIntercept,
#                   PHI_in,
#                   PHI0,
#                   priorPHI_in,
#                   priorSigma_in,
#                   i_mat,
#                   i_vec,
#                   progressbar,
#                   PHI_tol
#   )
#
#   #Rcpp timer is in nanoseconds
#   #conversion to secs per iteration
#   res$bench <- diff(res$bench)/(10^(9)*(draws+burnin))
#   attributes(res$bench) <- list("names" = "secs/itr")
#   if(p>0){
#     dimnames(res$PHI)[2] <- list(colnames(X))
#     dimnames(res$PHI)[3] <- list(colnames(Y))
#     phinames <- as.vector((vapply(seq_len(M), function(i) paste0(colnames(Y)[i], "~", colnames(X[,1:(ncol(X)-intercept)])), character(K))))
#     if(priorPHI$prior %in% c("DL","DL_h")){
#       if(priorPHI_in$DL_method==1L){
#         colnames(res$phi_hyperparameter) <- c(paste0("a",1:priorPHI_in$n_groups), paste0("psi: ", phinames), paste0("theta: ", phinames), paste0("zeta",1:priorPHI_in$n_groups))#
#       }else if(priorPHI_in$DL_method==2L){
#         colnames(res$phi_hyperparameter) <- c(paste0("a",1:priorPHI_in$n_groups), paste0("psi: ", phinames), paste0("lambda: ", phinames), paste0("xi",1:priorPHI_in$n_groups))
#       }
#
#     }else if(priorPHI$prior == "GT"){
#
#       colnames(res$phi_hyperparameter) <- c(paste0("a",1:priorPHI_in$n_groups),
#                                             paste0("xi",1:priorPHI_in$n_groups),
#                                             paste0("psi: ", phinames),
#                                             paste0("lambda: ", phinames))
#
#     }else if(priorPHI$prior == "R2D2"){
#       colnames(res$phi_hyperparameter) <- c(paste0("zeta",1:priorPHI_in$n_groups),
#                                             paste0("psi: ", phinames),
#                                             paste0("theta: ", phinames),
#                                             paste0("xi",1:priorPHI_in$n_groups),
#                                             paste0("b",1:priorPHI_in$n_groups),
#                                             paste0("api",1:priorPHI_in$n_groups))#
#
#     }else if(priorPHI$prior == "HS"){
#       colnames(res$phi_hyperparameter) <- c(paste0("zeta", 1:priorPHI_in$n_groups),
#                                             paste0("theta: ", phinames),
#                                             paste0("varpi", 1:priorPHI_in$n_groups),
#                                             paste0("nu: ", phinames))
#     }else if(priorPHI$prior == "NG"){
#       colnames(res$phi_hyperparameter) <- c(paste0("zeta", 1:priorPHI_in$n_groups),
#                                             paste0("a",1:priorPHI_in$n_groups),
#                                             paste0("theta_tilde: ", phinames))
#
#     }else if(priorPHI$prior == "SSVS"){
#       colnames(res$phi_hyperparameter) <- c(paste0("gamma: ", phinames), paste0("p_i: ", phinames))#
#
#     }else if(priorPHI$prior == "HMP"){
#       colnames(res$phi_hyperparameter) <- c("lambda_1", "lambda_2")
#     }
#   }
#   res$phi_hyperparameter <- as.data.frame(res$phi_hyperparameter)
#   res$Y <- Y
#   res$X <- X
#   res$p <- p
#   res$intercept <- ifelse(intercept>0, TRUE, FALSE)
#   res$Yraw <- Y_tmp
#   res$Traw <- Traw
#   res$datamat <- data.frame(cbind(Y, X))
#   class(res$PHI) <- "PHI"
#   class(res) <- "bvar_fsv"
#   res
# }
#
