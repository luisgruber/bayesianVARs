#' @export
bvar_fast <- function(Yraw,
                      p,
                      intercept = FALSE,
                      standardize = TRUE,
                      persistence = 0,
                      priorPHI,
                      priorL,
                      draws,
                      burnin,
                      SV=TRUE,
                      sv_spec,
                      progressbar=TRUE
                      ){

# Data preliminaries ------------------------------------------------------


  # M: number of variables, T: number of observations used for estimation,
  # K: number of covariates per equation, n: number of VAR coefficients,
  # n_L: number of free off-diagonal elements in L
  # Y: Yraw without first p observations, X: lagged values of Yraw

  M <- ncol(Yraw)
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
  if(intercept){
    cbind(X,1)
  }
  colnames(X) <- paste0(colnames(Y_tmp), ".l", sort(rep(1:p,M)))
  Y <- Y_tmp[-c(1:p), ]

  mu_Y <- colMeans(Y)
  sd_Y <- apply(Y, 2 ,sd)
  if(standardize==TRUE){

    Y <- scale(Y)
    if(intercept){
      X[,-ncol(X)] <- scale(X[,-ncol(X)])
    }else X <- scale(X)

  }

  T <- Traw - p
  K <- ncol(X)
  if(T!=nrow(Y) | T!=nrow(X)){
    stop("Something went wrong: T != nrow(Y). \n")
  }
  n <- K*M
  n_L <- (M^2 - M)/2

# Hyperparameter settigns -------------------------------------------------

  if(persistence == 0) {

    PHI0 <- matrix(0, K, M)

  }else {

    PHI0 <- matrix(0, K, M)
    PHI0[1:M, 1:M] <- diag(M)*persistence

  }


  V_i_L <- rep(1, n_L)

  if(SV == TRUE & is.null(sv_spec)){
    sv_spec <- list(priormu = c(0,100),
                    priorphi = c(20, 1.5),
                    priorsigma2 = c(0.5,0.5)#,
                    #priorh0 = -1 #h0 from stationary distribution
    )
  }


  #specify_priorPHI <- function(prior, DL_a,
  #                             SSVS_c0, SSVS_c1, SSVS_semiautomatic,
  #                             HMP_lambda1, HMP_lambda2,
  #                             V_i = NULL){

  #  if(!(prior %in% c("DL", "HMP", "SSVS", "normal"))){
  #    stop("Argument 'prior' must be one of
  #         'DL', 'SSVS', 'HMP' or 'normal'.")
  #  }

  #}

  if(!(priorPHI$prior %in% c("DL", "HMP", "SSVS", "normal"))){
    stop("Argument 'priorPHI$prior' must be one of
           'DL', 'SSVS', 'HMP' or 'normal'. \n")
  }

  if(is.null(priorPHI$V_i)){
    V_i <- rep(1, n)}
  a <- 0.5

  if(priorPHI$prior == "normal"){
    if(is.null(priorPHI$V_i)){
      warning("Prior variances for VAR coefficients not specified.
              Using 'rep(1,n)', where 'n' is the number of coefficients. \n")
    }else{
      if(length(priorPHI$V_i) > 1 & (length(priorPHI$V_i) < n | length(priorPHI$V_i) > n) ){
        stop("Length of 'priorPHI$V_i' does not match number of VAR coefficients:
        Specify either one prior variance (single numeric value) that will be
             recycled for all coefficients, or a numeric vector which length
             equals the number of VAR coefficients. \n")
      }
      V_i <- rep(priorPHI$V_i, length = n)
      if(any(V_i<=0)){
        stop("Prior variances (elements of 'priorPHI$V_i') must be strictly
             positive. \n")
      }
    }
  }else if(priorPHI$prior == "DL"){
    if(is.numeric(priorPHI$DL_a)){
      a <-  priorPHI$DL_a
    }else if(priorPHI$DL_a == "hyperprior"){


    }

  }else if(priorPHI$prior == "SSVS"){
    c0 <- priorPHI$SSVS_c0
  }

  if(priorL$prior == "normal"){
    if(is.null(priorL$V_i)){
      warning("Prior variances for VAR coefficients not specified.
              Using 'rep(1,n)', where 'n' is the number of coefficients. \n")
    }else{
      if(length(priorL$V_i) > 1 & (length(priorL$V_i) < n_L | length(priorL$V_i) > n_L) ){
        stop("Length of 'priorL$V_i' does not match number of free off-diagonal elements in L:
        Specify either one prior variance (single numeric value) that will be
             recycled for all coefficients, or a numeric vector which length
             equals the number of free off-diagonal elements in L \n")
      }
      V_i_L <- rep(priorL$V_i, length = n_L)
      if(any(V_i_L<=0)){
        stop("Prior variances (elements of 'priorL$V_i') must be strictly
             positive. \n")
      }
    }
  }else if(priorL$prior == "DL"){
    if(is.numeric(priorL$DL_b)){
      b <-  priorL$DL_b
    }else if(priorL$DL_b == "hyperprior"){


    }

  }else if(priorL$prior == "SSVS"){
    c0 <- priorL$SSVS_c0
  }


# Initialize --------------------------------------------------------------

  # Posterior mean of a flat conjugate Normal inverse Wishart prior
  # exists even when OLS estimate does not exist (in situations where T < K)
  # N(0, 10^3) on PHI, and invWish(I, M+2) on Sigma
  XX <- crossprod(X)
  V_post_flat <- tryCatch(solve(diag(1/rep(10^3, K)) + XX),
                          error = function(e) chol2inv(chol(diag(1/rep(10^3, K)) + XX)))
  PHI <- V_post_flat %*% (diag(1/rep(10^3, K))%*%PHI0 + t(X)%*%Y)
  S_post <- diag(M) + crossprod(Y - X%*%PHI) + t(PHI - PHI0) %*%
    diag(1/rep(10^3, K)) %*% (PHI - PHI0)
  Sigma <- (S_post)/(M +2 + T - M - 1)
  U <- chol(Sigma)
  D <- diag(U)^2
  L_inv <- U/sqrt(D)
  L <- backsolve(L_inv, diag(M))


  sv_para_init <- matrix(data= c(rep(-10,M), rep(0.9,M), rep(0.2,M), rep(-10,M)),
                                nrow = 4, ncol = M, byrow = TRUE)
  rownames(sv_para_init) <- c("mu", "phi", "sigma", "h0")

  h_init <- matrix(rep(-10, T*M), T,M)

  res <- bvar_cpp(Y,
                  X,
                  M,
                  T,
                  K,
                  draws,
                  burnin,
                  PHI,
                  PHI0,
                  priorPHI$prior,
                  a,
                  priorL$prior,
                  b,
                  L,
                  V_i,
                  V_i_L,
                  #expert,
                  sv_spec,
                  h_init,
                  sv_para_init,
                  progressbar
                  )

  #Rcpp timer is in nanoseconds
  #conversion to secs per iteration
  res$bench <- diff(res$bench)/(10^(9)*(draws+burnin))
  attributes(res$bench) <- list("names" = "secs/itr")
  dimnames(res$PHI_draws)[2] <- list(colnames(X))
  dimnames(res$PHI_draws)[3] <- list(colnames(Y))
  dimnames(res$L_draws)[2] <- dimnames(res$L_draws)[3] <- list(colnames(Y))
  res$Y <- Y
  res$X <- X

  res
}
