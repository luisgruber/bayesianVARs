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

# Indicator matrix --------------------------------------------------------

  if(intercept){
    i_intercept <- rep(0,M)
  }else i_intercept <- NULL

  i_mat_1 <- diag(M)
  i_mat_1[upper.tri(i_mat_1)] <-
    i_mat_1[lower.tri(i_mat_1)] <- -1
  i_mat <- NULL
  for (j in seq_len(p)) {
    i_mat <- rbind(i_mat, i_mat_1 * j)
  }
  i_mat <- rbind(i_mat, i_intercept)
  i_vec <- as.vector(i_mat)

# Hyperparameter settigns -------------------------------------------------

  if(persistence == 0) {

    PHI0 <- matrix(0, K, M)

  }else {

    PHI0 <- matrix(0, K, M)
    PHI0[1:M, 1:M] <- diag(M)*persistence

  }

  ##V_i_L <- rep(1, n_L) cpp

  if(!(priorPHI$prior %in% c("DL", "DL_h", "HMP", "SSVS", "normal"))){
    stop("Argument 'priorPHI$prior' must be one of
           'DL', 'SSVS', 'HMP' or 'normal'. \n")
  }

  ##if(is.null(priorPHI$V_i)){ cpp
  ##  V_i <- rep(1, n)}
  if(priorPHI$prior == "DL"){
    if(priorPHI$DL_a == "1/K") priorPHI$DL_a <- 1/K
    if(priorPHI$DL_a == "1/n") priorPHI$DL_a <- 1/n

  }else if(priorPHI$prior == "SSVS"){
    priorPHI$SSVS_tau0 <- rep(priorPHI$SSVS_c0, n)
    priorPHI$SSVS_tau1 <- rep(priorPHI$SSVS_c1, n)
  }else if(priorPHI$prior == "normal"){
    priorPHI$V_i <- rep(priorPHI$V_i, length = n)
  }else if(priorPHI$prior == "HMP"){
    sigma_sq <- MP_sigma_sq(Yraw, 6, standardize)
    # prepare prior variances down to lambdas
    priorPHI$V_i_prep <- MP_V_prior_prep(sigma_sq, K, M, intercept)
  }

  if(priorL$prior == "DL"){
    if(priorL$DL_b == "1/n") priorL$DL_b <- 1/n_L
  }else if(priorL$prior == "SSVS"){
    priorL$SSVS_tau0 <- rep(priorL$SSVS_c0, n)
    priorL$SSVS_tau1 <- rep(priorL$SSVS_c1, n)
  }else if(priorL$prior == "normal"){
    priorL$V_i <- rep(priorL$V_i, length = n_L)
  }


  if(SV == TRUE & is.null(sv_spec)){
    sv_spec <- list(priormu = c(0,100),
                    priorphi = c(20, 1.5),
                    priorsigma2 = c(0.5,0.5)#,
                    #priorh0 = -1 #h0 from stationary distribution
    )
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
                  priorPHI,
                  priorL,
                  L,
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
  res$Y <- Y
  res$X <- X

  res
}
