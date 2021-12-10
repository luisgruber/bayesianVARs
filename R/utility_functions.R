# specify Priors ----------------------------------------------------------

#' @export
specify_priorPHI <- function(prior, DL_a = "1/K", R2D2_b = 0.5,
                             SSVS_c0 = 0.01, SSVS_c1 = 100, SSVS_semiautomatic = TRUE, SSVS_sa = 0.5, SSVS_sb = 0.5,
                             HMP_lambda1 = c(0.01,0.01), HMP_lambda2 = c(0.01,0.01),
                             V_i = NULL){
  if(!(prior %in% c("DL", "HMP", "SSVS", "normal", "R2D2"))){
    stop("Argument 'prior' must be one of 'DL', 'SSVS', 'HMP' or 'normal'. \n")
  }

  if(prior == "DL"){
    text <- c("Argument 'DL_a' must be either a single positive numeric or one of 'hyperprior',
           '1/K' or '1/n'. \n ")
    if(is.numeric(DL_a) & DL_a <= 0) stop(text)
    if(length(DL_a)>1) stop(text)
    if(is.character(DL_a) & !(DL_a %in% c("hyperprior", "1/K", "1/n"))){
      stop(text)
    }
    if(DL_a == "hyperprior"){
      prior <- "DL_h"
      DL_a <- 0.5 # initial value
    }
    out <- list(prior = prior, DL_a = DL_a)

  }else if(prior == "R2D2"){

    out <- list(prior = prior, R2D2_b = R2D2_b)

  }else if(prior == "SSVS"){
    if(!(SSVS_c0>0 & SSVS_c1>0)){
      stop("'SSVS_c0' and 'SSVS_c1' must be positive numeric values.")
    }
    out <- list(prior = prior, SSVS_c0=SSVS_c0, SSVS_c1=SSVS_c1,
                semiautomatic=SSVS_semiautomatic, SSVS_s_a=SSVS_sa, SSVS_s_b=SSVS_sb)
  }else if(prior == "normal"){
    if(is.null(V_i)){
      V_i <- 10
    }
    out <- list(prior=prior, V_i=V_i)
  }else if(prior == "HMP"){
    out <- list(prior = prior, lambda_1 = HMP_lambda1, lambda_2 = HMP_lambda2)
  }
  out
}

#' @export
specify_priorL <- function(prior, DL_b = "1/n", R2D2_b = 0.5,
                             SSVS_c0 = 0.001, SSVS_c1 = 1, SSVS_sa = 0.5, SSVS_sb = 0.5,
                           HMP_lambda3 = c(0.01,0.01),
                             V_i = NULL){
  if(!(prior %in% c("DL", "HMP", "SSVS", "normal", "R2D2"))){
    stop("Argument 'prior' must be one of 'DL', 'SSVS', 'HMP' or 'normal'. \n")
  }

  if(prior == "DL"){
    text <- c("Argument 'DL_b' must be either a single positive numeric or one of 'hyperprior',
           or '1/n'. \n ")
    if(is.numeric(DL_b) & DL_b <= 0) stop(text)
    if(length(DL_b)>1) stop(text)
    if(is.character(DL_b) & !(DL_b %in% c("hyperprior", "1/n"))){
      stop(text)
    }
    if(DL_b == "hyperprior"){
      prior <- "DL_h"
      DL_b <- 0.5 # initial value
    }
    out <- list(prior = prior, DL_b = DL_b)

  }else if(prior == "R2D2"){

    out <- list(prior = prior, R2D2_b = R2D2_b)

  }else if(prior == "SSVS"){
    if(!(SSVS_c0>0 & SSVS_c1>0)){
      stop("'SSVS_c0' and 'SSVS_c1' must be positive numeric values.")
    }
    out <- list(prior = prior, SSVS_c0=SSVS_c0, SSVS_c1=SSVS_c1, SSVS_s_a=SSVS_sa, SSVS_s_b=SSVS_sb)
  }else if(prior == "normal"){
    if(is.null(V_i)){
      V_i <- 1
    }
    out <- list(prior=prior, V_i=V_i)
  }else if(prior == "HMP"){
    out <- list(prior = prior, lambda_3 = HMP_lambda3)
  }
  out
}

#' Specify hyperparameters for Dirichlet-Laplace prior on VAR coefficients
#'
#' @param a single non-negative number, indicating the concentration parameter of
#' the DL prior. Only necessary if \code{hyperhyper=FALSE}. Good properties somewhere in
#' the interval \code{[1/n,1/2]}, where \code{n} is the number of VAR coefficients.
#' Smaller values imply heavier regularization towards zero.
#' @param hyperhyper logical. \code{TRUE} imposes a discrete uniform hyperprior on
#' the concentration parameter on the interval \code{[1/n,1/2]} with 1000 support points.
#' If set to \code{FALSE}, \code{a} has to be specified.
#'
#' @return list
#' @export
specify_PHI_DLprior <- function(a=0.5, hyperhyper = FALSE){
  theta_dirichlet <-  c("a" = a)
  return(list(theta_dirichlet = theta_dirichlet,
              hyperhyper=hyperhyper))
}


#' @export
specify_L_DLprior <- function(theta_dirichlet = c("b" = 0.5),
                              hyperhyper = FALSE){
  return(list(theta_dirichlet = theta_dirichlet,
              hyperhyper=hyperhyper))
}

specify_PHI_MPprior <- function(lambda_1_gamma = c("shape" = 0.01, "rate" = 0.01),
                                lambda_2_gamma = c("shape" = 0.01, "rate" = 0.01)){
  return(list(lambda_1_gamma = lambda_1_gamma,
              lambda_2_gamma = lambda_2_gamma))
}

specify_L_MPprior <- function(lambda_4_gamma = c("shape" = 0.01, "rate" = 0.01)){
  return(list(lambda_4_gamma = lambda_4_gamma))
}

specify_PHI_SSVSprior <- function(tau_0 = 0.01, tau_1 = 100,
                                  p_i_beta = c("shape_a" = 0.5, "shape_b" = 0.5),
                                  semi_automatic = TRUE, df=NULL, rate=1/3,
                                  k_SA=5, k_AA=5){
  #df: if PHI_prior == "SSVS-mix", df are the degrees of freedom,
  # either numeric if wished to be fixed or "hyper" if wished to be treated as random variable
  #rate: if df == "hyper", rate is the rate of the exponential prior on the degrees of freedom
  return(list(tau_0 = tau_0,
              tau_1 = tau_1,
              p_i_beta = p_i_beta,
              semi_automatic = semi_automatic,
              df = df,
              rate = rate,
              k_SA = k_SA,
              k_AA = k_AA))
}

specify_L_SSVSprior <- function(kappa_0 = 0.001, kappa_1 = 1,
                                q_i_beta = c("shape_a" = 0.5,
                                             "shape_b" = 0.5),
                                df=NULL, rate=1/3,
                                k_SA=5, k_AA=5){
  return(list(kappa_0 = kappa_0,
              kappa_1 = kappa_1,
              q_i_beta = q_i_beta,
              df = df,
              rate = rate,
              k_SA = k_SA,
              k_AA = k_AA))
}

specify_PHI_nonhierarchical <- function(V_prior = 10) {
  return(V_prior)
}

specify_L_nonhierarchical <- function(V_prior = 10) {
  return(V_prior)
}

#' @export
pred_eval <- function(s, Y_obs, mod, VoI){
  # Y_obs: ex post observed data for evaluation
  # mod: model object estimated via BVAR_*
  # VoI: variables of interest for joint & marginal predictive likelihoods and MSFE

  # relevant mod settings
  SV <- mod$SV
  intercept <- mod$intercept

  # data preparation
  variables <- colnames(mod$Y)
  draws <- dim(mod$PHI)[1]
  M <- ncol(mod$Y)
  p <- mod$p

  if(nrow(Y_obs)!=s | ncol(Y_obs) !=M){
    stop("Y_obs has wrong dimensions! \n")
  }

  ## X_fore1: regressors for one-step ahead forecasts
  X_fore1 <- as.vector(t(mod$Yraw[mod$Traw:(mod$Traw-p+1),])) # Y_t:Y_(t-p+1) for one-step ahead

  if(intercept) X_fore1 <- c(X_fore1, 1)

  if(SV==TRUE) {
    sv_mu <- mod$sv_para[,1,]
    sv_phi <- mod$sv_para[,2,]
    sv_sigma <- mod$sv_para[,3,]
    sv_h_T <- mod$sv_latent[, dim(mod$sv_latent)[2],]
  }else if(SV == FALSE){
    D_draws <- mod$DRAWS$D$D_draws
  }
  Y_obs <- matrix(Y_obs, s, M)
  colnames(Y_obs) <- variables

  # storage
  PL <- rep(as.numeric(NA), draws)
  PL_joint <- matrix(as.numeric(NA), draws, s)
  PL_marginal <- array(as.numeric(NA), c(draws, s, length(VoI)), dimnames = list(NULL, NULL, VoI))

  predictions <- array(as.numeric(NA), c(draws, s, M), dimnames = list(NULL, paste0("s: ", 1:s), variables))

  h_fore <- matrix(as.numeric(NA), M, s)
  Sigma <- array(as.numeric(NA), c(M,M,s), dimnames = list(variables, variables, NULL))

  for (i in seq.int(draws)) {

    L_inv <- backsolve(mod$L[i,,], diag(M))

    if(SV==TRUE) {
      # compute k-step ahead forecasts of latent log volas
      h_fore <- sv_h_T[i, ]
      for (k in seq_len(s)) {
        mu_h <- sv_mu[i,] + sv_phi[i,]*(h_fore - sv_mu[i,])
        sigma_h <- sv_sigma[i, ]
        h_fore <- stats::rnorm(M, mean = mu_h, sd= sigma_h)

        # compute SIGMA[t+s]
        Sigma[,,k] <- t(L_inv) %*% diag(exp(h_fore)) %*% L_inv
      }

    }else if(SV==FALSE) {
      # compute SIGMA
      Sigma[,,1:s] <- t(L_inv) %*% diag(D_draws[i,]) %*% L_inv
    }

    ## one-steap ahead predictive likelihoods and predictions
    # predictive mean
    pred_mean_temp <- as.vector(X_fore1%*%mod$PHI[i,,])
    names(pred_mean_temp) <- variables

    # Predictive likelihoods
    PL[i] <- mvtnorm::dmvnorm(as.vector(Y_obs[1,]),pred_mean_temp,Sigma[,,1])
    PL_joint[i,1] <-  mvtnorm::dmvnorm(as.vector(Y_obs[1,VoI]),pred_mean_temp[VoI],Sigma[VoI,VoI,1])
    PL_marginal[i,1,] <-  stats::dnorm(as.vector(Y_obs[1, VoI]), pred_mean_temp[VoI], sqrt(diag(Sigma[VoI,VoI,1])))

    # Predictions
    predictions[i,1,] <- tryCatch(pred_mean_temp + t(chol(Sigma[,,1]))%*%stats::rnorm(M), error=function(e) MASS::mvrnorm(1, pred_mean_temp, Sigma[,,1]))

    ## s-step ahead
    # X_fore_s holds regressors for s-step ahead predictions
    # will be recursively updated
    X_fore_s <- X_fore1

    if(s>1){

      for (kk in seq_len(s-1)) {

        if(intercept){
          if(p == 1){

              X_fore_s <- c(predictions[i,kk,],1)

          }else {

              X_fore_s <- c(predictions[i,kk,], X_fore_s[1:(length(X_fore_s)-M-1)],1)

          }

        }else {

          if(p == 1){

              X_fore_s <- predictions[i,kk,]

          }else{

            X_fore_s <- c(predictions[i,kk,], X_fore_s[1:(length(X_fore_s)-M)])

          }


        }

        pred_mean_temp <- as.vector(X_fore_s%*%mod$PHI[i,,])
        names(pred_mean_temp) <- variables

        PL_joint[i, kk+1] <-  mvtnorm::dmvnorm(as.vector(Y_obs[kk+1, VoI]),pred_mean_temp[VoI],Sigma[VoI,VoI,kk+1])
        PL_marginal[i, kk+1,] <-  stats::dnorm(as.vector(Y_obs[kk+1, VoI]), pred_mean_temp[VoI], sqrt(diag(Sigma[VoI,VoI,kk+1])))

        predictions[i, kk+1,] <- tryCatch(pred_mean_temp + t(chol(Sigma[,,kk+1]))%*%stats::rnorm(M), error=function(e) MASS::mvrnorm(1, pred_mean_temp, Sigma[,,kk+1]))

      }

    }

  }

  # MSFE
  #MSFE <- (colMeans(predictions[,VoI]) - Y_obs[,VoI])^2
  MSFE <- (apply(predictions[,,VoI, drop=FALSE], 2:3, mean) - Y_obs[,VoI])^2

  return(list(predictions=predictions,
              LPL=log(mean(PL)),
              LPL_joint=log(colMeans(PL_joint)),
              LPL_marginal=log(apply(PL_marginal, 2:3, mean)),
              MSFE=MSFE,
              PL = PL,
              PL_joint = PL_joint,
              PL_marginal=PL_marginal,
              X_fore1=X_fore1))
}
