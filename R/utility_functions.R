# specify Priors ----------------------------------------------------------

#' Specify prior on PHI
#'
#' Specify
#'
#' @param prior
#'
#' @param DL_a
#' @param R2D2_b
#' @param SSVS_c0
#' @param SSVS_c1
#' @param SSVS_semiautomatic
#' @param SSVS_sa
#' @param SSVS_sb
#' @param HMP_lambda1
#' @param HMP_lambda2
#' @param V_i
#' @param ...
#'
#' @export
specify_priorPHI <- function(prior, DL_a = "1/K", R2D2_b = 0.5,
                             SSVS_c0 = 0.01, SSVS_c1 = 100, SSVS_semiautomatic = TRUE, SSVS_sa = 0.5, SSVS_sb = 0.5,
                             HMP_lambda1 = c(0.01,0.01), HMP_lambda2 = c(0.01,0.01),
                             V_i = NULL,...){
  if(!(prior %in% c("DL", "HMP", "SSVS", "normal", "R2D2", "SL"))){
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
  }else if(prior == "SL"){
    out <- list(prior = prior, ...)
  }
  out
}

#' Specify prior on L
#'
#' Specify
#'
#' @param prior
#'
#' @param DL_b
#' @param R2D2_b
#' @param SSVS_c0
#' @param SSVS_c1
#' @param SSVS_sa
#' @param SSVS_sb
#' @param HMP_lambda3
#' @param V_i
#' @param ...
#'
#' @export
specify_priorL <- function(prior, DL_b = "1/n", R2D2_b = 0.5,
                             SSVS_c0 = 0.001, SSVS_c1 = 1, SSVS_sa = 0.5, SSVS_sb = 0.5,
                           HMP_lambda3 = c(0.01,0.01),
                             V_i = NULL,
                           ...){
  if(!(prior %in% c("DL", "HMP", "SSVS", "normal", "R2D2", "SL"))){
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
  }else if(prior == "SL"){
    out <- list(prior = "SL", ...)
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
pred_eval <- function(s, Y_obs, mod, VoI, new = TRUE){
  # Y_obs: ex post observed data for evaluation
  # mod: model object estimated via BVAR_*
  # VoI: variables of interest for joint & marginal predictive likelihoods

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

  ## X_fore1: predictors for one-step ahead forecasts
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
  LPL_draws <- rep(as.numeric(NA), draws)
  LPL_joint_draws <- matrix(as.numeric(NA), draws, s)
  PL_marginal_draws <- array(as.numeric(NA), c(draws, s, M), dimnames = list(NULL, NULL, variables))

  predictions <- array(as.numeric(NA), c(draws, s, M), dimnames = list(NULL, paste0("s: ", 1:s), variables))

  h_fore <- matrix(as.numeric(NA), M, s)
  Sigma <- array(as.numeric(NA), c(M,M,s), dimnames = list(variables, variables, NULL))
 # Sigma_fore <- matrix(NA, M,M)


  for (i in seq.int(draws)) {

    L_inv <- backsolve(mod$L[i,,], diag(M))

      X_fore_k <- X_fore1

      if(!SV){
        # compute SIGMA
        Sigma_fore <- t(L_inv) %*% diag(D_draws[i,]) %*% L_inv
        rownames(Sigma_fore) <- colnames(Sigma_fore) <- variables
      }else if(SV){
        # initialize latent vola
        h_fore <- sv_h_T[i, ]
      }

      for(k in seq.int(s)){

        pred_mean_temp <- as.vector(X_fore_k%*%mod$PHI[i,,])
        names(pred_mean_temp) <- variables

        # compute prediction of variance-covariance matrix
        if(SV){
          # compute k-step ahead forecasts of latent log volas
          mu_h <- sv_mu[i,] + sv_phi[i,]*(h_fore - sv_mu[i,])
          sigma_h <- sv_sigma[i, ]
          h_fore <- stats::rnorm(M, mean = mu_h, sd= sigma_h)

          # compute SIGMA[t+s]
          Sigma_fore <- t(L_inv) %*% diag(exp(h_fore)) %*% L_inv
          rownames(Sigma_fore) <- colnames(Sigma_fore) <- variables
        }

        predictions[i,k,] <- tryCatch(pred_mean_temp + t(chol(Sigma_fore))%*%stats::rnorm(M),
                                      error=function(e) MASS::mvrnorm(1, pred_mean_temp, Sigma_fore))

        if(k==1){
          LPL_draws[i] <- mvtnorm::dmvnorm(as.vector(Y_obs[1,]),pred_mean_temp,Sigma_fore, log = TRUE)
        }

        LPL_joint_draws[i, k] <-  mvtnorm::dmvnorm(as.vector(Y_obs[k, VoI]),pred_mean_temp[VoI],Sigma_fore[VoI,VoI], log = TRUE)
        PL_marginal_draws[i, k,] <-  stats::dnorm(as.vector(Y_obs[k, ]), pred_mean_temp, sqrt(diag(Sigma_fore)))

        if(k<s){
          X_fore_k <- c(predictions[i,k,], X_fore_k[1:((p-1)*M)])
          if(intercept){
            X_fore_k <- c(X_fore_k,1)
          }
        }

      }# end 1:k
  }# end 1:draws

  numericalnormalizer <- max(LPL_draws) - 700
  LPL <- log(mean(exp(LPL_draws - numericalnormalizer))) + numericalnormalizer
  LPL0 <- log(mean(exp(LPL_draws)))

  numericalnormalizer2 <- apply(LPL_joint_draws,2,max) - 700
  LPL_joint <- log(colMeans(exp( t(t(LPL_joint_draws) - numericalnormalizer2)))) + numericalnormalizer2
  LPL_joint0 <- log(colMeans(exp(LPL_joint_draws)))

  return(list(predictions=predictions,
              LPL=LPL,
              LPL0=LPL0,
              LPL_joint=LPL_joint,
              LPL_joint0=LPL_joint0,
              LPL_marginal=log(apply(PL_marginal_draws, 2:3, mean)),
              LPL_draws = LPL_draws,
              LPL_joint_draws = LPL_joint_draws,
              PL_marginal_draws=PL_marginal_draws,
              X_fore1=X_fore1))
}

#' @export
predict_bvar <- function(s,mod, LPL = FALSE, Y_obs = NA, LPL_VoI = NA){
  # Y_obs: ex post observed data for evaluation
  # mod: model object estimated via BVAR_*
  # VoI: variables of interest for joint & marginal predictive likelihoods

  # relevant mod settings
  SV <- mod$SV
  intercept <- mod$intercept

  # data preparation
  variables <- colnames(mod$Y)
  draws <- dim(mod$PHI)[1]
  M <- ncol(mod$Y)
  K <- ncol(mod$X)
  p <- mod$p

  if(nrow(Y_obs)!=s | ncol(Y_obs) !=M){
    stop("Y_obs has wrong dimensions! \n")
  }


  if(SV==TRUE) {
    sv_mu <- mod$sv_para[,1,]
    sv_phi <- mod$sv_para[,2,]
    sv_sigma <- mod$sv_para[,3,]
    sv_h_T <- mod$sv_latent[, dim(mod$sv_latent)[2],]
  }else if(SV == FALSE){
    D_draws <- mod$DRAWS$D$D_draws
  }

  # storage
  predictions <- array(as.numeric(NA), c(draws, s, M), dimnames = list(NULL, paste0("s: ", 1:s), variables))

  if(LPL){
    Y_obs <- matrix(Y_obs, s, M)
    colnames(Y_obs) <- variables

    LPL_draws <- matrix(as.numeric(NA), draws, s)
    PL_marginal_draws <- array(as.numeric(NA), c(draws, s, M), dimnames = list(NULL, NULL, variables))
    if(!any(is.na(LPL_VoI))){
      LPL_joint_draws <- matrix(as.numeric(NA), draws, s)
    }
  }

  # initialization and placeholders
  h_fore <- matrix(as.numeric(NA), M, s)
  Sigma_padding <- matrix(0, M*p + intercept,M)
  Sigma_padding[1:M,1:M] <- diag(M)
  z_fore0 <- as.matrix(mod$datamat[nrow(mod$datamat), 1:(p*M)])
  if(intercept){
    z_fore0 <- cbind(z_fore0,1)
  }

  for (i in seq.int(draws)) {

    # Companion form
    # F: matrix of coefficients
    # z: endogenous variables
    F <- get_companion(mod$PHI[i,,], p, intercept)
    z_fore <- z_fore0
    SIGMA_k <- matrix(0,K,K)

    L_inv <- backsolve(mod$L[i,,], diag(M))

    if(!SV){
      # compute SIGMA
      Sigma_fore <- t(L_inv) %*% diag(D_draws[i,]) %*% L_inv
    }else if(SV){
      # initialize latent vola
      h_fore <- sv_h_T[i, ]
    }

    for(k in seq_len(s)){

      # recursively update the predictive mean
      z_fore <- z_fore %*% F
      pred_mean <- z_fore[1:M]
      names(pred_mean) <- variables

      # compute prediction of variance-covariance matrix
      if(SV){
        # compute k-step ahead forecasts of latent log volas
        mu_h <- sv_mu[i,] + sv_phi[i,]*(h_fore - sv_mu[i,])
        sigma_h <- sv_sigma[i, ]
        h_fore <- stats::rnorm(M, mean = mu_h, sd = sigma_h)

        # compute SIGMA[t+s]
        Sigma_fore <- t(L_inv) %*% diag(exp(h_fore)) %*% L_inv
      }
      SIGMA_k <- crossprod(F,SIGMA_k) %*% F + Sigma_padding %*% tcrossprod(Sigma_fore, Sigma_padding)
      SIGMA_k_crop <- SIGMA_k[1:M,1:M]
      rownames(SIGMA_k_crop) <- colnames(SIGMA_k_crop) <- variables

      predictions[i,k,] <- tryCatch(pred_mean + t(chol(SIGMA_k_crop ))%*%stats::rnorm(M),
               error=function(e) MASS::mvrnorm(1, pred_mean, SIGMA_k_crop))

      if(LPL){
        LPL_draws[i,k] <- mvtnorm::dmvnorm(as.vector(Y_obs[k,]),pred_mean,SIGMA_k_crop, log = TRUE)
        PL_marginal_draws[i, k,] <-  stats::dnorm(as.vector(Y_obs[k, ]), pred_mean, sqrt(diag(SIGMA_k_crop)))
        if(!any(is.na(LPL_VoI))){
          LPL_joint_draws[i, k] <-  mvtnorm::dmvnorm(as.vector(Y_obs[k, LPL_VoI]),pred_mean[LPL_VoI],SIGMA_k_crop[LPL_VoI,LPL_VoI], log = TRUE)
        }
      }

    }# end 1:s step ahead

  }# end 1:draws

  out <- list(predictions = predictions)

  if(LPL){
    numericalnormalizer <- apply(LPL_draws, 2, function(x) max(x) - 700)
    LPL <- apply(t(t(LPL_draws) - numericalnormalizer), 2, function(x) log(mean(exp(x)))) +
      numericalnormalizer
    out$LPL <- LPL
    out$LPL_draws <- LPL_draws
    out$LPL_marginal = log(apply(PL_marginal_draws, 2:3, mean))
    out$PL_marginal_draws <- PL_marginal_draws
  }

  if(!any(is.na(LPL_VoI))){
    numericalnormalizer2 <- apply(LPL_joint_draws, 2, function(x) max(x) - 700)
    LPL_joint <- apply(t(t(LPL_draws) - numericalnormalizer), 2, function(x) log(mean(exp(x)))) +
      numericalnormalizer
    out$LPL_joint <- LPL_joint
    out$LPL_joint_draws <- LPL_joint_draws
  }

  return(out)
}
