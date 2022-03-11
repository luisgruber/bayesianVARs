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
                             V_i = 10,...){
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

    out <- list(prior = prior, R2D2_b = R2D2_b,...)

  }else if(prior == "SSVS"){
    if(!(all(SSVS_c0>0) & all(SSVS_c1>0))){
      stop("'SSVS_c0' and 'SSVS_c1' must be positive numeric values. \n")
    }
    out <- list(prior = prior, SSVS_c0=SSVS_c0, SSVS_c1=SSVS_c1,
                semiautomatic=SSVS_semiautomatic, SSVS_s_a=SSVS_sa, SSVS_s_b=SSVS_sb)
  }else if(prior == "normal"){
    if(!(all(V_i>0))){
      stop("'V_i' must be positive. \n")
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
                             V_i = 10,
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
    if(!(all(V_i>0))){
      stop("'V_i' must be positive. \n")
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
predict.bvar <- function(mod, nsteps, LPL = FALSE, Y_obs = NA, LPL_VoI = NA,...){

  # relevant mod settings
  SV <- mod$SV
  intercept <- mod$intercept

  # data preparation
  variables <- colnames(mod$Y)
  draws <- dim(mod$PHI)[1]
  p <- mod$p
  M <- ncol(mod$Y)
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
    Y_obs <- matrix(Y_obs, nsteps, M)
    colnames(Y_obs) <- variables
  }
  if(SV==TRUE) {
    # extract sv parameters
    sv_mu <- mod$sv_para[,1,]
    sv_phi <- mod$sv_para[,2,]
    sv_sigma <- mod$sv_para[,3,]
    # extract current state of log-vola
    sv_h_T <- mod$sv_latent[, dim(mod$sv_latent)[2],]
  }else if(SV == FALSE){
    D_sqrt_draws <- exp(mod$sv_latent[, dim(mod$sv_latent)[2],]/2)
  }

  # storage
  predictions <- array(as.numeric(NA), c(draws, nsteps, M),
                       dimnames = list(NULL, paste0("t+", 1:nsteps), variables))
  if(LPL){
    LPL_draws <- matrix(as.numeric(NA), draws, nsteps)
    colnames(LPL_draws) <- paste0("t+", 1:nsteps)
    PL_univariate_draws <- array(as.numeric(NA), c(draws, nsteps, M),
                                 dimnames = list(NULL, paste0("t+", 1:nsteps), variables))
    if(LPL_subset){
      LPL_sub_draws <- matrix(as.numeric(NA), draws, nsteps)
      colnames(LPL_sub_draws) <- paste0("t+", 1:nsteps)
    }
  }

  ## X_fore1: predictors for one-step ahead forecasts
  X_fore1 <- as.vector(t(mod$Yraw[mod$Traw:(mod$Traw-p+1),]))

  if(intercept) X_fore1 <- c(X_fore1, 1)

  for (i in seq.int(draws)) {

    L_inv <- backsolve(mod$L[i,,], diag(M))

      X_fore_k <- X_fore1

      if(!SV){
        # compute SIGMA
        #Sigma_fore <- crossprod(L_inv, diag(D_draws[i,])) %*% L_inv #???
        #t(L_inv) %*% diag(D_draws[i,]) %*% L_inv #???
        Sigma_chol_fore <- diag(D_sqrt_draws[i,]) %*% L_inv
        Sigma_fore <- crossprod(Sigma_chol_fore)
      }else if(SV){
        # initialize latent log-vola at current state
        h_fore <- sv_h_T[i, ]
      }

      for(k in seq.int(nsteps)){

        mean_fore <- as.vector(X_fore_k%*%mod$PHI[i,,])

        # compute prediction of variance-covariance matrix
        if(SV){
          # compute k-step ahead forecast of latent log-volas
          h_fore <- sv_mu[i,] + sv_phi[i,]*(h_fore - sv_mu[i,]) +
            sv_sigma[i,]*stats::rnorm(M, mean = 0, sd = 1)

          # compute SIGMA[t+k]
          #Sigma_fore <- crossprod(L_inv, diag(exp(h_fore))) %*% L_inv #???
          #t(L_inv) %*% diag(exp(h_fore)) %*% L_inv
          #???
          Sigma_chol_fore <- diag(exp(h_fore/2)) %*% L_inv
          Sigma_fore <- crossprod(Sigma_chol_fore)
        }

        predictions[i,k,] <- tryCatch(
          mean_fore + stats::rnorm(M) %*% Sigma_chol_fore, #??? mean_fore + t(chol(Sigma_fore)) %*% stats::rnorm(M)
          error = function(e) MASS::mvrnorm(1, mean_fore, Sigma_fore)
          )

        if(LPL){
          LPL_draws[i,k] <- mydmvnorm(Y_obs[k,], mean_fore,Sigma_chol_fore, log = TRUE) #??? mvtnorm::dmvnorm(as.vector(Y_obs[k,]),mean_fore,Sigma_fore, log = TRUE)
          PL_univariate_draws[i, k,] <-  stats::dnorm(as.vector(Y_obs[k, ]), mean_fore, sqrt(diag(Sigma_fore)))
          if(LPL_subset){
            LPL_sub_draws[i, k] <-  mydmvnorm(Y_obs[k, VoI], mean_fore[VoI], chol(Sigma_fore[VoI,VoI, drop = FALSE]), log = TRUE)#??? mvtnorm::dmvnorm(as.vector(Y_obs[k, VoI]),mean_fore[VoI],Sigma_fore[VoI,VoI, drop = FALSE], log = TRUE)
          }
        }

        if(k<nsteps){
          if(p==1){
            X_fore_k <- predictions[i,k,]
          }else{
            X_fore_k <- c(predictions[i,k,], X_fore_k[1:((p-1)*M)])
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
    names(LPL) <- paste0("t+", 1:nsteps)
    out$LPL <- LPL
    out$LPL_draws <- LPL_draws

    out$LPL_univariate <- log(apply(PL_univariate_draws, 2:3, mean))
    rownames(out$LPL_univariate) <- paste0("t+", 1:nsteps)
    out$PL_univariate_draws <- PL_univariate_draws

    if(LPL_subset){
      numericalnormalizer2 <- apply(LPL_sub_draws,2,max) - 700
      LPL_VoI <- log(colMeans(exp( t(t(LPL_sub_draws) - numericalnormalizer2)))) +
        numericalnormalizer2
      names(LPL_VoI) <- paste0("t+", 1:nsteps)
      out$LPL_VoI <- LPL_VoI
      out$LPL_VoI_draws <- LPL_sub_draws
      out$VoI <- variables[VoI]
    }
  }
  class(out) <- "bvar_predict"
  return(out)
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
  F <- get_companion(mod$PHI[1,,], p, intercept)
  for (i in seq.int(draws)) {

    # Companion form
    # F: matrix of coefficients
    # z: endogenous variables
    F[,1:M] <- mod$PHI[i,,]
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

#' @export
summary.bvar_predict <- function(object, ...){
  out <- list()
  if(!is.null(object$LPL)){
    out$LPL <- object$LPL
    out$LPL_univariate <- object$LPL_univariate
  }
  if(!is.null(object$LPL_VoI)){
    out$LPL_VoI <- object$LPL_VoI
    out$VoI <- object$VoI
  }
  out$prediction_quantiles <- apply(object$predictions, 2:3, quantile, c(.05,.5,.95))
  class(out) <- "summary.bvar_predict"
  out
}

#' @export
print.summary.bvar_predict <- function(x, ...){
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

  cat("\nPrediction quantiles:\n")
  print(x$prediction_quantiles, digits = digits)

  invisible(x)
}
