# specify Priors ----------------------------------------------------------

#' Specify hyperparameters for Dirichlet-Laplace prior on VAR coefficients
#'
#' @param a numeric constant. Concentration parameter, that dictates the shrinkage. Only necessary if ´hyperhyper = FALSE'.
#' @param hyperhyper bool. TRUE: Discrete uniform hyperprior on a. FALSE: `a' has to be defined.
#'
#' @return list
#' @export
specify_PHI_DLprior <- function(a=0.5, hyperhyper = FALSE){
  theta_dirichlet <-  c("a" = a)
  return(list(theta_dirichlet = theta_dirichlet,
              hyperhyper=hyperhyper))
}

#' @export
specify_PHI_DLprior_new <- function(type = "own/cross", a_1=0.5, a_2=1/10,
                                    hyperhyper = FALSE){
  #type: "global": one shrinkage parameter for all coefficients;
  #"own/cross": seperate shrinkage parameters for own- and crosslags;
  #"lagwise-own/cross": seperate shrinkage parameters for own- and crosslags for every lag
  #a_1: if type=="global": global shrinkage parameter;
  #if type=="own/cross: shrinkage parameter for own lags;
  #if type=="lagwise-own/cross": shrinkage parameters for own lags, vector of length p;
  #a_2: if type=="global": redundant (NULL);
  #if type=="own/cross: shrinkage parameter for cross lags;
  #if type=="lagwise-own/cross": shrinkage parameters for cross lags, vector of length p;
  #hyperhyper: discrete uniform hyperprior on the shrinkage hyperparameters: no need to specify a_1 and a_2;

  return(list(type=type,
              a_1=a_1,
              a_2=a_2,
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

pred_eval <- function(h, Y_obs, mod, PHI_draws, L_draws, VoI, SV = TRUE){
  # Y_obs: ex post observed data for evaluation
  # mod: model object estimated via BVAR_*
  # PHI_draws: either posterior draws or savs draws
  # L_draws: either posterior or savs draws
  # VoI: variables of interest for joint & marginal predictive likelihoods and MSFE

  # data preparation
  intercept <- mod$intercept
  variables <- colnames(mod$Y)
  draws <- dim(PHI_draws)[1]
  M <- ncol(mod$Yraw)#############
  p <- mod$p

  X_fore1 <- mod$X_fore

  if(SV==TRUE) {
    mus <- mod$Posterior_draws$D$SVpara_draws[,1,]
    phis <- mod$Posterior_draws$D$SVpara_draws[,2,]
    sigs <- mod$Posterior_draws$D$SVpara_draws[,3,]
    h_ts <- mod$Posterior_draws$D$SVlatent_draws[, dim(mod$Posterior_draws$D$SVlatent_draws)[2],]
  }else if(SV == FALSE){
    D_draws <- mod$DRAWS$D$D_draws
  }
  if(mod$standardize == TRUE){
    Y_obs <- t((t(Y_obs) - mod$mu_Y)/mod$sd_Y)
  }
  Y_obs <- as.matrix(Y_obs, h, M)
  colnames(Y_obs) <- variables

  # storage
  PL <- rep(as.numeric(NA), draws)
  PL_joint <- matrix(as.numeric(NA), draws, h)
  PL_marginal <- array(as.numeric(NA), c(draws, h, length(VoI)), dimnames = list(NULL, NULL, VoI))

  predictions <- array(as.numeric(NA), c(draws, h, M), dimnames = list(NULL, NULL, variables))

  h_fore <- matrix(as.numeric(NA), M, h)
  Sigma <- array(as.numeric(NA), c(M,M,h), dimnames = list(variables, variables, NULL))

  for (i in seq.int(draws)) {

    L_inv <- backsolve(L_draws[i,,], diag(M))

    if(SV==TRUE) {
      # compute k-step ahead forecasts of latent log volas
      h_fore <- h_ts[i, ]
      for (k in seq_len(h)) {
        mu_h <- mus[i,] + phis[i,]*(h_fore - mus[i,])
        sigma_h <- sigs[i, ]
        h_fore <- stats::rnorm(M, mean = mu_h, sd= sigma_h)

        # compute SIGMA[t+h]
        Sigma[,,k] <- t(L_inv) %*% diag(exp(h_fore)) %*% L_inv
      }

    }else if(SV==FALSE) {
      # compute SIGMA
      Sigma[,,1:h] <- t(L_inv) %*% diag(D_draws[i,]) %*% L_inv
    }

    ## one-steap ahead predictive likelihoods and predictions
    # predictive mean
    pred_mean_temp <- as.vector(X_fore1%*%PHI_draws[i,,])
    names(pred_mean_temp) <- variables

    # Predictive likelihoods
    PL[i] <- mvtnorm::dmvnorm(as.vector(Y_obs[1,]),pred_mean_temp,Sigma[,,1])
    PL_joint[i,1] <-  mvtnorm::dmvnorm(as.vector(Y_obs[1,VoI]),pred_mean_temp[VoI],Sigma[VoI,VoI,1])
    PL_marginal[i,1,] <-  stats::dnorm(as.vector(Y_obs[1, VoI]), pred_mean_temp[VoI], sqrt(diag(Sigma[VoI,VoI,1])))

    # Predictions
    predictions[i,1,] <- tryCatch(pred_mean_temp + t(chol(Sigma[,,1]))%*%stats::rnorm(M), error=function(e) MASS::mvrnorm(1, pred_mean_temp, Sigma[,,1]))

    ## h-step ahead
    X_fore_h <- X_fore1
    if(h>1){

      for (kk in seq_len(h-1)) {

        if(intercept){
          if(p == 1){

            X_fore_h <- c(predictions[i,kk,],1)

          }else X_fore_h <- c(predictions[i,kk,], X_fore_h[1:(length(X_fore_h)-M-1)],1)

        }else {
          if(p == 1){

            X_fore_h <- predictions[i,kk,]

          }else X_fore_h <- c(predictions[i,kk,], X_fore_h[1:(length(X_fore_h)-M)])

        }

        pred_mean_temp <- as.vector(X_fore_h%*%PHI_draws[i,,])
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
              MSFE=MSFE))
}
