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
                             V_i = 10, global_local_grouping="global",...){
  if(!(prior %in% c("DL", "HMP", "SSVS", "normal", "R2D2", "SL"))){
    stop("Argument 'prior' must be one of 'DL', 'SSVS', 'HMP' or 'normal'. \n")
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

    if(is.character(global_local_grouping)){
      if(!(global_local_grouping %in% c("global", "equation-wise", "covariate-wise", "fol", "olcl-lagwise"))){
        stop("Argument 'global_local_grouping' must be one of 'global',
           'equation-wise', 'covariate-wise', 'olcl-lagwise' or 'fol'. \n")
      }
    }
    out <- list(prior = prior, DL_a = DL_a, global_local_grouping = global_local_grouping)

  }else if(prior == "R2D2"){
    if(is.character(global_local_grouping)){
      if(!(global_local_grouping %in% c("global", "equation-wise", "covariate-wise", "fol", "olcl-lagwise"))){
        stop("Argument 'global_local_grouping' must be one of 'global',
           'equation-wise', 'covariate-wise', 'olcl-lagwise' or 'fol'. \n")
      }
    }

    out <- list(prior = prior, R2D2_b = R2D2_b, global_local_grouping = global_local_grouping)

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

#' @export
summary.bvar <- function(object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),...){
  PHImedian <- apply(object$PHI, 2:3, median)
  PHIquantiles <- apply(object$PHI, 2:3, quantile, quantiles)
  PHIiqr <- apply(object$PHI, 2:3, IQR)
  out <- list(PHImedian = PHImedian,
             PHIquantiles = PHIquantiles,
             PHIiqr = PHIiqr)
  class(out) <- "summary.bvar"
  out
}

#' @export
print.summary.bvar <- function(x, ...){
  digits <- max(3, getOption("digits") - 3)
    cat("\nPosterior median of PHI:\n")
    print(x$PHImedian, digits = digits)
    cat("\nPosterior interquartile range of PHI:\n")
    print(x$PHIiqr, digits = digits)

  invisible(x)
}
