

# density of Dirichlet  ---------------------------------------------------

ddir <- function(x, alpha, log = FALSE) {
  if(!is.matrix(x)) x <- rbind(x)
  if(!is.matrix(alpha)) {
    alpha <- matrix(alpha, nrow = nrow(x), ncol = length(alpha), byrow = TRUE)
  }
  logd <- as.vector(rowSums((alpha - 1) * log(x)) + lgamma(rowSums(alpha)) - rowSums(lgamma(alpha)))
  if(log == TRUE) out <- logd else if(log == FALSE) out <- exp(logd)
  return(out)
}

ddir2 <- function(x, prep1, prep2, log = FALSE) {
  # this function is useful when a discrete uniform hyperprior is imposed on the concentration parameter of the DL prior
  # every draw of theta which depends upon the concentration parameter, has to be evaluated over a fixed grid of possible values in every iteration
  # pre calucations of prep1 and prep1 are done in the function BVAR before the posterior sampler starts
  logd <- colSums(t(prep1) * (log(x))) + prep2
  if(log == TRUE) out <- logd else if(log == FALSE) out <- exp(logd)
  return(out)
}

# Data preparation --------------------------------------------------------

data_prep <- function(Y, p, intercept){

  # Load data

  Yraw <- as.matrix(Y)
  Traw <- nrow(Yraw) # number of raw observations

  # Create design matrix X and adjust Y
  # constant (1) plus regressors (lagged values of Y)
  Ylag <- mlag(Yraw,p)          # Ylag: lagged Y matrix, will be part of X

  if(intercept){

    X <- as.matrix(cbind(Ylag, 1))

  }else X <- Ylag

  Y <- Yraw[(p+1):Traw,]        # delete first lags to match dimensions of X

  return(list(Y=Y,X=X))
}

mlag <- function(X,p){
  M<-ncol(X)
  nr <- nrow(X)-p
  nc <- p*M
  Xlag <- matrix(0, nrow = nr,ncol = nc)
  for (i in 1:(ncol(Xlag)/M)) {
    Xlag[,(i*M-M+1):(i*M)] <- X[(p-i+1):(nr+p-i),]
  }
  Xlag
}



# MP Prior functions ------------------------------------------------------

MP_sigma_sq <- function(Y,p){ # X,

  T <- nrow(Y)
  M <- ncol(Y)
  sigma_sq <- rep(as.numeric(NA), M)

  for(i in 1:M){

    Y_i <- as.matrix(Y[,i])
    X_i <- mlag(Y_i, p)
    Y_i <- Y_i[-(1:p),]
    phi_i <- tryCatch(solve(crossprod(X_i))%*%crossprod(X_i, Y_i), error = function(e) solve(qr(X_i, LAPACK = TRUE), Y_i)) # OLS estimates of i-th equation
    sigma_sq[i] <- (crossprod(Y_i-X_i%*%phi_i))/(T-p)

  }
  return(sigma_sq)
}

MP_V_prior_prep <- function(sigma_sq ,K, M, intercept){
  # sigma_sq: vector of variances of univariate AR(p) models of the individual variables estimated via "MP_sigma_sq"
  # K: number of regressors per equation
  # M: number of variables

  V_i <- matrix(0,nrow = K, ncol = M)

  for(i in 1:M){ #for each i-th equation
    for(j in 1:K){ #for each j-th RHS variable
      if(j==K & intercept){
        V_i[j,i] <- sigma_sq[i] #variance on constant
      }else if((j-i)%%M==0){
        r <- ceiling((j)/M)
        V_i[j,i] <- 1/(r^2) #variance on own lags
      }else{
        r <- ceiling((j)/M)
        if(j%%M==0){
          V_i[j,i] <- (sigma_sq[i])/((r^2)*sigma_sq[(M)])
        }else if(j%%M==1){
          V_i[j,i] <- (sigma_sq[i])/((r^2)*sigma_sq[1])
        }else {
          V_i[j,i] <- (sigma_sq[i])/((r^2)*sigma_sq[(j%%M)])
        }
      }
    }
  }
  return(V_i_prep=as.vector(V_i))
}

get_MP_V_prior <- function(lambda1=0.04, lambda2=0.0016, lambda3=10^3, sigma_sq ,K, M, intercept){
  # sigma_sq: vector of variances of univariate AR(p) models of the individual variables estimated via "MP_sigma_sq"
  # K: number of regressors per equation
  # M: number of variables

  V_i <- matrix(0,nrow = K, ncol = M)

  for(i in 1:M){ #for each i-th equation
    for(j in 1:K){ #for each j-th RHS variable
      if(j==K & intercept){
        V_i[j,i] <- lambda3 * sigma_sq[i] #variance on constant
      }else if((j-i)%%M==0){
        r <- ceiling((j)/M)
        V_i[j,i] <- lambda1/(r^2) #variance on own lags
      }else{
        r <- ceiling((j)/M)
        if(j%%M==0){
          V_i[j,i] <- (lambda2 * sigma_sq[i])/((r^2)*sigma_sq[(M)])
        }else if(j%%M==1){
          V_i[j,i] <- (lambda2 * sigma_sq[i])/((r^2)*sigma_sq[1])
        }else {
          V_i[j,i] <- (lambda2 * sigma_sq[i])/((r^2)*sigma_sq[(j%%M)])
        }
      }
    }
  }
  return(V_i=as.vector(V_i))
}


# Companion ---------------------------------------------------------------

get_companion <- function(PHI, p, intercept=TRUE){

  M <- ncol(PHI)
  K <- nrow(PHI)

  companion <- matrix(0, K, K)
  companion[,1:M] <- PHI
  if(p>1){
    companion[1:((p-1)*M), (M+1):(M+(p-1)*M)] <- diag((p-1)*M)
  }
  if(intercept){
    companion[K, K] <- 1
  }

  #cc <- matrix(0, M*p + intercept,M)
  #cc[1:M,1:M] <- diag(M)

  return(companion
    #list(companion = companion,
          #    cc = cc)
    )

}


# Density of multivariate Normal (Cholesky) -------------------------------

mydmvnorm <- function (x, mean , cholsigma, log = FALSE) {
  p <- ncol(cholsigma)
  #x <- matrix(x, ncol = p)
  mean <- as.vector(mean)
  x <- as.vector(x)

  tmp <- backsolve(cholsigma, x - mean, transpose = TRUE)
  rss <- sum(tmp^2)
  logretval <- -sum(log(diag(cholsigma))) - 0.5 * p * log(2 *pi) - 0.5 * rss

  names(logretval) <- rownames(x)
  if (log)
    logretval
  else exp(logretval)
}

