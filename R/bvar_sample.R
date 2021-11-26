#' Title
#'
#' @param Yraw data
#' @param p lag-length
#' @param intercept logical
#' @param persistence single value in the interval \code{[0,1]} indicating persistence
#' of first own-lag covariates. I.e. the prior mean for the coefficients representing
#' first own-lags.
#' @param PHI_prior character. Prior on VAR coefficients: \code{"MP"} for the Hierarchical Minnesota Prior,
#' \code{"SSVS"} for Stochastic Search Variable Selection or \code{"DL"} for the
#' Dirichlet-Laplace prior.
#' @param PHI_hyper list. Specify prior settings on PHI. Use \code{specify_PHI_*prior},
#' where \code{*} is one of {MP, SSVS} or \code{DL}.
#' @param L_prior character. Prior on impact matrix L (LDL' decomposition of
#' variance-covariance matrix): \code{"MP"} for the Hierarchical Minnesota Prior,
#' \code{"SSVS"} for Stochastic Search Variable Selection or \code{"DL"} for the
#' Dirichlet-Laplace prior.
#' @param L_hyper list. Specify prior settings on L. Use \code{specify_L_*prior},
#' where \code{*} is one of {MP, SSVS} or \code{DL}.
#' @param draws Integer scalar. Number of posterior draws to be saved.
#' @param burnin Integer scalar. Number of burnin MCMC iterations, where draws
#' will not be saved.
#' @param SV logical. If \code{TRUE} stochastic volatility specification for orthogonalized
#' errors, i.e. time varying variance covariance matrix \eqn{\Sigma_t=LD_t L^\prime} is assumed.
#' If \code{FALSE} constant variance-covariance matrix is assumed with
#' inverse Gamma prior with \code{scale=0.01} and \code{rate=0.01}
#' on diagonal elements of D.
#' @param SV_hyper list. Use \code{specify_priors}; imported helper function from
#' package \code{stochvol}.
#' @param standardize logical. If \code{TRUE}, then data is standardized to have
#' zero mean and unit variance. That is the default, because priors are not invariant
#' to different scales of the data.
#'
#' @return Something
#' @export
bvar <- function(Yraw, p, intercept = FALSE, persistence = 0, PHI_prior = "DL", PHI_hyper = NULL,
                 L_prior = "DL", L_hyper = NULL, draws, burnin, SV = TRUE,
                 SV_hyper = NULL, standardize = TRUE, PHI_in = NULL, L_in = NULL) {

  #standardize: scales data (after Y (matrix of responses) and X (design matrix) are created) to have zero mean and variance 1.

  tot <- draws + burnin
  pb <- txtProgressBar(min = 1, max = tot, initial = 0, style = 3)

  # Data preparation --------------------------------------------------------

  data <- data_prep(Yraw, p, intercept)

  Y <- data$Y
  X <- data$X
  # before scaling, store mu and sd of Y
  # is needed for predictive evaluation to scale Y_observed
  mu_Y <- colMeans(Y)
  sd_Y <- apply(Y, 2 ,sd)

  if(standardize==TRUE){

    Y <- scale(Y)
    if(intercept){
      X[,-ncol(X)] <- scale(X[,-ncol(X)])
    }else X <- scale(X)

  }

  XX <- crossprod(X)
  M <- ncol(Y)
  K <- ncol(X)
  T <- nrow(Y)
  n <- K*M
  n_l <- (M*(M-1))/2 # number of elements in upper tri of Sigma_t
  variables <- colnames(Yraw)

  # Storage -----------------------------------------------------------------

  DRAWS <- vector("list", length = 4L)
  names(DRAWS) <- c("PHI", "L", "D", "Ytilde_draws")
  DRAWS$Ytilde_draws <- array(as.numeric(NA), c(draws, T, M))
  V_i_draws <- matrix(as.numeric(NA), draws, n)

  if(PHI_prior == "MP") {

    DRAWS$PHI <- list(PHI_draws = array(as.numeric(NA), c(draws, K, M)),
                      lambda_1_draws = rep(as.numeric(NA), draws),
                      lambda_2_draws = rep(as.numeric(NA), draws))

  }else if(PHI_prior == "DL") {

    DRAWS$PHI <- list(PHI_draws = array(as.numeric(NA), c(draws, K, M)),
                      psi_draws = matrix(as.numeric(NA), draws, n),
                      zeta_draws = rep(as.numeric(NA), draws),
                      theta_draws = matrix(as.numeric(NA), draws, n),
                      a_draws = rep(as.numeric(NA), draws))

  }else if(PHI_prior == "SSVS") {

    DRAWS$PHI <- list(PHI_draws = array(as.numeric(NA), c(draws, K, M)),
                      gamma_draws = matrix(0, draws, n),
                      p_i_draws = matrix(as.numeric(NA), draws, n))

  }else if(PHI_prior == "non_hierarchical") {

    DRAWS$PHI <- list(PHI_draws = array(as.numeric(NA), c(draws, K, M)))

  }

  if(L_prior == "non_hierarchical") {

    DRAWS$L <- list(L_draws = array(as.numeric(NA), c(draws, M, M)))

  }else if(L_prior == "DL") {

    DRAWS$L <- list(L_draws = array(as.numeric(NA), c(draws, M, M)),
                    psi_l_draws = matrix(as.numeric(NA), draws, n_l),
                    zeta_l_draws = rep(as.numeric(NA), draws),
                    theta_l_draws = matrix(as.numeric(NA), draws, n_l),
                    b_draws = rep(as.numeric(NA), draws))

  }else if(L_prior == "SSVS") {

    DRAWS$L <- list(L_draws = array(as.numeric(NA), c(draws, M, M)),
                    omega_draws = matrix(0, draws, n_l),
                    q_i_draws = matrix(as.numeric(NA), draws, n_l))

  }else if(L_prior == "MP") {

    DRAWS$L <- list(L_draws = array(as.numeric(NA), c(draws, M, M)),
                    lambda_4_draws = rep(as.numeric(NA), draws))

  }

  if(SV == TRUE) {

    DRAWS$D <- list(SVpara_draws = array(as.numeric(NA), c(draws,4,M)),
                    SVlatent_draws = array(as.numeric(NA),c(draws,T,M)))

  }else if (SV == FALSE) {

    DRAWS$D <- list(D_draws = matrix(as.numeric(NA), draws, M))

  }

  # Prior hyperparameters ---------------------------------------------------

  # Prior mean of PHI
  if(persistence == 0) {

    PHI0 <- matrix(0, K, M)

  }else {

    PHI0 <- matrix(0, K, M)
    PHI0[1:M, 1:M] <- diag(M)*persistence

  }
  phi_prior <- as.vector(PHI0)

  if(PHI_prior == "MP") {

    # ind_matrix: indicates whether coefficient is
    ## intercept (constant): 3,
    ## own lag: 1,
    ## cross lag: 2.
    if(intercept){
      ind_intercept <- rep(3,M)
    }else ind_intercept <- NULL

    ind_matrix_small <- diag(M)
    ind_matrix_small[upper.tri(ind_matrix_small)] <-
      ind_matrix_small[lower.tri(ind_matrix_small)] <- 2
    ind_matrix <- rbind(matrix(rep(ind_matrix_small, p), p*M, M, byrow = TRUE),
                        ind_intercept)
    ind_vector <- as.vector(ind_matrix)
    n_1 <- length(which(ind_vector==1)) # number of own lag coefficients
    n_2 <- length(which(ind_vector==2)) # number of cross lag coefficients
    phi_prior_1 <- phi_prior[which(ind_vector==1)]
    phi_prior_2 <- phi_prior[which(ind_vector==2)]

    #vector of variances of univariate AR(6) models of the individual variables
    #scaling factors of prior variances (see Litterman 1986)
    sigma_sq <- MP_sigma_sq(Y = Yraw, p = 6, standardize = standardize)

    #preliminary prior variances of the coefficients (only scaling factors)
    V_i_p <- MP_V_prior_prep(sigma_sq=sigma_sq, K, M, intercept)

    # shape and rate parameters for Gamma hyperprior
    if(is.null(PHI_hyper)) {

      s1 <- s2 <- r1 <- r2 <- 0.01

    }else {

      s1 <- PHI_hyper$lambda_1_gamma[1]
      r1 <- PHI_hyper$lambda_1_gamma[2]
      s2 <- PHI_hyper$lambda_2_gamma[1]
      r2 <- PHI_hyper$lambda_2_gamma[2]

    }


  }else if(PHI_prior == "DL") {

    # concentration parameter of Dirichlet prior
    if(is.null(PHI_hyper)) {
      a <- 1/K
    }else if(PHI_hyper$hyperhyper == FALSE){
      a <- PHI_hyper$theta_dirichlet
    }else if(PHI_hyper$hyperhyper == TRUE){
      grid <- 1000
      a_tilde <- seq(1/(n),1/2,length.out = grid)
      a_tilde_mat <- matrix(rep(a_tilde,n), nrow = grid)
      prep1 <- a_tilde_mat - 1
      prep2 <- lgamma(rowSums(a_tilde_mat)) - rowSums(lgamma(a_tilde_mat))
      a <- 1/2
      #accept <- 0
      #sc <- 0.1
    }

  }else if(PHI_prior == "SSVS") {

    if(is.null(PHI_hyper)) {

      tau_0 <- 0.01 #spike
      tau_1 <- 100 #slab

      # shape parameters of BETA hyperprior on inclusion probabilities
      s_a <- 0.5
      s_b <- 0.5
      semi_automatic <- TRUE

    }else {

      tau_0 <- rep(PHI_hyper$tau_0, length.out = n)
      tau_1 <- rep(PHI_hyper$tau_1, length.out = n)
      s_a <- PHI_hyper$p_i_beta[1]
      s_b <- PHI_hyper$p_i_beta[2]
      semi_automatic <- PHI_hyper$semi_automatic

    }

    if(semi_automatic) {

      # Posterior mean of a flat Normal inverse Wishart prior
      V_prior_flat <- rep(10^3, K)
      V_post_flat <- tryCatch(solve(diag(1/V_prior_flat) + XX), error = function(e) chol2inv(chol(diag(1/V_prior_flat) + XX)))
      PHI_flat <- V_post_flat %*% (diag(1/V_prior_flat)%*%PHI0 + t(X)%*%Y)
      S_post <- diag(M) + crossprod(Y - X%*%PHI_flat) + t(PHI_flat - PHI0) %*%
        diag(1/V_prior_flat) %*% (PHI_flat - PHI0)
      Sigma_flat <- (S_post)/(M + 2 + T - M - 1)

      # variances of flat phi posterior estimates
      sigma_phi <- sqrt(diag(Sigma_flat %x% V_post_flat))

      # scale tau with variances
      tau_0 <- tau_0*sigma_phi
      tau_1 <- tau_1*sigma_phi
    }

  }else if(PHI_prior=="non_hierarchical"){

    if(is.null(PHI_hyper)) {

      V_i <- rep(10, n)

    }else V_i <- rep_len(PHI_hyper, length.out = n)
  }

  if(L_prior == "non_hierarchical") {

    if(is.null(L_hyper)) {

      V_i_L <- rep(10, n_l)

    }else V_i_L <- rep(L_hyper, n_l)

  }else if(L_prior == "DL") {

    # concentration parameter of Dirichlet prior
    if(is.null(L_hyper)) {
      b <- 1/n_l
    }else if(L_hyper$hyperhyper == FALSE){
      b <- L_hyper$theta_dirichlet
    }else if(L_hyper$hyperhyper == TRUE){
      b <- 1/2
      #accept_b <- 0
      #sc_b <- 0.1
      grid_b <- 1000
      b_tilde <- seq(1/(n_l),1/2,length.out = grid_b)
      b_tilde_mat <- matrix(rep(b_tilde,n_l), nrow = grid_b)
      prep1_b <- b_tilde_mat - 1
      prep2_b <- lgamma(rowSums(b_tilde_mat)) - rowSums(lgamma(b_tilde_mat))
    }

  }else if(L_prior == "SSVS" ) {

    if(is.null(L_hyper)) {

      kappa_0 <- rep(0.001, length.out=n_l) #spike
      kappa_1 <- rep(1, length.out=n_l) #slab

      # shape parameters of BETA hyperprior on inclusion probabilities
      sq_a <- 0.5
      sq_b <- 0.5

    }else {

      kappa_0 <- rep(L_hyper$kappa_0, length.out=n_l)
      kappa_1 <- rep(L_hyper$kappa_1, length.out=n_l)
      sq_a <- L_hyper$q_i_beta[1]
      sq_b <- L_hyper$q_i_beta[2]

    }

  }else if(L_prior == "MP") {

    if(is.null(L_hyper)) {
      s4 <- r4 <- 0.01
    }else {

      s4 <- L_hyper$lambda_4_gamma[1]
      r4 <- L_hyper$lambda_4_gamma[2]

    }


  }

  if(SV==FALSE){

    S0 <- 0.01   # prior scale of inverse Gamma (sigma2 of linear regressions)
    v0 <- 0.01   # prior shape parameter

  }else if(SV ==TRUE & is.null(SV_hyper)){

    SV_hyper <- stochvol::specify_priors(
      mu = stochvol::sv_normal(mean = 0, sd = 100),
      phi = stochvol::sv_beta(shape1 = 20, shape2 = 1.5),
      sigma2 = stochvol::sv_gamma(0.5,0.5),
      latent0_variance = "stationary"
    )

  }

  # Initial draws -----------------------------------------------------------

  # OLS estimates for PHI
  if(is.null(PHI_in)){
    PHI <- PHI_OLS <- tryCatch(solve(XX)%*%t(X)%*%Y,
                               error=function(e) MASS::ginv(XX)%*%t(X)%*%Y)
    PHI[PHI>-1e-100 & PHI <=0] <- -1e-100
    PHI[PHI < 1e-100 & PHI >=0] <- 1e-100

  }else{
    PHI <- PHI_in
  }
  phi <- as.vector(PHI)

  Ytilde <- Y - X %*% PHI
  if(is.null(L_in)){
    # OLS estimates for Sigma_t (and  L)
    SSR_OLS <- crossprod(Ytilde)
    Sigma <- Sigma_OLS <- if(T>K) SSR_OLS/(T-K) else SSR_OLS
    U <- try(chol(Sigma_OLS), silent = TRUE)
    if(!inherits(U, "try-error")){
      D <- diag(U)^2
      L_i <- U/sqrt(D)
      L <- backsolve(L_i, diag(M))
    }else{
      L <- diag(M)
      L[upper.tri(L)] <- rnorm(n_l)
    }
  }else{
    L <- L_in
  }


  l <- L[upper.tri(L)]
  if(SV == FALSE) {

    d <- matrix(1, T, M)

  }else if(SV == TRUE) {

    para <- matrix(data= c(rep(-10,M), rep(0.9,M), rep(0.2,M), rep(Inf, M),
                           rep(0,M), rep(NA, M), rep(-10, M)),
                   nrow = 7, ncol = M, byrow=TRUE)
    rownames(para) <- c("mu", "phi", "sigma", "nu", "rho", "beta", "latent0")
    h <- matrix(rep(-10, T*M), T,M)
    d <- exp(h)

  }

  if(PHI_prior == "MP") {

    lambdas <- rep(as.numeric(NA),n)
    # for numerical reasons (especially in large dimensions when K>>T)
    # lambda_1 and _2 will be kept fixed for the first few draws during the burnin
    lambda_1 <- lambdas[which(ind_vector==1)] <- 0.04 #rgamma(1,s1, r1)
    lambda_2 <- lambdas[which(ind_vector==2)] <- 0.0016 #rgamma(1, s2, r2)
    lambdas[which(ind_vector==3)] <- 10^3
    V_i <- lambdas * V_i_p

  }else if(PHI_prior == "DL") {

    theta <- rep(1/n, n)
    zeta <- rgamma(1, n*0.5, .5)

  }else if(PHI_prior == "SSVS" ){

    p_i <- 0.5 #rbeta(n, s_a, s_b)
    gammas <- rep(0,n) #rbinom(n, 1, p_i)

    d_i <- rep(as.numeric(NA),n) # d_i is tau_0 if gamma=0 and tau_1 if gamma=1
    d_i[which(gammas==0)] <- tau_0[which(gammas==0)]
    d_i[which(gammas==1)] <- tau_1[which(gammas==1)]

    V_i <-  d_i^2

  }

  if(L_prior == "DL") {

    theta_l <- rep(1/n_l, n_l)
    zeta_l <- rgamma(1, n_l*b, .5)

  }else if(L_prior =="SSVS") {

    q_i <- 0.5 #rbeta(n_l, sq_a, sq_b)
    omegas <- rep(0,n_l) #rbinom(n_l, 1, q_i)

    d_i_L <- rep(as.numeric(NA),n_l)
    d_i_L[which(omegas==0)] <- kappa_0[which(omegas==0)]
    d_i_L[which(omegas==1)] <- kappa_1[which(omegas==1)]

    V_i_L <-  d_i_L^2

  }

  start_time <- Sys.time()
  # Sampler -----------------------------------------------------------------

  for (rep in seq_len(tot)) {#seq_len(tot)
    setTxtProgressBar(pb,rep)

    ## Draw hyperparameters wrt PHI (get prior variances V_i)

    if(PHI_prior == "MP" & rep > 0.1*burnin) {

      # own lags
      V_i_1 <- V_i_p[which(ind_vector==1)]
      phi_1 <- phi[which(ind_vector==1)]

      lambda_1 <- GIGrvg::rgig(n=1, lambda = s1 - n_1/2,
                               chi = sum((phi_1 - phi_prior_1)^2/V_i_1),
                               psi = 2*r1 )

      # cross lags
      V_i_2 <- V_i_p[which(ind_vector==2)]
      phi_2 <- phi[which(ind_vector==2)]

      lambda_2 <- GIGrvg::rgig(n=1, lambda = s2 - n_2/2,
                               chi = sum((phi_2 - phi_prior_2)^2/V_i_2),
                               psi = 2*r2 )

      lambdas[which(ind_vector==1)] <- lambda_1
      lambdas[which(ind_vector==2)] <- lambda_2

      V_i <- lambdas * V_i_p

    }else if(PHI_prior == "DL") {


      if(!is.null(PHI_hyper)){

        if(PHI_hyper$hyperhyper==TRUE) {

          log_probs <- ddir2(theta, prep1 = prep1, prep2 = prep2, log = TRUE) +
            dgamma(zeta, shape = n*a_tilde, rate = 1/2, log = TRUE)

          w_i <- exp(log_probs - max(log_probs))
          weights <- w_i/sum(w_i)
          atmp <- as.vector(rmultinom(1, 1, weights))
          a <- a_tilde[which(atmp == 1)]

        }
      }

      psi <- as.vector(1/my_gig(n = 1, lambda = -.5, chi = 1,
                                psi = 1/(theta*zeta/abs(phi - phi_prior))^2))
      zeta<- GIGrvg::rgig(1, lambda = n*(a-1), psi = 1,
                          chi = 2*sum(abs(phi-phi_prior)/as.vector(theta)))
      w <- as.vector(my_gig(n = 1, lambda = a-1, chi = 2*abs(phi-phi_prior), psi = 1))
      theta <- as.vector(w/sum(w))

      V_i <- zeta^2*psi*theta^2

    }else if(PHI_prior == "SSVS"){

      if(rep > 0.1*burnin){

        #conditional posterior inclusion probabilities gst
        u_i1 <- dnorm(phi, phi_prior, (tau_1), log = TRUE) + log(p_i)
        u_i2 <- dnorm(phi, phi_prior, (tau_0), log = TRUE) + log(1 - p_i)
        logdif <- u_i2 - u_i1
        gst <- 1/(1 + exp(logdif))

        #update auxilliary dummy variable gamma
        gammas <- rbinom(n,1,gst)

        n0 <- length(which(gammas==0))
        selected <- which(gammas==1)
        notselected <- which(gammas == 0)

        #create vector of prior variances
        d_i[notselected] <- tau_0[notselected]
        d_i[selected] <- tau_1[selected]
        V_i <-  d_i^2

        #update prior inclusion probabilities
        p_i <- rbeta(n, s_a + gammas, s_b + 1 - gammas)

      }

    }

    ## Draw hyperparameters wrt L
    if(L_prior == "DL") {

      if(!is.null(L_hyper)){

        if(L_hyper$hyperhyper==TRUE) {

          log_probs_l <- ddir2(theta_l, prep1 = prep1_b, prep2 = prep2_b, log = TRUE) +
            dgamma(zeta_l, shape = n_l*b_tilde, rate = 1/2, log = TRUE)

          w_i_l <- exp(log_probs_l - max(log_probs_l))
          weights_l <- w_i_l/sum(w_i_l)
          btmp <- as.vector(rmultinom(1, 1, weights_l))
          b <- b_tilde[which(btmp == 1)]
        }
      }

      psi_l <- as.vector(1/my_gig(n = 1, lambda = -.5, chi = 1,
                                  psi = 1/(theta_l*zeta_l/abs(l))^2))
      zeta_l <- GIGrvg::rgig(1, lambda = n_l*(b-1), psi = 1,
                             chi = 2*sum(abs(l)/as.vector(theta_l)))
      w_l <- as.vector(my_gig(n = 1, lambda = b-1, chi = 2*abs(l), psi = 1))
      theta_l <- as.vector(w_l/sum(w_l))

      V_i_L <- zeta_l^2*psi_l*theta_l^2

    }else if(L_prior == "SSVS" & rep > 0.1*burnin){

      uu_i1 <- dnorm(l, 0, (kappa_1), log = TRUE) + log(q_i)
      uu_i2 <- dnorm(l, 0, (kappa_0), log = TRUE) + log(1 - q_i)

      logdif_l <- uu_i2 - uu_i1
      wght <- 1/(1 + exp(logdif_l))

      omegas <- rbinom(n_l,1,wght)

      selected_o <- which(omegas==1)
      notselected_o <- which(omegas == 0)
      n00 <- length(which(omegas==0))

      d_i_L[notselected_o] <- kappa_0[notselected_o]
      d_i_L[selected_o] <- kappa_1[selected_o]
      q_i <- rbeta(n_l, sq_a + omegas, sq_b + 1 - omegas)

      V_i_L <-  d_i_L^2

    }else if(L_prior == "MP") {

      lambda_4 <- GIGrvg::rgig(n=1, lambda = s4 - n_l/2,
                               chi = sum(l^2),
                               psi = 2*r4 )

      V_i_L <- rep(lambda_4, n_l)

    }

    ## Draw PHI
    PHI <- draw_PHI(PHI = PHI, PHI_prior = PHI0, Y = Y, X = X, L = L,
                    d = d, V_i = V_i, M = M, K = K)

    if(PHI_prior == "DL") {

      PHI[PHI >= 0 & PHI < 1e-100] <- 1e-100
      PHI[PHI <= 0 & PHI > -1e-100] <- -1e-100

    }

    phi <- as.vector(PHI)

    Ytilde <- Y - X %*% PHI

    ## Draw L
    L <- draw_L(Ytilde = Ytilde, V_i = V_i_L, d = d)

    if(L_prior == "DL") {

      L[upper.tri(L)][L[upper.tri(L)] >= 0 & L[upper.tri(L)] < 1e-100] <- 1e-100
      L[upper.tri(L)][L[upper.tri(L)] <= 0 & L[upper.tri(L)] > -1e-100] <- -1e-100

    }

    l <- L[upper.tri(L, diag = FALSE)]
    L_i <- backsolve(L, diag(M))

    Ytilde_L <- Ytilde %*% L

    ## Draw D_t
    if(SV == TRUE){

      for (i in seq_len(M)) {

        svdraw <- stochvol::svsample_fast_cpp(Ytilde_L[,i], draws = 1,
                                              burnin = 0, startpara = para[,i],
                                              startlatent = h[,i],
                                              priorspec = SV_hyper )
        h[,i] <- svdraw$latent
        para[1:5,i] <- svdraw$para
        para["latent0",i] <- svdraw$latent0
        d[,i] <- exp(h[,i])

      }


    }else if (SV == FALSE) {

      S_post <- S0 + 0.5*colSums(Ytilde_L^2)
      D <- 1/rgamma(M, shape=(v0+T)/2, rate =S_post)
      d <- matrix(rep(D, T), T, M, byrow = TRUE)

    }

    if(rep > burnin) {

      DRAWS$PHI$PHI_draws[rep-burnin,,] <- PHI
      DRAWS$L$L_draws[rep-burnin,,] <- L
      DRAWS$Ytilde_draws[rep-burnin,,] <- Ytilde
      V_i_draws[rep-burnin,] <- V_i

      if(PHI_prior == "MP") {

        DRAWS$PHI$lambda_1_draws[rep-burnin] <- lambda_1
        DRAWS$PHI$lambda_2_draws[rep-burnin] <- lambda_2

      }else if(PHI_prior == "DL") {

        DRAWS$PHI$psi_draws[rep-burnin,] <- psi
        DRAWS$PHI$zeta_draws[rep-burnin] <- zeta
        DRAWS$PHI$theta_draws[rep-burnin,] <- theta
        DRAWS$PHI$a_draws[rep-burnin] <- a

      }else if(PHI_prior =="SSVS" ) {

        DRAWS$PHI$gamma_draws[rep-burnin,] <- gammas
        DRAWS$PHI$p_i_draws[rep-burnin,] <- p_i

      }

      if(L_prior == "DL") {

        DRAWS$L$psi_l_draws[rep-burnin,] <- psi_l
        DRAWS$L$zeta_l_draws[rep-burnin] <- zeta_l
        DRAWS$L$theta_l_draws[rep-burnin,] <- theta_l
        DRAWS$L$b_draws[rep-burnin] <- b

      }else if(L_prior == "SSVS" ) {

        DRAWS$L$omega_draws[rep-burnin,] <- omegas
        DRAWS$L$q_i_draws[rep-burnin,] <- q_i

      }else if(L_prior == "MP") {

        DRAWS$L$lambda_4_draws[rep-burnin] <- lambda_4

      }

      if(SV==FALSE){

        DRAWS$D$D_draws[rep-burnin,] <- D

      }else if(SV==TRUE){

        DRAWS$D$SVpara_draws[rep-burnin,,] <- para[c("mu", "phi", "sigma", "latent0"),]
        DRAWS$D$SVlatent_draws[rep-burnin,,] <- h

      }
    }
  }##end of sampler

  timer <- Sys.time() - start_time
  close(pb)
  cat("Finished MCMC after ", format(round(timer, 2)), ".\n", sep = "")
  bench <- as.numeric(timer/tot)
  attributes(bench) <- list("names" = "secs/itr")
  # covariables for one-step ahead predictions
  if(intercept){
    X_fore <- c(as.vector(t(Y[T:(T-p+1),])), 1)
  }else X_fore <- as.vector(t(Y[T:(T-p+1),]))

  return(list(Posterior_draws = DRAWS,
              Y=Y,
              X=X,
              M=M,
              K=K,
              T=T,
              Yraw=Yraw,
              Traw=nrow(Yraw),
              p=p,
              variables=variables,
              X_fore=X_fore,
              #PHI_OLS=PHI_OLS,
              #Sigma_OLS=Sigma_OLS,
              mu_Y=mu_Y,
              sd_Y=sd_Y,
              standardize=standardize,
              intercept = intercept,
              bench=bench
  ))
}
