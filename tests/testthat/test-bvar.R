test_that("setup mainly PHI", {
  data("exrates", package = "stochvol")
  y <- factorstochvol::logret(exrates[1:15, 1:3], demean = TRUE)
  draws <- 27L
  burnin <- 10L
  
  thin_values <- c(1,3)
  lag_values <- c(0L,1L,3L)
  PHI_priors <- c("normal","HMP", "SSVS", "NG", "DL", "R2D2", "HS") #"HMP"
  intercept_values <- c(0,10)
  sv_keep <- c("all","last")
  semi_groups <- c("global", "olcl-lagwise")
  keep <- "last"
  error_spec <- c("cholesky", "factor") 
  ssvs_specs <- c(0.5, 1)
  gl_specs <- c(0.5, 1)
  # test prior_phi options
  res <- list(bvar(y,draws = draws, burnin = burnin,quiet=TRUE))
  for(thin in thin_values){
    #for(keep in sv_keep){
    for(lags in lag_values){
      for(errs in error_spec){
        prior_sigma <- specify_prior_sigma(data = y, type = errs, quiet = TRUE)
        for(intercept in intercept_values){
          myintercept <- if(intercept==0) FALSE else intercept
          for(prior in PHI_priors){
            for(group in semi_groups){
              if(!(prior %in% c("SSVS", "NG", "DL", "R2D2", "HS")) & which(semi_groups==group)>1) {
                next
              }
              for(ssvs_p in ssvs_specs){
                if(prior != "SSVS" & which(ssvs_specs==ssvs_p)>1) {
                  next
                }
                for(a_shrink in gl_specs){
                  if(!(prior %in% c("DL", "NG", "R2D2")) & which(gl_specs==a_shrink)>1) {
                    next
                  }
                  if(prior == "normal"){
                    prior_phi <- specify_prior_phi(data = y, lags = lags, prior = prior, normal_sds = 10)
                  }else if(prior == "HMP"){
                    prior_phi <- specify_prior_phi(data = y, lags = lags, prior = prior)
                  }else if(prior == "SSVS" || prior == "NG" || prior == "DL" ||
                           prior == "R2D2" || prior == "HS"){
                    if(prior == "SSVS"){
                      if(ssvs_p == 1){
                        ssvs_p = c(1,1)
                      }
                      prior_phi <- specify_prior_phi(data = y, lags = lags, prior = prior,
                                                     global_grouping = group,
                                                     SSVS_a = ssvs_p)
                    }else if(prior == "NG" || prior == "DL" || prior == "R2D2" ){
                      if(a_shrink==1){
                        a_shrink <- cbind(seq(0.1,1,by=.1), rep(1/10,10))
                      }
                      if(prior == "NG"){
                        prior_phi <- specify_prior_phi(data = y, lags = lags, prior = prior,
                                                       global_grouping = group,
                                                       NG_a = a_shrink)
                      }else if(prior == "R2D2"){
                        prior_phi <- specify_prior_phi(data = y, lags = lags, prior = prior,
                                                       global_grouping = group,
                                                       R2D2_a = a_shrink)
                      }else if(prior == "DL"){
                        prior_phi <- specify_prior_phi(data = y, lags = lags, prior = prior,
                                                       global_grouping = group,
                                                       DL_a = a_shrink)
                      }
                      
                    }else if(prior == "HS"){
                      prior_phi <- specify_prior_phi(data = y, lags = lags, prior = prior,
                                                     global_grouping = group)
                    }
                  }
                  expect_no_condition(bvar(y , lags = lags, prior_intercept = myintercept,
                                           draws = draws, burnin = burnin, thin = thin,
                                           prior_phi = prior_phi, prior_sigma = prior_sigma,
                                           sv_keep = keep,
                                           quiet = TRUE))
                }# a_shrink end
              }# ssvs_p end 
            }# group end
          }# prior end
        }# intercept end
      }# errs end
    }# lags end
  }# thin end
})

test_that("setup sigma cholesky", {
  data("exrates", package = "stochvol")
  y <- factorstochvol::logret(exrates[1:15, 1:3], demean = TRUE)
  draws <- 27L
  burnin <- 10L
  
  SIGMA_priors <- c("normal","HMP", "SSVS", "NG", "DL", "R2D2", "HS")
  heteroscedastic <- c(TRUE,FALSE)
  gl_specs <- c(0.5,1)
  ssvs_specs <- c(0.5,1)
  thin_values <- c(1,3)
  keep <- "last"
  
  for(thin in thin_values){
    for(scedastic in heteroscedastic){
      for(prior in SIGMA_priors){
        for(a_shrink in gl_specs){
          if(!(prior %in% c("NG", "DL", "R2D2")) & which(gl_specs==a_shrink)>1L) {
            next
          }
          for(ssvs_p in ssvs_specs){
            if(prior != "SSVS" & which(ssvs_specs == ssvs_p)>1L) {
              next
            }
            if(prior == "normal" |prior == "HMP" | prior =="HS"){
              prior_sigma <- specify_prior_sigma(data = y, type = "cholesky",
                                                 cholesky_U_prior = prior,
                                                 cholesky_heteroscedastic = scedastic,
                                                 quiet = TRUE)
            }else if(prior == "NG" | prior == "R2D2" | prior == "DL"){
              if(a_shrink==1){
                a_shrink <- cbind(seq(0.1,1,by=.1), rep(1/10,10))
              }
              if(prior == "NG"){
                prior_sigma <- specify_prior_sigma(data = y, type = "cholesky",
                                                   cholesky_U_prior = prior,
                                                   cholesky_heteroscedastic = scedastic,
                                                   cholesky_NG_a = a_shrink,
                                                   quiet = TRUE)
              }else if(prior == "R2D2"){
                prior_sigma <- specify_prior_sigma(data = y, type = "cholesky",
                                                   cholesky_U_prior = prior,
                                                   cholesky_heteroscedastic = scedastic,
                                                   cholesky_R2D2_a = a_shrink,
                                                   quiet = TRUE)
              }else if(prior == "DL"){
                prior_sigma <- specify_prior_sigma(data = y, type = "cholesky",
                                                   cholesky_U_prior = prior,
                                                   cholesky_heteroscedastic = scedastic,
                                                   cholesky_DL_a =  a_shrink,
                                                   quiet = TRUE)
              }
            }else if(prior == "SSVS" ){
              if(ssvs_p == 1){
                ssvs_p = c(1,1)
              }
              prior_sigma <- specify_prior_sigma(data = y, type = "cholesky",
                                                 cholesky_U_prior = prior,
                                                 cholesky_heteroscedastic = scedastic,
                                                 cholesky_SSVS_a = ssvs_p,
                                                 quiet = TRUE)
              
            }
            expect_no_condition(bvar(y , draws = draws, burnin = burnin, thin = thin,
                                     prior_sigma = prior_sigma,
                                     sv_keep = keep,
                                     quiet = TRUE))
          }
        }
      }
    }
  }
  
})

test_that("setup sigma factor", {
  data("exrates", package = "stochvol")
  y <- factorstochvol::logret(exrates[1:15, 1:3], demean = TRUE)
  draws <- 27L
  burnin <- 10L
  
  thin_values <- c(1,3)
  factors_values <- c(0, 1, 3)
  restrict_values <- list("upper", "none")
  priorfacloadtype_values <- c("normal", "rowwiseng", "colwiseng")
  priorhomoskedastic <- matrix(c(1.1, 1.1), nrow = NCOL(y),
                               ncol = 2, byrow = TRUE)
  heteroskedastic_values <- list(TRUE, c(FALSE, FALSE))
  huge_algo <- c(TRUE, FALSE)
  
  for (th in thin_values) {
    for(ha in huge_algo){
      for (pflt in priorfacloadtype_values) {
        for (hsk in heteroskedastic_values) {
          for (fs in factors_values) {
            if (fs > 1) {
              for (rst in restrict_values) {
                prior_sigma <- specify_prior_sigma(data = y, type = "factor", 
                                                   factor_factors = fs,
                                                   factor_priorfacloadtype = pflt,
                                                   factor_restrict = rst,
                                                   factor_heteroskedastic = hsk,
                                                   factor_priorhomoskedastic = if (!isTRUE(hsk)) priorhomoskedastic else NA,
                                                   factor_interweaving = if (!isTRUE(hsk)) 0 else 4,
                                                   quiet = TRUE)
              }
            }else{
              prior_sigma <- specify_prior_sigma(data = y, type = "factor", 
                                                 factor_factors = fs,
                                                 factor_priorfacloadtype = pflt,
                                                 factor_restrict = "none",
                                                 factor_heteroskedastic = hsk,
                                                 factor_priorhomoskedastic = if (!isTRUE(hsk)) priorhomoskedastic else NA,
                                                 factor_interweaving = if (!isTRUE(hsk)) 0 else 4,
                                                 quiet = TRUE)
            }
            expect_no_condition(bvar(y,draws = draws, burnin = burnin, thin = th, 
                                     prior_sigma = prior_sigma, quiet=TRUE,
                                     expert_huge = ha))
          }
        }
      } 
    }
  }
  
})

test_that("flat prior cholesky", {
  data <- usmacro_growth[,c("GDPC1","CPIAUCSL","FEDFUNDS")]
  phi <- specify_prior_phi(data = data, lags = 2L, prior = "normal",
                           normal_sds = 1e6)
  sigma <- specify_prior_sigma(data = data, type = "cholesky",
                               cholesky_U_prior = "normal",
                               cholesky_normal_sds = 1e6,
                               cholesky_heteroscedastic = FALSE,
                               cholesky_priorhomoscedastic = matrix(c(1e-6,1e-6), ncol(data), 2),
                               quiet = TRUE)
  set.seed(123)
  mod <- bvar(data, lags = 2, prior_intercept = 1e6, prior_phi = phi,
              prior_sigma = sigma, draws = 10000, quiet = TRUE)
  phi_post_mean <- apply(mod$PHI, 1:2, mean)
  ols <- solve(crossprod(mod$X), crossprod(mod$X,mod$Y))

  expect_lt(max(abs(ols-phi_post_mean)),0.01)
})

test_that("geweke test sample_phi_cholesky", {

  # distribution test function
  mydist.test <- function(x,y){
    # x: iid draws
    # y: autocorrelated draws
    x.n <- length(x)
    y.n <- length(y)
    x.mean <- mean(x)
    y.mean <- mean(y)
    x.var <- var(x)
    y.var <- coda::spectrum0.ar(y)$spec
    statistic <- (x.mean - y.mean)/sqrt(x.var/x.n + y.var/y.n)
    pnorm(-abs(statistic))
  }

  # settings
  set.seed(123)
  draws <- 10000
  n <- 50 # observations
  M <- 5 # time-series
  K <- 10 # coefficients per equation

  # 'known' variance-covariance
  U_inv <- diag(M)
  U_inv[upper.tri(U_inv)] <- seq(.1, by = .1, len = (M^2-M)/2)
  U <- backsolve(U_inv, diag(M))
  d <- rep(20, M)
  d_sqrt <- sqrt(d)
  d_sqrt_mat <- matrix(d_sqrt, n, M, byrow = TRUE)
  SIGMA <- crossprod(U_inv*d_sqrt)
  SIGMA_chol <- diag(d_sqrt)%*%U_inv
  #cov2cor(SIGMA)

  # simulate regressors
  X <- matrix(rnorm(n*K), n, K)

  # Prior: coefficients ~ (iid) N(0,1)
  V_prior <- matrix(1, K, M)
  PHI_prior <- matrix(0, K, M)

  # independent draws
  PHI_ind <- array(rnorm(K*M*draws,as.vector(PHI_prior), as.vector(sqrt(V_prior))),c(K,M,draws))

  # MCMC draws
  # initialize with draw from prior
  PHI_draws <- array(as.numeric(NA), c( K, M, draws))
  PHI <- matrix(rnorm(K*M, 0, sqrt(as.vector(V_prior))), K, M)
  # alternately draw from observables|unobservables and unobservables|oberservables
  for(r in seq.int(draws)){
    # simulate observables|unobservables
    Y <- X%*%PHI + matrix(rnorm(n*M), n, M, byrow = TRUE)%*%SIGMA_chol

    # simulate unobservables|oberservables
    PHI_draws[,,r] <- PHI <- bayesianVARs:::sample_PHI_cholesky(PHI, PHI_prior, Y,
                                                                X, U, d_sqrt_mat,
                                                                V_prior,
                                                                0,
                                                                FALSE,
                                                                1)
  }

  # compare distribution of quadratic length in terms of euclidean distance of iid and autocorrelated draws
  test1 <- mydist.test(apply(PHI_ind, 3, function(x) (sum(x^2))),apply(PHI_draws, 3, function(x) (sum(x^2))))
  #qqplot(apply(PHI_ind, 3, function(x) (sum(x^2))), apply(PHI_draws, 3, function(x) (sum(x^2))));abline(0,1)


  test2 <-  ks.test(apply(PHI_draws, 3, function(x) (sum(x^2))), "pchisq", df=K*M)$p.value
  #qqplot(qchisq(ppoints(draws), df = K*M), apply(PHI_draws[,,], 3, function(x) (sum(x^2))),
  #main = expression("Q-Q plot for" ~~ {chi^2}[nu == K*M]));abline(0,1)
  expect_gt(test1, 0.01)
  expect_gt(test2, 0.01)

})

test_that("miss-specified input", {
  data <- usmacro_growth[,c("GDPC1","CPIAUCSL","FEDFUNDS")]
  phi <- specify_prior_phi(data = data, lags = 4L)
  expect_error(bvar(data, lags = 1L, prior_phi = phi))
})

test_that("flat prior factor", {
  data <- usmacro_growth[,c("GDPC1","CPIAUCSL","FEDFUNDS")]
  phi <- specify_prior_phi(data = data, lags = 2L, prior = "normal",
                           normal_sds = 1e6)
  factor_prior_homoscedastic <- matrix(c(1e-6,1e-6), ncol(data), 2)
  facload_priorsd <- 1e6
  sigma <- specify_prior_sigma(data = data, type = "factor",
                               factor_factors = 3,
                               factor_priorfacloadtype = "normal",
                               factor_priorfacload = facload_priorsd,
                               factor_heteroskedastic = c(FALSE, FALSE),
                               factor_priorhomoskedastic = factor_prior_homoscedastic,
                               quiet = TRUE)
  # factor_prior_homoscedastic
  expect_equal(sigma$prior_sigma_cpp$factor_priorhomoskedastic, factor_prior_homoscedastic)
  # check that prior standard deviations are forwarded correctly
  expect_true(all(sigma$general_settings$factor_starttau2 == facload_priorsd^2))
  # normal prior implies that factor_ngprior==FALSE
  expect_false(sigma$prior_sigma_cpp$factor_ngprior)
  skip_on_cran()
  # the following caused a valgrind issue on CRAN which I could not replicate
  set.seed(123)
  mod <- bvar(data, lags = 2, prior_intercept = 1e6, prior_phi = phi,
              prior_sigma = sigma, draws = 10000, quiet = TRUE)
  phi_post_mean <- apply(mod$PHI, 1:2, mean)
  ols <- solve(crossprod(mod$X), crossprod(mod$X,mod$Y))

  expect_lt(max(abs(ols-phi_post_mean)),0.015)
  # check that factor logvars are all zero
  expect_true(all(mod$logvar[,-seq_len(ncol(mod$Y)),]==0))
})

test_that("iterative_use", {
  data <- 100*usmacro_growth[,c("GDPC1","CPIAUCSL","FEDFUNDS")]
  p <- 2L
  
  draws <-5L
  burnin <- 5L
  # "HS", "DL", "NG", "R2D2", "SSVS", "normal", "HMP"
  for(shrink in c("HS", "DL", "NG", "R2D2", "SSVS", "normal", "HMP")) {
    for(grouping in c("global", "olcl-lagwise")) {
      if(shrink %in% c("normal", "HMP") & grouping == "olcl-lagwise"){
        next
      }
      for(afix in c(TRUE, FALSE)){
        if(shrink %in% c("normal", "HMP", "HS") & base::isFALSE(afix)) {
          next
        }
        thea <- numeric(1L)
        if(afix) {
          thea <- 0.5
        } else {
          if(shrink == "SSVS") {
            thea <- c(1,1)
          } else if (shrink %in% c("DL", "NG", "R2D2")) {
            thea <- cbind(seq(0.1,1, length.out = 10L), 1)
          }
        }
      prior_phi <- specify_prior_phi(data, lags = p, prior = shrink,
                                     global_grouping = grouping,
                                     SSVS_a = thea,
                                     DL_a = thea,
                                     NG_a = thea,
                                     R2D2_a = thea
                                     )
        for(type in c("factor", "cholesky")){
          use_update <- sample(c(TRUE,FALSE), 1)
          prior_sigma <- specify_prior_sigma(data, type = type,
                                             quiet = TRUE,
                                             cholesky_U_prior = shrink,
                                             cholesky_NG_a = thea,
                                             cholesky_DL_a = thea,
                                             cholesky_R2D2_a = thea)
          set.seed(1)
          res0 <- bvar(data, lags = p, draws = draws, burnin = burnin,
                       prior_phi = prior_phi, prior_sigma = prior_sigma,
                       quiet = TRUE)
          testseed0 <- .Random.seed
          set.seed(1)
          res2_burn <- bvar(data, lags = p, draws = 0L, burnin = burnin,
                       prior_phi = prior_phi, prior_sigma = prior_sigma, quiet = TRUE)
          res2 <- if(use_update) {
            update_bvar(res2_burn, draws = draws, quiet = TRUE)
          } else {
            bvar(data, lags = p, draws = draws, burnin = 0L,
                 prior_phi = prior_phi, prior_sigma = prior_sigma,
                 startvals = res2_burn[["lastdraw"]], 
                 quiet = TRUE) 
          }
          testseed2 <- .Random.seed
          
          expect_equal(res2[["lastdraw"]], res0[["lastdraw"]])
          expect_identical(testseed2, testseed0)
        }
      }
    }  
  }
})

test_that("stable VAR", {
  data <- 100*usmacro_growth[,c("GDPC1","CPIAUCSL","FEDFUNDS")]
  lags <- c(1L,3L)
  for(p in lags) {
    prior_phi <- specify_prior_phi(data = data, lags = p, stable = TRUE, stable_lim = 0.99)
    res <- bvar(data = data, lags = p, prior_phi = prior_phi, burnin = 0L, draws = 100L,
                quiet = TRUE)
    # stable_bvar is not expected to detect draws that do not fulfill the stable
    # criterion
    res_stable_post <- stable_bvar(res, lim = 0.99, quiet = TRUE)
    expect_equal(res, res_stable_post) 
  }
})
