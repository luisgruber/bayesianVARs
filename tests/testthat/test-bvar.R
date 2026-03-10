
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
                                                                V_prior)
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

  expect_lt(max(abs(ols-phi_post_mean)),0.01)
  # check that factor logvars are all zero
  expect_true(all(mod$logvar[,-seq_len(ncol(mod$Y)),]==0))
})

test_that("expert_huge", {
  n <- 20
  m <- 10
  y <- matrix(rnorm(n*m), n, m)
  colnames(y) <- paste0("y", 1:m)
  expect_warning(bvar(y, draws = 2L, burnin = 2L, quiet = TRUE,
                      prior_sigma = specify_prior_sigma(y, type = "cholesky"),
                      expert_huge = TRUE))
})
