test_that("miss-specified input", {
  data <- usmacro_growth[,c("GDPC1","CPIAUCSL","FEDFUNDS")]
  phi <- specify_prior_phi(data = data, lags = 4L)
  expect_error(bvar(data, lags = 1L, prior_phi = phi))
})

test_that("flat prior factor", {
  data <- usmacro_growth[,c("GDPC1","CPIAUCSL","FEDFUNDS")]
  phi <- specify_prior_phi(data = data, lags = 2L, prior = "normal",
                           normal_sds = 1e6)
  sigma <- specify_prior_sigma(data = data, type = "factor",
                               factor_factors = 2,
                               factor_priorfacloadtype = "normal",
                               factor_priorfacload = 1e6,
                               factor_heteroskedastic = c(FALSE, FALSE),
                               factor_priorhomoskedastic = matrix(c(0.01,.01), ncol(data), 2))
 set.seed(123)
  mod <- bvar(data, lags = 2, prior_intercept = 1e6, prior_phi = phi,
              prior_sigma = sigma, draws = 10000)
  phi_post_mean <- apply(mod$PHI, 1:2, mean)
  ols <- solve(crossprod(mod$X), crossprod(mod$X,mod$Y))

  expect_lt(max(abs(ols-phi_post_mean)),0.01)
})

test_that("flat prior cholesky", {
  data <- usmacro_growth[,c("GDPC1","CPIAUCSL","FEDFUNDS")]
  phi <- specify_prior_phi(data = data, lags = 2L, prior = "normal",
                           normal_sds = 1e6)
  sigma <- specify_prior_sigma(data = data, type = "cholesky",
                               cholesky_U_prior = "normal",
                               cholesky_normal_sds = 1e6,
                               cholesky_heteroscedastic = FALSE,
                               cholesky_priorhomoscedastic = matrix(c(0.01,0.01), ncol(data), 2))
  set.seed(123)
  mod <- bvar(data, lags = 2, prior_intercept = 1e6, prior_phi = phi,
              prior_sigma = sigma, draws = 10000)
  phi_post_mean <- apply(mod$PHI, 1:2, mean)
  ols <- solve(crossprod(mod$X), crossprod(mod$X,mod$Y))

  expect_lt(max(abs(ols-phi_post_mean)),0.01)
})
