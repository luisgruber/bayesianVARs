
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

test_that("R2D2", {
  M <-  3L
  lags <- 2L
  sup <- seq(.1,1,.1)
  a_fix <- 0.1
  b_fix <- 0.4
  a_mat <- cbind(sup, dexp(sup))

  # R2D2, a fixed, global
  r2d2_afix_global <- specify_prior_phi(M = M, lags = lags, prior = "R2D2",
                                        R2D2_a = a_fix, R2D2_b = b_fix,
                                        global_grouping = "global")
  # R2D2, a random, global
  r2d2_ah_global <- specify_prior_phi(M = M, lags = lags, prior = "R2D2",
                                      R2D2_a = a_mat, R2D2_b = b_fix,
                                      global_grouping = "global")
  # R2D2, a fixed, semi_global
  r2d2_afix_semiglobal <- specify_prior_phi(M = M, lags = lags, prior = "R2D2",
                                            R2D2_a = a_fix, R2D2_b = b_fix,
                                            global_grouping = "olcl-lagwise")
  # R2D2, a random, semi-global
  r2d2_ah_semiglobal <- specify_prior_phi(M = M, lags = lags, prior = "R2D2",
                                          R2D2_a = a_mat, R2D2_b = b_fix,
                                          global_grouping = "olcl-lagwise")

  # GT_vs == 1/2
  expect_equal(r2d2_afix_global[["prior_phi_cpp"]][["GT_vs"]], 1/2)
  expect_equal(r2d2_afix_semiglobal[["prior_phi_cpp"]][["GT_vs"]], 1/2)
  expect_equal(r2d2_ah_global[["prior_phi_cpp"]][["GT_vs"]], 1/2)
  expect_equal(r2d2_ah_semiglobal[["prior_phi_cpp"]][["GT_vs"]], 1/2)
  # GT_priorkernel == "exponential"
  expect_identical(r2d2_afix_global[["prior_phi_cpp"]][["GT_priorkernel"]], "exponential")
  expect_identical(r2d2_afix_semiglobal[["prior_phi_cpp"]][["GT_priorkernel"]], "exponential")
  expect_identical(r2d2_ah_global[["prior_phi_cpp"]][["GT_priorkernel"]], "exponential")
  expect_identical(r2d2_ah_semiglobal[["prior_phi_cpp"]][["GT_priorkernel"]], "exponential")
  # a = a_fix
  expect_equal(r2d2_afix_global[["prior_phi_cpp"]][["a"]], a_fix)
  expect_equal(r2d2_afix_semiglobal[["prior_phi_cpp"]][["a"]], rep(a_fix,r2d2_afix_semiglobal[["prior_phi_cpp"]][["n_groups"]]))
  # b = b_fix
  expect_equal(r2d2_afix_global[["prior_phi_cpp"]][["b"]], b_fix)
  expect_equal(r2d2_afix_semiglobal[["prior_phi_cpp"]][["b"]], rep(b_fix,r2d2_afix_semiglobal[["prior_phi_cpp"]][["n_groups"]]))
  expect_equal(r2d2_afix_global[["prior_phi_cpp"]][["b"]], b_fix)
  expect_equal(r2d2_afix_semiglobal[["prior_phi_cpp"]][["b"]], rep(b_fix,r2d2_afix_semiglobal[["prior_phi_cpp"]][["n_groups"]]))
  # c_rel_a == TRUE
  expect_true(r2d2_afix_global[["prior_phi_cpp"]][["c_rel_a"]])
  expect_true(r2d2_afix_semiglobal[["prior_phi_cpp"]][["c_rel_a"]])
  expect_true(r2d2_ah_global[["prior_phi_cpp"]][["c_rel_a"]])
  expect_true(r2d2_ah_semiglobal[["prior_phi_cpp"]][["c_rel_a"]])
  # c == 0.5*a
  expect_equal(2*r2d2_afix_global[["prior_phi_cpp"]][["c"]], r2d2_afix_global[["prior_phi_cpp"]][["a"]])
  expect_equal(2*r2d2_afix_semiglobal[["prior_phi_cpp"]][["c"]], r2d2_afix_semiglobal[["prior_phi_cpp"]][["a"]])
  expect_equal(2*r2d2_ah_global[["prior_phi_cpp"]][["c"]], r2d2_ah_global[["prior_phi_cpp"]][["a"]])
  expect_equal(2*r2d2_ah_semiglobal[["prior_phi_cpp"]][["c"]], r2d2_ah_semiglobal[["prior_phi_cpp"]][["a"]])
  # c_vec = 0.5*a_vec (only for "a" random)
  expect_equal(2*r2d2_ah_global[["prior_phi_cpp"]][["c_vec"]], r2d2_ah_global[["prior_phi_cpp"]][["a_vec"]])
  expect_equal(2*r2d2_ah_semiglobal[["prior_phi_cpp"]][["c_vec"]], r2d2_ah_semiglobal[["prior_phi_cpp"]][["a_vec"]])
  # a,b,c for each group
  expect_equal(r2d2_afix_global[["prior_phi_cpp"]][["n_groups"]], 1L)
  expect_equal(r2d2_afix_semiglobal[["prior_phi_cpp"]][["n_groups"]], length(r2d2_afix_semiglobal[["prior_phi_cpp"]][["groups"]]))
  expect_equal(r2d2_afix_semiglobal[["prior_phi_cpp"]][["n_groups"]], length(unique(as.vector(r2d2_afix_semiglobal[["general_settings"]][["i_mat"]]))))
  expect_equal(length(r2d2_afix_global[["prior_phi_cpp"]][["a"]]), r2d2_afix_global[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(r2d2_afix_semiglobal[["prior_phi_cpp"]][["a"]]), r2d2_afix_semiglobal[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(r2d2_ah_global[["prior_phi_cpp"]][["a"]]), r2d2_ah_global[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(r2d2_ah_semiglobal[["prior_phi_cpp"]][["a"]]), r2d2_ah_semiglobal[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(r2d2_afix_global[["prior_phi_cpp"]][["b"]]), r2d2_afix_global[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(r2d2_afix_semiglobal[["prior_phi_cpp"]][["b"]]), r2d2_afix_semiglobal[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(r2d2_ah_global[["prior_phi_cpp"]][["b"]]), r2d2_ah_global[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(r2d2_ah_semiglobal[["prior_phi_cpp"]][["b"]]), r2d2_ah_semiglobal[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(r2d2_afix_global[["prior_phi_cpp"]][["c"]]), r2d2_afix_global[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(r2d2_afix_semiglobal[["prior_phi_cpp"]][["c"]]), r2d2_afix_semiglobal[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(r2d2_ah_global[["prior_phi_cpp"]][["c"]]), r2d2_ah_global[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(r2d2_ah_semiglobal[["prior_phi_cpp"]][["c"]]), r2d2_ah_semiglobal[["prior_phi_cpp"]][["n_groups"]])
  # GT_hyper
  expect_true(r2d2_ah_global[["prior_phi_cpp"]][["GT_hyper"]])
  expect_true(r2d2_ah_semiglobal[["prior_phi_cpp"]][["GT_hyper"]])
  expect_false(r2d2_afix_global[["prior_phi_cpp"]][["GT_hyper"]])
  expect_false(r2d2_afix_semiglobal[["prior_phi_cpp"]][["GT_hyper"]])
})

test_that("NG", {
  M <-  3L
  lags <- 2L
  sup <- seq(.1,1,.1)
  a_fix <- 0.1
  b_fix <- 0.4
  c_fix <- 0.2
  c_rel_a <- "0.3*a"
  a_mat <- cbind(sup, dexp(sup))

  # NG, a fixed, global
  ng_afix_global <- specify_prior_phi(M = M, lags = lags, prior = "NG",
                                      NG_a = a_fix, NG_b = b_fix, NG_c = c_fix,
                                      global_grouping = "global")
  # NG, a random, global
  ng_ah_global <- specify_prior_phi(M = M, lags = lags, prior = "NG",
                                    NG_a = a_mat, NG_b = b_fix, NG_c = c_fix,
                                    global_grouping = "global")
  # NG, a fixed, semi_global
  ng_afix_semiglobal <- specify_prior_phi(M = M, lags = lags, prior = "NG",
                                          NG_a = a_fix, NG_b = b_fix, NG_c = c_fix,
                                          global_grouping = "olcl-lagwise")
  # NG, a random, semi-global
  ng_ah_semiglobal <- specify_prior_phi(M = M, lags = lags, prior = "NG",
                                        NG_a = a_mat, NG_b = b_fix, NG_c = c_fix,
                                        global_grouping = "olcl-lagwise")
  # NG, a fixed, global, c_rel_a
  ng_afix_global_crela <- specify_prior_phi(M = M, lags = lags, prior = "NG",
                                      NG_a = a_fix, NG_b = b_fix, NG_c = c_rel_a,
                                      global_grouping = "global")
  # NG, a random, global, c_rel_a
  ng_ah_global_crela <- specify_prior_phi(M = M, lags = lags, prior = "NG",
                                    NG_a = a_mat, NG_b = b_fix, NG_c = c_rel_a,
                                    global_grouping = "global")
  # NG, a fixed, semi_global, c_rel_a
  ng_afix_semiglobal_crela <- specify_prior_phi(M = M, lags = lags, prior = "NG",
                                          NG_a = a_fix, NG_b = b_fix, NG_c = c_rel_a,
                                          global_grouping = "olcl-lagwise")
  # NG, a random, semi-global, c_rel_a
  ng_ah_semiglobal_c_rel_a <- specify_prior_phi(M = M, lags = lags, prior = "NG",
                                        NG_a = a_mat, NG_b = b_fix, NG_c = c_rel_a,
                                        global_grouping = "olcl-lagwise")

  # GT_vs == 1
  expect_equal(ng_afix_global[["prior_phi_cpp"]][["GT_vs"]], 1)
  expect_equal(ng_afix_semiglobal[["prior_phi_cpp"]][["GT_vs"]], 1)
  expect_equal(ng_ah_global[["prior_phi_cpp"]][["GT_vs"]], 1)
  expect_equal(ng_ah_semiglobal[["prior_phi_cpp"]][["GT_vs"]], 1)
  expect_equal(ng_afix_global_crela[["prior_phi_cpp"]][["GT_vs"]], 1)
  expect_equal(ng_afix_semiglobal_crela[["prior_phi_cpp"]][["GT_vs"]], 1)
  expect_equal(ng_ah_global_crela[["prior_phi_cpp"]][["GT_vs"]], 1)
  expect_equal(ng_ah_semiglobal_c_rel_a[["prior_phi_cpp"]][["GT_vs"]], 1)
  # GT_priorkernel == "normal"
  expect_identical(ng_afix_global[["prior_phi_cpp"]][["GT_priorkernel"]], "normal")
  expect_identical(ng_afix_semiglobal[["prior_phi_cpp"]][["GT_priorkernel"]], "normal")
  expect_identical(ng_ah_global[["prior_phi_cpp"]][["GT_priorkernel"]], "normal")
  expect_identical(ng_ah_semiglobal[["prior_phi_cpp"]][["GT_priorkernel"]], "normal")
  expect_identical(ng_afix_global_crela[["prior_phi_cpp"]][["GT_priorkernel"]], "normal")
  expect_identical(ng_afix_semiglobal_crela[["prior_phi_cpp"]][["GT_priorkernel"]], "normal")
  expect_identical(ng_ah_global_crela[["prior_phi_cpp"]][["GT_priorkernel"]], "normal")
  expect_identical(ng_ah_semiglobal_c_rel_a[["prior_phi_cpp"]][["GT_priorkernel"]], "normal")
  # a = a_fix
  expect_equal(ng_afix_global[["prior_phi_cpp"]][["a"]], a_fix)
  expect_equal(ng_afix_semiglobal[["prior_phi_cpp"]][["a"]], rep(a_fix,ng_afix_semiglobal[["prior_phi_cpp"]][["n_groups"]]))
  expect_equal(ng_afix_global_crela[["prior_phi_cpp"]][["a"]], a_fix)
  expect_equal(ng_afix_semiglobal_crela[["prior_phi_cpp"]][["a"]], rep(a_fix,ng_afix_semiglobal_crela[["prior_phi_cpp"]][["n_groups"]]))
  # b = b_fix
  expect_equal(ng_afix_global[["prior_phi_cpp"]][["b"]], b_fix)
  expect_equal(ng_afix_semiglobal[["prior_phi_cpp"]][["b"]], rep(b_fix,ng_afix_semiglobal[["prior_phi_cpp"]][["n_groups"]]))
  expect_equal(ng_afix_global[["prior_phi_cpp"]][["b"]], b_fix)
  expect_equal(ng_afix_semiglobal[["prior_phi_cpp"]][["b"]], rep(b_fix,ng_afix_semiglobal[["prior_phi_cpp"]][["n_groups"]]))
  expect_equal(ng_afix_global_crela[["prior_phi_cpp"]][["b"]], b_fix)
  expect_equal(ng_afix_semiglobal_crela[["prior_phi_cpp"]][["b"]], rep(b_fix,ng_afix_semiglobal_crela[["prior_phi_cpp"]][["n_groups"]]))
  expect_equal(ng_afix_global_crela[["prior_phi_cpp"]][["b"]], b_fix)
  expect_equal(ng_afix_semiglobal_crela[["prior_phi_cpp"]][["b"]], rep(b_fix,ng_afix_semiglobal_crela[["prior_phi_cpp"]][["n_groups"]]))
  # c_rel_a == TRUE
  expect_false(ng_afix_global[["prior_phi_cpp"]][["c_rel_a"]])
  expect_false(ng_afix_semiglobal[["prior_phi_cpp"]][["c_rel_a"]])
  expect_false(ng_ah_global[["prior_phi_cpp"]][["c_rel_a"]])
  expect_false(ng_ah_semiglobal[["prior_phi_cpp"]][["c_rel_a"]])
  expect_true(ng_afix_global_crela[["prior_phi_cpp"]][["c_rel_a"]])
  expect_true(ng_afix_semiglobal_crela[["prior_phi_cpp"]][["c_rel_a"]])
  expect_true(ng_ah_global_crela[["prior_phi_cpp"]][["c_rel_a"]])
  expect_true(ng_ah_semiglobal_c_rel_a[["prior_phi_cpp"]][["c_rel_a"]])
  # c == c_rel_a
  q <- as.numeric(strsplit(c_rel_a, "\\*")[[1]][1])
  expect_equal(ng_afix_global_crela[["prior_phi_cpp"]][["c"]], q*ng_afix_global_crela[["prior_phi_cpp"]][["a"]])
  expect_equal(ng_afix_semiglobal_crela[["prior_phi_cpp"]][["c"]], q*ng_afix_semiglobal_crela[["prior_phi_cpp"]][["a"]])
  expect_equal(ng_ah_global_crela[["prior_phi_cpp"]][["c"]], q*ng_ah_global_crela[["prior_phi_cpp"]][["a"]])
  expect_equal(ng_ah_semiglobal_c_rel_a[["prior_phi_cpp"]][["c"]], q*ng_ah_semiglobal_c_rel_a[["prior_phi_cpp"]][["a"]])
  # c_vec = q*a_vec (only for "a" random)
  expect_equal(ng_ah_global_crela[["prior_phi_cpp"]][["c_vec"]], q*ng_ah_global_crela[["prior_phi_cpp"]][["a_vec"]])
  expect_equal(ng_ah_semiglobal_c_rel_a[["prior_phi_cpp"]][["c_vec"]], q*ng_ah_semiglobal_c_rel_a[["prior_phi_cpp"]][["a_vec"]])
  # a,b,c for each group
  expect_equal(ng_afix_global[["prior_phi_cpp"]][["n_groups"]], 1L)
  expect_equal(ng_afix_semiglobal[["prior_phi_cpp"]][["n_groups"]], length(ng_afix_semiglobal[["prior_phi_cpp"]][["groups"]]))
  expect_equal(ng_afix_semiglobal[["prior_phi_cpp"]][["n_groups"]], length(unique(as.vector(ng_afix_semiglobal[["general_settings"]][["i_mat"]]))))
  expect_equal(length(ng_afix_global[["prior_phi_cpp"]][["a"]]), ng_afix_global[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_afix_semiglobal[["prior_phi_cpp"]][["a"]]), ng_afix_semiglobal[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_ah_global[["prior_phi_cpp"]][["a"]]), ng_ah_global[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_ah_semiglobal[["prior_phi_cpp"]][["a"]]), ng_ah_semiglobal[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_afix_global[["prior_phi_cpp"]][["b"]]), ng_afix_global[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_afix_semiglobal[["prior_phi_cpp"]][["b"]]), ng_afix_semiglobal[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_ah_global[["prior_phi_cpp"]][["b"]]), ng_ah_global[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_ah_semiglobal[["prior_phi_cpp"]][["b"]]), ng_ah_semiglobal[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_afix_global[["prior_phi_cpp"]][["c"]]), ng_afix_global[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_afix_semiglobal[["prior_phi_cpp"]][["c"]]), ng_afix_semiglobal[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_ah_global[["prior_phi_cpp"]][["c"]]), ng_ah_global[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_ah_semiglobal[["prior_phi_cpp"]][["c"]]), ng_ah_semiglobal[["prior_phi_cpp"]][["n_groups"]])

  expect_equal(length(ng_afix_global_crela[["prior_phi_cpp"]][["a"]]), ng_afix_global_crela[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_afix_semiglobal_crela[["prior_phi_cpp"]][["a"]]), ng_afix_semiglobal_crela[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_ah_global_crela[["prior_phi_cpp"]][["a"]]), ng_ah_global_crela[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_ah_semiglobal_c_rel_a[["prior_phi_cpp"]][["a"]]), ng_ah_semiglobal_c_rel_a[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_afix_global_crela[["prior_phi_cpp"]][["b"]]), ng_afix_global_crela[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_afix_semiglobal_crela[["prior_phi_cpp"]][["b"]]), ng_afix_semiglobal_crela[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_ah_global_crela[["prior_phi_cpp"]][["b"]]), ng_ah_global_crela[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_ah_semiglobal_c_rel_a[["prior_phi_cpp"]][["b"]]), ng_ah_semiglobal_c_rel_a[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_afix_global_crela[["prior_phi_cpp"]][["c"]]), ng_afix_global_crela[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_afix_semiglobal_crela[["prior_phi_cpp"]][["c"]]), ng_afix_semiglobal_crela[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_ah_global_crela[["prior_phi_cpp"]][["c"]]), ng_ah_global_crela[["prior_phi_cpp"]][["n_groups"]])
  expect_equal(length(ng_ah_semiglobal_c_rel_a[["prior_phi_cpp"]][["c"]]), ng_ah_semiglobal_c_rel_a[["prior_phi_cpp"]][["n_groups"]])
  #GT_hyper
  expect_true(ng_ah_global[["prior_phi_cpp"]][["GT_hyper"]])
  expect_true(ng_ah_global_crela[["prior_phi_cpp"]][["GT_hyper"]])
  expect_true(ng_ah_semiglobal[["prior_phi_cpp"]][["GT_hyper"]])
  expect_true(ng_ah_semiglobal_c_rel_a[["prior_phi_cpp"]][["GT_hyper"]])
  expect_false(ng_afix_global[["prior_phi_cpp"]][["GT_hyper"]])
  expect_false(ng_afix_global_crela[["prior_phi_cpp"]][["GT_hyper"]])
  expect_false(ng_afix_semiglobal[["prior_phi_cpp"]][["GT_hyper"]])
  expect_false(ng_afix_semiglobal_crela[["prior_phi_cpp"]][["GT_hyper"]])
})

test_that("DL", {
  M <- 3L
  lags <- 2L
  a_fix1 <- 0.2
  a_fix2 <- "1/n"
  a_fix3 <- "1/K"
  sup <- seq(.1,1,.1)
  a_mat <- cbind(sup, dexp(sup))

  # a fixed, global
  dl <- specify_prior_phi(M = M, lags = lags, prior = "DL",
                          DL_a = a_fix1,
                          global_grouping = "global")
  dl2 <- specify_prior_phi(M = M, lags = lags, prior = "DL",
                          DL_a = a_fix2,
                          global_grouping = "global")
  dl3 <- specify_prior_phi(M = M, lags = lags, prior = "DL",
                          DL_a = a_fix3,
                          global_grouping = "global")
  # a fixed, semi-global
  dl_semiglobal <- specify_prior_phi(M = M, lags = lags, prior = "DL",
                          DL_a = a_fix1,
                          global_grouping = "olcl-lagwise")
  dl2_semiglobal <- specify_prior_phi(M = M, lags = lags, prior = "DL",
                                     DL_a = a_fix2,
                                     global_grouping = "olcl-lagwise")
  dl3_semiglobal <- specify_prior_phi(M = M, lags = lags, prior = "DL",
                                     DL_a = a_fix3,
                                     global_grouping = "olcl-lagwise")
  # a random, global
  dl_a <- specify_prior_phi(M = M, lags = lags, prior = "DL",
                          DL_a = a_mat,
                          global_grouping = "global")
  # a random, semi-global
  dl_a_semi <- specify_prior_phi(M = M, lags = lags, prior = "DL",
                          DL_a = a_mat,
                          global_grouping = "olcl-lagwise")

  # a == a_fix, and n_groups (implicitly)
  K <- M*lags
  n <- K*M
  expect_equal(dl$prior_phi_cpp$a, a_fix1)
  expect_equal(dl2$prior_phi_cpp$a, 1/n)
  expect_equal(dl3$prior_phi_cpp$a, 1/K)
  expect_equal(dl_semiglobal$prior_phi_cpp$a, rep(a_fix1, dl_semiglobal$prior_phi_cpp$n_groups))
  expect_equal(dl2_semiglobal$prior_phi_cpp$a, rep(1/n, dl2_semiglobal$prior_phi_cpp$n_groups))
  expect_equal(dl3_semiglobal$prior_phi_cpp$a, rep(1/K, dl3_semiglobal$prior_phi_cpp$n_groups))
  # check support points and weights
  expect_equal(dl_a$prior_phi_cpp$a_vec, a_mat[,1])
  expect_equal(dl_a$prior_phi_cpp$a_weight, a_mat[,2])
  expect_equal(dl_a_semi$prior_phi_cpp$a_vec, a_mat[,1])
  expect_equal(dl_a_semi$prior_phi_cpp$a_weight, a_mat[,2])
  # check n_groups for a random
  expect_equal(dl_a$prior_phi_cpp$n_groups, length(dl_a$prior_phi_cpp$a))
  expect_equal(dl_a_semi$prior_phi_cpp$n_groups, length(dl_a_semi$prior_phi_cpp$a))
  #DL_hyper
  expect_true(dl_a$prior_phi_cpp$DL_hyper)
  expect_true(dl_a_semi$prior_phi_cpp$DL_hyper)
  expect_false(dl$prior_phi_cpp$DL_hyper)
  expect_false(dl2$prior_phi_cpp$DL_hyper)
  expect_false(dl3$prior_phi_cpp$DL_hyper)
  expect_false(dl2_semiglobal$prior_phi_cpp$DL_hyper)
  expect_false(dl3_semiglobal$prior_phi_cpp$DL_hyper)
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

test_that("factor", {
  M <- 10L
  factors <- 3L
  ngcol <- specify_prior_sigma(M = M, type = "factor", factor_factors = factors,
                               factor_priorfacloadtype = "colwiseng",
                               quiet = TRUE)
  expect_true(ngcol$prior_sigma_cpp$factor_columnwise)
  expect_true(ngcol$prior_sigma_cpp$factor_ngprior)
  expect_true(all(ngcol$prior_sigma_cpp$factor_restrinv==1L))
  expect_equal(ngcol$prior_sigma_cpp$factor_factors, factors)

  ngrow_upper <- specify_prior_sigma(M = M, type = "factor", factor_factors = factors,
                                     factor_priorfacloadtype = "rowwiseng",
                                     factor_restrict = "upper", quiet = TRUE)
  upper_ind <- upper.tri(ngrow_upper$prior_sigma_cpp$factor_restrinv)
  expect_true(all(ngrow_upper$prior_sigma_cpp$factor_restrinv[upper_ind])==0L)
  expect_true(all(ngrow_upper$prior_sigma_cpp$factor_restrinv[!upper_ind])==1L)

  heteroscedastic <- specify_prior_sigma(M = M, type = "factor", factor_factors = factors,
                                         factor_heteroskedastic = TRUE,
                                         quiet = TRUE)
  expect_identical(length(heteroscedastic$prior_sigma_cpp$sv_heteroscedastic), M+factors)
  expect_true(all(heteroscedastic$prior_sigma_cpp$sv_heteroscedastic==TRUE))

  homoscedastic <- specify_prior_sigma(M = M, type = "factor", factor_factors = factors,
                                         factor_heteroskedastic = FALSE,
                                       quiet = TRUE)
  expect_identical(length(homoscedastic$prior_sigma_cpp$sv_heteroscedastic), M+factors)
  expect_true(all(homoscedastic$prior_sigma_cpp$sv_heteroscedastic==FALSE))

  idisv_fachomo <- specify_prior_sigma(M = M, type = "factor", factor_factors = factors,
                                      factor_heteroskedastic = c(TRUE,FALSE),
                                      quiet = TRUE)
  expect_identical(length(idisv_fachomo$prior_sigma_cpp$sv_heteroscedastic), M+factors)
  expect_true(all(idisv_fachomo$prior_sigma_cpp$sv_heteroscedastic[-c(1:M)]==FALSE))
  expect_true(all(idisv_fachomo$prior_sigma_cpp$sv_heteroscedastic[c(1:M)]==TRUE))

  idihomo_facsv <- specify_prior_sigma(M = M, type = "factor", factor_factors = factors,
                                       factor_heteroskedastic = c(FALSE, TRUE),
                                       quiet = TRUE)
  expect_identical(length(idihomo_facsv$prior_sigma_cpp$sv_heteroscedastic), M+factors)
  expect_true(all(idihomo_facsv$prior_sigma_cpp$sv_heteroscedastic[-c(1:M)]==TRUE))
  expect_true(all(idihomo_facsv$prior_sigma_cpp$sv_heteroscedastic[c(1:M)]==FALSE))

  heter <- rep_len(c(TRUE,FALSE), M + factors)
  svmix <- specify_prior_sigma(M = M, type = "factor", factor_factors = factors,
                                  factor_heteroskedastic = heter,
                               quiet = TRUE)
  expect_equal(svmix$prior_sigma_cpp$sv_heteroscedastic, heter)
})

test_that("cholesky", {
  M <- 3L
  a_fix <- 0.2
  sup <- seq(.1,1,.1)
  a_mat <- cbind(sup, dexp(sup))
  b_fix <- 0.4
  c_fix <- 0.8
  crela <- "0.3*a"

  # check hetero/homoscedasticity settings
  heter <- specify_prior_sigma(M = M, type = "cholesky", cholesky_heteroscedastic = TRUE,
                               quiet = TRUE)
  expect_equal(heter$prior_sigma_cpp$sv_heteroscedastic, rep(TRUE, M))
  homo <- specify_prior_sigma(M = M, type = "cholesky", cholesky_heteroscedastic = FALSE,
                               quiet = TRUE)
  expect_equal(homo$prior_sigma_cpp$sv_heteroscedastic, rep(FALSE, M))
  expect_identical(dim(homo$prior_sigma_cpp$cholesky_priorhomoscedastic), c(M, 2L))

  # DL
  dl <- specify_prior_sigma(M = M, type = "cholesky", cholesky_U_prior = "DL",
                            cholesky_DL_a = a_fix,
                            quiet = TRUE)
  expect_equal(dl$prior_sigma_cpp$cholesky_a, a_fix)
  expect_false(dl$prior_sigma_cpp$cholesky_DL_hyper)

  dl_a <- specify_prior_sigma(M = M, type = "cholesky", cholesky_U_prior = "DL",
                              cholesky_DL_a = a_mat,
                              quiet = TRUE)
  expect_equal(dl_a$prior_sigma_cpp$cholesky_a_vec, a_mat[,1])
  expect_equal(dl_a$prior_sigma_cpp$cholesky_a_weight, a_mat[,2])
  expect_true(dl_a$prior_sigma_cpp$cholesky_DL_hyper)

  # R2D2
  r2d2 <- specify_prior_sigma(M = M, type = "cholesky", cholesky_U_prior = "R2D2",
                              cholesky_R2D2_a = a_fix, cholesky_R2D2_b = b_fix,
                              quiet = TRUE)
  expect_equal(r2d2$prior_sigma_cpp$cholesky_a, a_fix)
  expect_equal(r2d2$prior_sigma_cpp$cholesky_b, b_fix)
  expect_equal(r2d2$prior_sigma_cpp$cholesky_c, 0.5*a_fix)
  expect_false(r2d2$prior_sigma_cpp$cholesky_GT_hyper)
  expect_true(r2d2$prior_sigma_cpp$cholesky_c_rel_a)
  expect_equal(r2d2$prior_sigma_cpp$cholesky_GT_priorkernel, "exponential")
  expect_equal(r2d2$prior_sigma_cpp$cholesky_GT_vs, 0.5)

  r2d2_a <- specify_prior_sigma(M = M, type = "cholesky", cholesky_U_prior = "R2D2",
                              cholesky_R2D2_a = a_mat, cholesky_R2D2_b = b_fix,
                              quiet = TRUE)
  expect_equal(r2d2_a$prior_sigma_cpp$cholesky_b, b_fix)
  expect_equal(r2d2_a$prior_sigma_cpp$cholesky_c, 0.5*r2d2_a$prior_sigma_cpp$cholesky_a)
  expect_equal(r2d2_a$prior_sigma_cpp$cholesky_c_vec, 0.5*r2d2_a$prior_sigma_cpp$cholesky_a_vec)
  expect_true(r2d2_a$prior_sigma_cpp$cholesky_GT_hyper)
  expect_true(r2d2_a$prior_sigma_cpp$cholesky_c_rel_a)
  expect_equal(r2d2_a$prior_sigma_cpp$cholesky_GT_priorkernel, "exponential")
  expect_equal(r2d2_a$prior_sigma_cpp$cholesky_GT_vs, 0.5)

  # NG
  ng <- specify_prior_sigma(M = M, type = "cholesky", cholesky_U_prior = "NG",
                            cholesky_NG_a = a_fix, cholesky_NG_b = b_fix,
                            cholesky_NG_c = c_fix,
                            quiet = TRUE)
  expect_equal(ng$prior_sigma_cpp$cholesky_a, a_fix)
  expect_equal(ng$prior_sigma_cpp$cholesky_b, b_fix)
  expect_equal(ng$prior_sigma_cpp$cholesky_c, c_fix)
  expect_false(ng$prior_sigma_cpp$cholesky_GT_hyper)
  expect_false(ng$prior_sigma_cpp$cholesky_c_rel_a)
  expect_equal(ng$prior_sigma_cpp$cholesky_GT_priorkernel, "normal")
  expect_equal(ng$prior_sigma_cpp$cholesky_GT_vs, 1)

  ng_a <- specify_prior_sigma(M = M, type = "cholesky", cholesky_U_prior = "NG",
                              cholesky_NG_a = a_mat, cholesky_NG_b = b_fix,
                              cholesky_NG_c = c_fix,
                              quiet = TRUE)
  expect_equal(ng_a$prior_sigma_cpp$cholesky_a_vec, a_mat[,1])
  expect_equal(ng_a$prior_sigma_cpp$cholesky_a_weight, a_mat[,2])
  expect_equal(ng_a$prior_sigma_cpp$cholesky_b, b_fix)
  expect_equal(ng_a$prior_sigma_cpp$cholesky_c, c_fix)
  expect_true(ng_a$prior_sigma_cpp$cholesky_GT_hyper)
  expect_false(ng_a$prior_sigma_cpp$cholesky_c_rel_a)
  expect_equal(ng_a$prior_sigma_cpp$cholesky_GT_priorkernel, "normal")
  expect_equal(ng_a$prior_sigma_cpp$cholesky_GT_vs, 1)

  ng_a_c <- specify_prior_sigma(M = M, type = "cholesky", cholesky_U_prior = "NG",
                              cholesky_NG_a = a_mat, cholesky_NG_b = b_fix,
                              cholesky_NG_c = crela,
                              quiet = TRUE)
  q <- as.numeric(strsplit(crela, "\\*")[[1]][1])
  expect_equal(ng_a_c$prior_sigma_cpp$cholesky_a_vec, a_mat[,1])
  expect_equal(ng_a_c$prior_sigma_cpp$cholesky_a_weight, a_mat[,2])
  expect_equal(ng_a_c$prior_sigma_cpp$cholesky_b, b_fix)
  expect_equal(ng_a_c$prior_sigma_cpp$cholesky_c, q*ng_a_c$prior_sigma_cpp$cholesky_a)
  expect_equal(ng_a_c$prior_sigma_cpp$cholesky_c_vec, q*ng_a_c$prior_sigma_cpp$cholesky_a_vec)
  expect_true(ng_a_c$prior_sigma_cpp$cholesky_GT_hyper)
  expect_true(ng_a_c$prior_sigma_cpp$cholesky_c_rel_a)
  expect_equal(ng_a_c$prior_sigma_cpp$cholesky_GT_priorkernel, "normal")
  expect_equal(ng_a_c$prior_sigma_cpp$cholesky_GT_vs, 1)
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
