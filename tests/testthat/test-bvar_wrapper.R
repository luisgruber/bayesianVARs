data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]

test_that("stable_bvar", {
  for(type in c("cholesky", "factor")){
    prior_sigma <- specify_prior_sigma(data, type = type, quiet = TRUE)
    mod <- bvar(data, sv_keep = "all", draws = 10L, burnin = 2L, 
                prior_sigma = prior_sigma, quiet = TRUE)
    expect_no_condition(stable_bvar(mod))
    expect_no_condition(stable_bvar(mod, quiet = TRUE))
    expect_no_condition(stable_bvar(mod, lim = 0.99))
    expect_error(stable_bvar(mod, lim = 1.1))
    expect_error(stable_bvar(mod, lim = -1)) 
  }
})
