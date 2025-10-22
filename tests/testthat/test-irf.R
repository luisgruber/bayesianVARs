train_data <- 100 * usmacro_growth[,c("GDPC1", "GDPCTPI", "GS1", "M2REAL", "CPIAUCSL")]
choices_sigma_type <- list(
	cholesky = list(
		prior=specify_prior_sigma(M=ncol(train_data), type="cholesky", cholesky_heteroscedastic=FALSE),
		n_shocks = ncol(train_data),
		restrictions_B0=rbind(
			c(1 ,NA,0 ,NA,NA),
			c(0 ,1 ,0 ,NA,NA),
			c(0 ,NA,1 ,NA,NA),
			c(0 ,0 ,NA,1 ,NA),
			c(0 ,0 ,0 ,0 ,1 )
  		),
  		restrictions_facload = NULL
	),
	factor = list(
		prior=specify_prior_sigma(M=ncol(train_data), type="factor", factor_factors=2L, factor_heteroskedastic=FALSE),
		n_shocks = 2,
		restrictions_B0 = NULL,
  		restrictions_facload = rbind(
			c(1,NA),
			c(0,1),
			c(NA,NA),
			c(NA,NA),
			c(NA,NA)
  		)
	)
)

for (sigma_type in choices_sigma_type) {
	mod <- bvar(train_data, lags=5L, burnin=0, draws=12, prior_sigma=sigma_type$prior, quiet=TRUE)
	for (do_include_restrictions in c(FALSE, TRUE))
	for (hairy in c(FALSE, TRUE))
	for (solver in c("randomized", "lp"))
		test_that(
			paste("irf completes with",
				"sigma type:",sigma_type$prior$type,
				"restr:",do_include_restrictions,
				"hair:", hairy,
				"solver", solver
			),
		{
			structural_restrictions <- NULL
			if (do_include_restrictions) {
				structural_restrictions <- specify_structural_restrictions(
					mod,
					restrictions_B0 = sigma_type$restrictions_B0,
					restrictions_facload = sigma_type$restrictions_facload
				)
			}
			result <- irf(
				mod, ahead=8,
				structural_restrictions=structural_restrictions,
				hairy=hairy,
				solver=solver
			)
	
			# do not test the sample counte, because some samples can be rejected
			expect_equal(dim(result)[-4], c(ncol(train_data), sigma_type$n_shocks, 1+8))
			if (hairy) {
				expect_length(attr(result, "hair_order"), dim(result)[4])
			} else {
				expect_null(attr(result, "hair_order"))
			}
			if (solver == "lp") {
				expect_equal(dim(result)[4], mod$config$draws)
			}
		}
	)
}
