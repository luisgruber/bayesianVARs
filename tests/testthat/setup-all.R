data("exrates", package = "stochvol")

y <- factorstochvol::logret(exrates[1:15, 1:3], demean = TRUE)
draws <- 27L
burnin <- 10L

# Mainly PHI specifications -----------------------------------------------

thin_values <- c(1,3)
lag_values <- c(0L,1L,2L,3L)
PHI_priors <- c("normal","HMP", "SSVS", "NG", "DL", "R2D2", "HS") #"HMP"
intercept_values <- c(0,10)
sv_keep <- c("all","last")
semi_groups <- c("global", "olcl-lagwise")
keep <- "last"
error_spec <- c("cholesky", "factor") 
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
          if(prior == "normal"){
            prior_phi <- specify_prior_phi(data = y, lags = lags, prior = prior, normal_sds = 10)
            res <- c(res, list(bvar(y , lags = lags, prior_intercept = myintercept,
                                    draws = draws, burnin = burnin, thin = thin,
                                    prior_phi = prior_phi, prior_sigma = prior_sigma,
                                    sv_keep = keep,
                                    quiet = TRUE)))
          }else if(prior == "HMP"){
            prior_phi <- specify_prior_phi(data = y, lags = lags, prior = prior)
            res <- c(res, list(bvar(y , lags = lags, prior_intercept = myintercept,
                                    draws = draws, burnin = burnin, thin = thin,
                                    prior_phi = prior_phi, prior_sigma = prior_sigma,
                                    sv_keep = keep,
                                    quiet = TRUE)))
          }else if(prior == "SSVS" || prior == "NG" || prior == "DL" ||
                   prior == "R2D2" || prior == "HS"){
            for(group in semi_groups){
              if(prior == "SSVS"){
                for(ssvs_p in c(0.5, 1)){
                  if(ssvs_p == 1){
                    ssvs_p = c(1,1)
                  }
                  prior_phi <- specify_prior_phi(data = y, lags = lags, prior = prior,
                                                 global_grouping = group,
                                                 SSVS_p = ssvs_p)
                  res <- c(res, list(bvar(y , lags = lags, prior_intercept = myintercept,
                                          draws = draws, burnin = burnin, thin = thin,
                                          prior_phi = prior_phi, prior_sigma = prior_sigma,
                                          sv_keep = keep,
                                          quiet = TRUE)))
                }
              }else if(prior == "NG" || prior == "DL" || prior == "R2D2" ){
                for(a_shrink in c(0.5,1)){
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
                  res <- c(res, list(bvar(y , lags = lags, prior_intercept = myintercept,
                                          draws = draws, burnin = burnin, thin = thin,
                                          prior_phi = prior_phi, prior_sigma = prior_sigma,
                                          sv_keep = keep,
                                          quiet = TRUE)))
                }
                
              }else if(prior == "HS"){
                prior_phi <- specify_prior_phi(data = y, lags = lags, prior = prior,
                                               global_grouping = group)
                res <- c(res, list(bvar(y , lags = lags, prior_intercept = myintercept,
                                        draws = draws, burnin = burnin, thin = thin,
                                        prior_phi = prior_phi, prior_sigma = prior_sigma,
                                        sv_keep = keep,
                                        quiet = TRUE)))
              }
            }
          }
        }
      }
    }
  }
  #}
}


# SIGMA cholesky ----------------------------------------------------------

SIGMA_priors <- c("normal","HMP", "SSVS", "NG", "DL", "R2D2", "HS")
heteroscedastic <- c(TRUE,FALSE)
for(thin in thin_values){
  for(scedastic in heteroscedastic){
    for(prior in SIGMA_priors){

      if(prior == "normal" |prior == "HMP" | prior =="HS"){
        prior_sigma <- specify_prior_sigma(data = y, type = "cholesky",
                                           cholesky_U_prior = prior,
                                           cholesky_heteroscedastic = scedastic,
                                           quiet = TRUE)
        res <- c(res, list(bvar(y , draws = draws, burnin = burnin, thin = thin,
                                prior_sigma = prior_sigma,
                                sv_keep = keep,
                                quiet = TRUE)))
      }else if(prior == "NG" | prior == "R2D2" | prior == "DL"){
        for(a_shrink in c(0.5,1)){
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
        }
        res <- c(res, list(bvar(y ,draws = draws, burnin = burnin, thin = thin,
                                prior_sigma = prior_sigma,
                                sv_keep = keep,
                                quiet = TRUE)))
      }else if(prior == "SSVS" ){
        for(ssvs_p in c(0.5, 1)){
          if(ssvs_p == 1){
            ssvs_p = c(1,1)
          }
          prior_sigma <- specify_prior_sigma(data = y, type = "cholesky",
                                           cholesky_U_prior = prior,
                                           cholesky_heteroscedastic = scedastic,
                                           cholesky_SSVS_p = ssvs_p,
                                           quiet = TRUE)
          res <- c(res, list(bvar(y , draws = draws, burnin = burnin, thin = thin,
                                  prior_sigma = prior_sigma,
                                  sv_keep = keep,
                                  quiet = TRUE)))
        }

      }
    }
  }
}


# SIGMA factor ------------------------------------------------------------

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
          res <- bvar(y,draws = draws, burnin = burnin, prior_sigma = prior_sigma ,quiet=TRUE,
                      expert_huge = ha)
        }
      }
    } 
  }
}
