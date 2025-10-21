data("exrates", package = "stochvol")

y <- factorstochvol::logret(exrates[1:15, 1:3], demean = TRUE)
draws <- 27
burnin <- 10

# Mainly PHI specifications -----------------------------------------------

thin_values <- c(1,3)
lag_values <- c(1L,2L,3L)
PHI_priors <- c("normal","HMP", "SSVS", "NG", "DL", "R2D2", "HS") #"HMP"
intercept_values <- c(0,10)
sv_keep <- c("all","last")
semi_groups <- c("global", "olcl-lagwise")
keep <- "last"
# test prior_phi options
res <- list(bvar(y,draws = draws, burnin = burnin,quiet=TRUE))
for(thin in thin_values){
  #for(keep in sv_keep){
    for(lags in lag_values){
      for(intercept in intercept_values){
        myintercept <- if(intercept==0) FALSE else intercept
        for(prior in PHI_priors){
          if(prior == "normal"){
            prior_phi <- specify_prior_phi(data = y, lags = lags, prior = prior, normal_sds = 10)
            res <- c(res, list(bvar(y , lags = lags, prior_intercept = myintercept,
                                    draws = draws, burnin = burnin, thin = thin,
                                    prior_phi = prior_phi,
                                    sv_keep = keep,
                                    quiet = TRUE)))
          }else if(prior == "HMP"){
            prior_phi <- specify_prior_phi(data = y, lags = lags, prior = prior)
            res <- c(res, list(bvar(y , lags = lags, prior_intercept = myintercept,
                                    draws = draws, burnin = burnin, thin = thin,
                                    prior_phi = prior_phi,
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
                                          prior_phi = prior_phi,
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
                                          prior_phi = prior_phi,
                                          sv_keep = keep,
                                          quiet = TRUE)))
                }

              }else if(prior == "HS"){
                prior_phi <- specify_prior_phi(data = y, lags = lags, prior = prior,
                                             global_grouping = group)
                res <- c(res, list(bvar(y , lags = lags, prior_intercept = myintercept,
                                        draws = draws, burnin = burnin, thin = thin,
                                        prior_phi = prior_phi,
                                        sv_keep = keep,
                                        quiet = TRUE)))
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
res <- list(bvar(y,draws = draws, burnin = burnin,quiet = TRUE))
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

# factors_values <- c(0, 1, 3)
# restrict_mat <- matrix(FALSE, nrow = NCOL(y), ncol = max(factors_values))
# restrict_mat[1, 1] <- TRUE  # restrict the upper left element to zero
# restrict_values <- list("upper", "none")
# priorfacloadtype_values <- c("normal", "rowwiseng", "colwiseng")
# priorhomoskedastic <- matrix(c(1.1, 1.1), nrow = NCOL(y),
#                              ncol = 2, byrow = TRUE)
# heteroskedastic_values <- list(TRUE, c(FALSE, FALSE))
#
#
# for (th in thin_values) {
#   for (pflt in priorfacloadtype_values) {
#     for (hsk in heteroskedastic_values) {
#       for (fs in factors_values) {
#         res <- c(res, list(fsvsample(y, draws = draws, burnin = burnin, quiet = TRUE,
#                                      factors = fs, thin = th, priorfacloadtype = pflt,
#                                      restrict = "none",
#                                      heteroskedastic = hsk,
#                                      priorhomoskedastic = if (!isTRUE(hsk)) priorhomoskedastic else NA,
#                                      interweaving = if (!isTRUE(hsk)) 0 else 4,
#                                      runningstore = if (fs == 0) 1 else 6)))
#         if (fs > 1) {
#           for (rst in restrict_values) {
#             res <- c(res, list(fsvsample(y, draws = draws, burnin = burnin, quiet = TRUE,
#                                          factors = fs, thin = th, priorfacloadtype = pflt,
#                                          heteroskedastic = hsk,
#                                          priorhomoskedastic = if (!isTRUE(hsk)) priorhomoskedastic else NA,
#                                          interweaving = if (!isTRUE(hsk)) 0 else 4,
#                                          restrict = rst)))
#           }
#         }
#       }
#     }
#   }
# }
