#' @export
bvar_fast <- function(y, draws, burnin, sv_spec){

  M <- ncol(y)
  T <- nrow(y)

  if(is.null(sv_spec)){
    sv_spec <- list(priormu = c(0,100),
                    priorphi = c(20, 1.5),
                    priorsigma2 = c(0.5,0.5)#,
                    #priorh0 = -1 #h0 from stationary distribution
                    )
  }

  sv_para_init <- matrix(data= c(rep(-10,M), rep(0.9,M), rep(0.2,M), rep(-10,M)),
                                nrow = 4, ncol = M, byrow = TRUE)
  rownames(sv_para_init) <- c("mu", "phi", "sigma", "h0")

  h_init <- matrix(rep(-10, T*M), T,M)

  res <- bvar_cpp(y,
                  draws,
                  burnin,
                  M,
                  T,
                  #expert,
                  sv_spec,
                  h_init,
                  sv_para_init
                  )
  res
}
