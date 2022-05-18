#' @useDynLib bayesianVARs, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom GIGrvg rgig
#' @importFrom stats rnorm dnorm rgamma dgamma rbeta rbinom rmultinom sd embed
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom MASS mvrnorm
#' @importFrom stochvol svsample validate_and_process_expert specify_priors svsample_fast_cpp
#' @importFrom colorspace sequential_hcl diverge_hcl
## specify priors and svsample_fast_cpp for bvar_R
NULL
