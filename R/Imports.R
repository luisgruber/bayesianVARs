#' @useDynLib bayesianVARs, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom GIGrvg rgig
#' @importFrom stats rnorm dnorm rgamma dgamma rbeta rbinom rmultinom sd embed cor fitted quantile ts.plot
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom MASS mvrnorm
#' @importFrom stochvol svsample validate_and_process_expert specify_priors svsample_fast_cpp
#' @importFrom factorstochvol fsvsim
#' @importFrom colorspace sequential_hcl diverge_hcl
#' @importFrom graphics abline axis image layout par plot.new rect text
#' @importFrom utils tail
#' @importFrom scales alpha
#' @importFrom graphics lines mtext pairs points polygon strheight strwidth
## specify priors and svsample_fast_cpp for bvar_R
NULL
