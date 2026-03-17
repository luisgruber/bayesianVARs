#' @useDynLib bayesianVARs, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom GIGrvg rgig
#' @importFrom stats rnorm dnorm rgamma embed cor ts.plot fitted quantile sd
#' @importFrom utils flush.console
#' @importFrom MASS mvrnorm
#' @importFrom stochvol svsample
#' @importFrom factorstochvol fsvsim
#' @importFrom colorspace sequential_hcl diverge_hcl
#' @importFrom graphics abline axis par rect text lines mtext pairs points polygon strheight strwidth
#' @importFrom scales alpha
NULL
