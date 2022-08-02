#' Plot method for BVAR coefficients
#'
#' @param object \code{PHI}-object obtained from \code{\link{bvar}}.
#' @param summary character indicating the posterior summary to be visualized.
#'   One of \code{"median"}, \code{"mean"}, \code{"IQR"}, \code{"sd"} or
#'   \code{"var"}.
#' @param ylabels \code{ylabels=""} indicates that no ylabels are displayed.
#'   \code{ylabels=NULL}, the default, indicates that the names of the dependent
#'   variables will be displayed.
#' @param xlabels \code{xlabels=NULL}, the default, indicates that the lags all
#'   dependent variables will be displayed. \code{xlabels="lags"} indicates that
#'   only the lags will be displayed. \code{xlabels=""} indicates that no
#'   ylabels are displayed.
#' @param add_numbers logical. \code{add_numbers=TRUE}, the default indicates
#'   that the actual values of \code{summary} will be displayed.
#' @param ...
#'
#' @export
#'
plot.PHI <- function(object, summary = "median", ylabels = NULL, xlabels = NULL,
                 add_numbers = TRUE, ...){

  PHI <- apply(object, 2:3, function(x) do.call(what = summary,
                                                 args = list(x)))
  PHI_star <- t(apply(PHI, 1, rev)) # image orders differently

  if(add_numbers){
    alpha <- .5
  }else{
    alpha <- 1
  }

  if(summary %in% c("median", "mean")){
    colspace <- colorspace::diverge_hcl(1000, alpha = alpha, palette = "Blue-Red")
    zlim <- c(-max(abs(PHI)),max(abs(PHI)))
  }else if(summary %in% c("sd", "var", "IQR")){
    colspace <- colorspace::sequential_hcl(1000, alpha = alpha, rev = TRUE,
                                           palette = "Reds 2")#colorspace::sequential_hcl(1000, alpha = alpha, rev = TRUE)#colorspace::sequential_hcl(5, h = 0, c = c(100, 0), l = 65, rev = TRUE, power = 1, alpha = alpha) #colorspace::sequential_hcl(1000, alpha = alpha, rev = TRUE)
    zlim <- c(0,max(abs(PHI)))
  }

  M <- ncol(PHI)
  Kp <- nrow(PHI)
  p <- floor(Kp/M)
  image((PHI_star), zlim = zlim ,xaxt = "n", yaxt="n", col = colspace,
        bty="n")
  if((Kp-1)>M){
    abline(v = seq(M,Kp-1, M)/(Kp-1) - 0.5/(Kp-1), xpd = FALSE)
  }
  if(is.null(xlabels)){
    axis(3, labels = rownames(PHI), at = seq(0,1, length.out=Kp),
         tick = FALSE, las=2, cex.axis=.75)
  }else if(xlabels=="lags"){
    text(seq(M,Kp, M)/(Kp-1) - 0.5/(Kp-1) - M/(2*(Kp-1)), 1 + 0.5/(M-1) ,
         pos = 3, xpd = TRUE, labels = paste0("Lag ", 1:p))
  }
  if(is.null(ylabels)){
    axis(2, at = seq(1,0, length.out=M), labels = colnames(PHI),
         tick = FALSE, las = 2, cex.axis=.75)
  }

  if(add_numbers){
    text(rep(seq(0,1, length.out=Kp),M),
         sort(rep(seq(1,0, length.out=M),Kp), decreasing = TRUE),
         round(as.vector(PHI),3), cex=.75)
  }

}

#'
#' @export
#'
plot.L <- function(object, summary = "median", ylabels = NULL, xlabels = NULL,
                     add_numbers = TRUE, ...){

  L <- apply(object, 2:3, function(x) do.call(what = summary,
                                                args = list(x)))
  M <- ncol(L)
    L[lower.tri(L, diag = TRUE)] <- NA

  L_star <- t(apply(L, 1, rev)) # image orders differently

  if(add_numbers){
    alpha <- .5
  }else{
    alpha <- 1
  }

  if(summary %in% c("median", "mean")){
    colspace <- colorspace::diverge_hcl(1000, alpha = alpha)
  }else if(summary %in% c("sd", "var", "IQR")){
    colspace <- colorspace::sequential_hcl(1000, alpha = alpha, rev = TRUE)
  }

  maxz <- max(abs(L), na.rm = TRUE)
  image((L_star), zlim=c(-maxz,maxz) ,xaxt = "n", yaxt="n", col = colspace, bty="n")
  if(is.null(xlabels)){
    axis(3, labels = rownames(L), at = seq(0,1, length.out=M),
         tick = FALSE, las=2, cex.axis=.75)
  }
  if(is.null(ylabels)){
    axis(2, at = seq(1,0, length.out=M), labels = colnames(L),
         tick = FALSE, las = 2, cex.axis=.75)
  }

  if(add_numbers){
    text(rep(seq(0,1, length.out=M),M),
         sort(rep(seq(1,0, length.out=M),M), decreasing = TRUE),
         round(as.vector(L),3), cex=.75)
  }

}

