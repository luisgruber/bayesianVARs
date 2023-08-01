#' Plot method for BVAR coefficients
#'
#' @param object \code{PHI}-object obtained from \code{\link{bvar}}.
#' @param summary character indicating the posterior summary to be visualized.
#'   One of \code{"median"}, \code{"mean"}, \code{"IQR"}, \code{"sd"} or
#'   \code{"var"}.
#' @param ylabels \code{ylabels=NULL}, the default, indicates that the names of
#'   the dependent variables will be displayed. \code{ylabels=""} indicates that
#'   no ylabels will be displayed.
#' @param xlabels \code{xlabels=NULL}, the default, indicates that the labels of
#'   all covariables (the lagged values of the dependent variables) will be
#'   displayed. \code{xlabels="lags"} indicates that only the lags will be
#'   marked. \code{xlabels=""} indicates that no ylabels are displayed.
#' @param add_numbers logical. \code{add_numbers=TRUE}, the default indicates
#'   that the actual values of \code{summary} will be displayed.
#' @param zlim numeric vector of length two indicating the minimum and maximum
#'   values for which colors should be plotted. By default this range is
#'   determined by the maximum of the absolute values of the selected summary.
#' @param ...
#'
#' @export
#'
plot.PHI <- function(object, summary = "median", colorbar = TRUE, ylabels = NULL, xlabels = NULL,
                     add_numbers = FALSE, zlim = NULL,main="",...){

  optionals <- list(...)

  zlim_in <- zlim
  if(add_numbers){
    alpha <- .5
  }else{
    alpha <- 1
  }

  if(colorbar | length(summary)>1){
    nrplots <- length(summary)
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    if(colorbar){
      #mat <- matrix(sort(c(rep(seq(1,3*nrplots,3),9), seq(2,3*nrplots,3),seq(3,3*nrplots,3))), nrow = 11*nrplots)
      mat <- matrix(sort(c(rep(seq(1,2*nrplots,2),9), seq(2,2*nrplots,2))), nrow = 10*nrplots)
      layout(mat)
    }else{
      par(mfrow=c(nrplots,1))
    }

  }

  for(i in seq.int(length(summary))){
    PHI <- apply(object, 2:3, function(x) do.call(what = summary[i],
                                                  args = list(x)))
    PHI_star <- t(apply(PHI, 1, rev)) # image orders differently

    if(summary[i] %in% c("median", "mean")){
      colspace <- colorspace::diverge_hcl(1001, alpha = alpha, palette = "Blue-Red")
      if(is.null(zlim_in)){
        zlim <- c(-max(abs(PHI)),max(abs(PHI)))
      }
      colbreaks <- seq(zlim[1]*1.001, zlim[2]*1.001, len=1002)#[-1]
    }else if(summary[i] %in% c("sd", "var", "IQR")){
      colspace <- colorspace::sequential_hcl(1001, alpha = alpha, rev = TRUE,
                                             palette = "Reds 2")#colorspace::sequential_hcl(1000, alpha = alpha, rev = TRUE)#colorspace::sequential_hcl(5, h = 0, c = c(100, 0), l = 65, rev = TRUE, power = 1, alpha = alpha) #colorspace::sequential_hcl(1000, alpha = alpha, rev = TRUE)
      if(is.null(zlim_in)){
        zlim <- c(0,max(abs(PHI)))
      }
    }

    M <- ncol(PHI)
    Kp <- nrow(PHI)
    p <- floor(Kp/M)

    if(colorbar){
      #      oldpar <- par(no.readonly = TRUE)
      #      on.exit(par(oldpar), add = TRUE)
      #      mat <- matrix(c(rep(1,9),2),nrow = 10)
      #      layout(mat)
      par(mar=c(0,4,4,2)+.1) ###
    }

    image((PHI_star), zlim = zlim ,xaxt = "n", yaxt="n", col = colspace,
          bty="n", main = main)
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

    if(colorbar){
      #    mymar2 <- oldpar$mar
      #    mymar2[3] <- .1
      #    par(mar=mymar2)
      par(mar=c(2,4,0,2)+.1)
      plot.new()
      if(summary[i] %in% c("median", "mean")){

        colbarlim <- length(which(colbreaks<min(PHI_star)))
        #colbarlim <- ceiling((min(PHI_star) - zlim[1]*1.001)/(diff(zlim)+2*zlim[2]/1000)*1001)
        if(colbarlim>1){
          adjustedcolspace <- colspace[-c(1:(colbarlim-1))]
        }else adjustedcolspace <- colspace

        if(sign(min(PHI_star))!=sign(max(PHI_star))){
          llabels <- round(c(colbreaks[colbarlim],0,tail(colbreaks,1)),2)
          #llabels <- round(sort(c(range(PHI_star),0)),2)
          #
          at <- sort(c(0,1,(0-colbreaks[colbarlim])/(tail(colbreaks,1)-colbreaks[colbarlim])))
        }else{
          llabels <- round(range(PHI_star),2)
          at <- c(0,1)
        }
      }else if(summary[i] %in% c("sd", "var", "IQR")){

        adjustedcolspace <- colspace
        llabels <- round(zlim,2)
        at <- c(0,1)

      }

      len <- length(adjustedcolspace)
      xxx <- seq(0,1,length=len+1)
      rect(xxx[1:len], rep(0, len), xxx[-1],
           rep(1, len), col = adjustedcolspace, border = adjustedcolspace)

      axis(1, labels = llabels, at = at,
           tick = FALSE, las=1, cex.axis=1.25, line=-.8)
      #plot.new()
    }

  }

}

#' Plot method for contemporaneous coefficients
#'
#' @param object \code{L}-object obtained from \code{\link{bvar}}.
#' @param summary character indicating the posterior summary to be visualized.
#'   One of \code{"median"}, \code{"mean"}, \code{"IQR"}, \code{"sd"} or
#'   \code{"var"}.
#' @param ylabels \code{NULL}, the default, indicates that ylabels will be
#'   displayed. \code{""} indicates that no ylabels will be displayed.
#' @param xlabels \code{NULL}, the default, indicates that xlabels will be
#'   displayed. \code{""} indicates that no xlabels will be displayed.
#' @param add_numbers logical. \code{add_numbers=TRUE}, the default indicates
#'   that the actual values of \code{summary} will be displayed.
#' @param zlim numeric vector of length two indicating the minimum and maximum
#'   values for which colors should be plotted. By default this range is
#'   determined by the maximum of the absolute values of the selected summary.
#' @param ...
#'
#' @export
#'
plot.L <- function(object, summary = "median", ylabels = NULL, xlabels = NULL,
                     add_numbers = TRUE, zlim = NULL, main="", ...){

  L <- apply(object, 2:3, function(x) do.call(what = summary,
                                                args = list(x)))
  M <- ncol(L)
  L[lower.tri(L, diag = TRUE)] <- NA
  maxz <- max(abs(L), na.rm = TRUE)

  L_star <- t(apply(L, 1, rev)) # image orders differently

  if(add_numbers){
    alpha <- .5
  }else{
    alpha <- 1
  }

  if(summary %in% c("median", "mean")){
    colspace <- colorspace::diverge_hcl(1000, alpha = alpha, palette = "Blue-Red")
    if(is.null(zlim)){
      zlim <- c(-maxz,maxz)
    }
  }else if(summary %in% c("sd", "var", "IQR")){
    colspace <- colorspace::sequential_hcl(1000, alpha = alpha, rev = TRUE,
                                           palette = "Reds 2")#colorspace::sequential_hcl(1000, alpha = alpha, rev = TRUE)#colorspace::sequential_hcl(5, h = 0, c = c(100, 0), l = 65, rev = TRUE, power = 1, alpha = alpha) #colorspace::sequential_hcl(1000, alpha = alpha, rev = TRUE)
    if(is.null(zlim)){
      zlim <- c(0,maxz)
    }
  }

  image((L_star), zlim= zlim ,xaxt = "n", yaxt="n", col = colspace, bty="n", main=main)
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

