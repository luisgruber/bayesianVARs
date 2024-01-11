#' Posterior heatmaps for VAR coefficients or variance-covariance matrices
#'
#' @param x An array of dimension \eqn{a \times b \times draws}, where \eqn{a
#'   \times b} is the dimension of the parameter to visualize and draws is the
#'   number of posterior draws.
#' @param FUN The summary function to be applied to margins `c(1,2)` of x. E.g.
#'   \code{"median"}, \code{"mean"}, \code{"IQR"}, \code{"sd"} or \code{"var"}.
#'   `apply(x, 1:2, FUN, ...)` must return a matrix!
#' @param ... optional arguments to `FUN`.
#' @param colorbar logical indicating whether to display a colorbar or not.
#'   Default is \code{TRUE}.
#' @param xlabels \code{ylabels=NULL}, the default, indicates that the names of
#'   the dependent variables will be displayed. \code{ylabels=""} indicates that
#'   no ylabels will be displayed.
#' @param ylabels \code{xlabels=NULL}, the default, indicates that the labels of
#'   all covariables (the lagged values of the dependent variables) will be
#'   displayed. \code{xlabels="lags"} indicates that only the lags will be
#'   marked. \code{xlabels=""} indicates that no ylabels are displayed.
#' @param add_numbers logical. \code{add_numbers=TRUE}, the default indicates
#'   that the actual values of \code{summary} will be displayed.
#' @param zlim numeric vector of length two indicating the minimum and maximum
#'   values for which colors should be plotted. By default this range is
#'   determined by the maximum of the absolute values of the selected summary.
#' @param colspace Optional argument.
#' @param main main title for the plot.
#' @param cex.axis The magnification to be used for y-axis annotation relative
#'   to the current setting of cex.
#' @param cex.colbar The magnification to be used for colorbar annotation
#'   relative to the current setting of cex.
#' @param cex.numbers The magnification to be used for the actual values (if
#'   `add_numbers=TRUE`) relative to the current setting of cex.
#' @param asp aspect ratio.
#'
#' @return Returns `x` invisibly.
#' @seealso Other plotting [`plot.bayesianVARs_bvar()`],
#'   [`plot.bayesianVARs_fitted()`], [`plot.bayesianVARs_predict()`],
#'   [`pairs.bayesianVARs_predict()`].
#' @export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Estimate a model
#' mod <- bvar(100*data, sv_keep = "all", quiet = TRUE)
#'
#' # Extract posterior draws of VAR coefficients
#' phi_post <- coef(mod)
#'
#' # Visualize posterior median of VAR coefficients
#' posterior_heatmap(phi_post, median)
#'
#' # Extract posterior draws of variance-covariance matrices (for each point in time)
#' sigma_post <- vcov(mod)
#' # Visualize posterior interquartile-range of variance-covariance matrix of the first observation
#' posterior_heatmap(sigma_post[1,,,], IQR)
posterior_heatmap <- function(x,
                              FUN,
                              ...,
                              colorbar = TRUE,
                              xlabels = NULL,
                              ylabels = NULL,
                              add_numbers = FALSE,
                              zlim = NULL,
                              colspace = NULL,
                              main="",
                              cex.axis = 0.75,
                              cex.colbar = 1,
                              cex.numbers = 1,
                              asp=NULL){



  # PHI <- apply(x, 1:2, function(x) do.call(what = summary,
  #                                          args = list(x)))
  PHI <- apply(x, 1:2, FUN, ...)
  PHI_range <- range(PHI)
  if(!is.matrix(PHI)){
    stop("Invalid 'FUN'. apply(x,1:2,FUN,...) must return a matrix!")
  }

  if(base::isFALSE(add_numbers)){
    alpha <- 1
  }else{
    alpha <- 0.5
  }

  if(is.null(colspace)){
    if(PHI_range[1]<0){
      colspace <- colorspace::diverge_hcl(1001, alpha = alpha, palette = "Blue-Red")
    }else{
      colspace <- colorspace::sequential_hcl(1001, alpha = alpha, rev = TRUE,
                                             palette = "Reds 2")
    }
  }
  colspace_length <- length(colspace)
  if(is.null(zlim)){
    if(PHI_range[1]<0){
      zlim <- c(-max(abs(PHI_range)),max(abs(PHI_range)))
    }else{
      zlim <- c(0,max(abs(PHI_range)))
    }
  }

  M <- ncol(PHI)
  K <- nrow(PHI)
  p <- floor(K/M)

  row_names <- rownames(PHI)
  col_names <- colnames(PHI)
  maxstrwidth_left <- max(strwidth(row_names, "inches", cex = cex.axis))
  maxstrheight_top <- max(strwidth(col_names, "inches", cex = cex.axis)) # because las=2

  colbar_labels <- if(PHI_range[1]<0 & PHI_range[2]>0) c(PHI_range[1], 0, PHI_range[2]) else PHI_range
  colbar_labels <- round(colbar_labels, 2)
  maxstrwidth_right <- max(strwidth(colbar_labels, units = "inches", cex = cex.colbar))
  strheight_right <- strheight(c(colbar_labels[1], colbar_labels[length(colbar_labels)]), units = "inches", cex = cex.colbar)

  cols_to_plot_indices <- round(((as.vector(PHI)-zlim[1])/diff(zlim))*(colspace_length-1)+1,0)
  colbar_range <- colspace[min(cols_to_plot_indices):max(cols_to_plot_indices)]
  cols_to_plot <- colspace[cols_to_plot_indices]

  # width plus space for colorbar
  width <- 30
  xlim_colorbar <- c(width-2, width)
  xlim_whitespace <- c(xlim_colorbar[1]-1,xlim_colorbar[1])
  xlim_heatmap <- c(0,xlim_whitespace[1])
  x_breaks <- seq(xlim_heatmap[1] + diff(xlim_heatmap)/M,xlim_heatmap[2], length.out = M)

  # height
  height <- K/M*diff(xlim_heatmap)
  y_breaks <- rev(seq(0+height/K, height, length.out = K))

  #par()$mar/par()$mai=5; convert inches to lines
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mgp=c(2,.1,0))
  mar <- rep(0.1,4)
  par(mar=mar)
  if(colorbar){
    mar[1] <- mar[1] + strheight_right[1]/2*5
    mar[3] <- mar[3] + if(is.null(xlabels) & colorbar) maxstrheight_top*5 +.3 else strheight_right[1]/2*5
    mar[4] <- mar[4] + maxstrwidth_right*5
    par(mar=mar)
  }
  if(is.null(ylabels)){
    mar[2] <- mar[2] + maxstrwidth_left*5+.1
    par(mar=mar)
  }
  #par(mar = c(strheight_right[1]/2*5,maxstrwidth_left*5+.1,strheight_right[1]/2*5,maxstrwidth_right*5)+.1)

  plot(c(0,width),rep(height,2), bty="n", xlim=c(0,width), ylim=c(0,height),
       type="n", xaxs="i",
       yaxs="i", xaxt="n", yaxt="n", xlab="", ylab="", asp=asp)

  rect(sort(rep(x_breaks-diff(x_breaks)[1],K), decreasing = FALSE),
       (rep(y_breaks + diff(y_breaks)[1],M)),
       sort(rep(x_breaks,K), decreasing = FALSE),
       (rep(y_breaks ,M)),
       col = cols_to_plot,
       border = cols_to_plot)

  if(add_numbers){
    text(x = sort(rep(x_breaks-diff(x_breaks)[1]/2,K), decreasing = FALSE),
         y = rep(y_breaks + diff(y_breaks)[1]/2,M),
         labels = round(PHI,2), cex = cex.numbers)
  }

  # seperate lags
  abline_at <- NULL
  if(p==1 & K>M){
    abline_at <- M
  }else if(p>1){
    abline_at <- M*(1:(p-1))
    if(K%%M != 0){
      abline_at <- c(abline_at, M*p)
    }
  }
  if(!is.null(abline_at)){
    for(i in seq.int(length(abline_at))){
      lines(xlim_heatmap, rep(height - abline_at[i]*height/K, 2))
    }
  }

  if(!is.null(row_names)  & is.null(ylabels)){
    axis(side = 2, at = y_breaks + diff(y_breaks)[1]/2, labels = row_names, las=2, tick = FALSE,
         cex.axis = cex.axis)
  }

  if(!is.null(col_names)  & is.null(xlabels)){
    axis(side = 3, at = x_breaks - diff(x_breaks)[1]/2, labels = col_names, las=2, tick = FALSE,
         cex.axis = cex.axis)
  }

  if(colorbar){
    len <- length(colbar_range)
    xxx <- seq(0,height,length=len+1)
    rect(rep(xlim_colorbar[1], len),xxx[1:len], rep(xlim_colorbar[2], len), xxx[-1],
         col = colbar_range, border = colbar_range)

    at <- if(PHI_range[1]<0 & PHI_range[2]>0) c(0,-PHI_range[1]/(diff(PHI_range)) * height,height) else c(0,height)
    axis(side = 4, tick = FALSE, labels = colbar_labels, at, las = 2, cex = cex.colbar)
  }

  invisible(x)
}

plot_predvals <- function(preds, quantiles, observed = NULL, var_names,
                          dates, n_col, nr_insample, nr_outsample){

  quantiles <- sort(quantiles)
  even <- (length(quantiles)) %% 2 == 0

  quant_store <- apply(preds, 1:2, quantile, quantiles, na.rm=TRUE)
  if(length(quantiles)==1L){
    quant_store <- array(quant_store, c(1,nrow(quant_store), ncol(quant_store)))
  }

  M <- dim(quant_store)[3]

  lags <- if(!is.null(observed)) nrow(observed) - dim(quant_store)[2] + nr_outsample else 0

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow=c(min(5,ceiling(M/n_col)),n_col), mar = c(2,2,2,1), mgp = c(2,.5,0))
  for(i in seq.int(M)){

    data_to_plot <- matrix(as.numeric(NA), length(dates), length(quantiles)+ if(!is.null(observed)) 1 else 0)
    if(!is.null(observed)){
      data_to_plot[1:nr_insample,ncol(data_to_plot)] <- observed[,i]
    }
    data_to_plot[(lags+1):nrow(data_to_plot),1:(length(quantiles))] <- t(quant_store[,,i])
    nr_intervals <- floor(length(quantiles)/2)

    # call to ts.plot defines plot region
    ts.plot(data_to_plot, gpars = list(xlab="", xaxt="n",
                                       lty = 0
                                       #lty = c(1,rep(0, nr_intervals), if(even) rep(0, nr_intervals) else c(1,rep(0, nr_intervals))),
                                       #col=c(3,rep(1, nr_intervals), if(even) rep(1, nr_intervals) else c(2,rep(1, nr_intervals))),
                                       #lwd=c(2, rep(1, nr_intervals), if(even) rep(1, nr_intervals) else c(2,rep(1, nr_intervals))))
    )
    )
    myaxis <- axis(side = 1, labels = FALSE, tick = FALSE) + 1
    equidist <- mean(diff(myaxis))
    myaxis <- seq(equidist/2, length(dates), by = equidist)
    axis(side=1, at = myaxis, labels = dates[myaxis])
    mtext(var_names[i], side = 3)
    if(nr_insample>0){
      if(nr_intervals>0){
        for(k in seq.int(nr_intervals)){
          alpha_upper <- 0.4
          alpha_lower <- 0.2
          alphas <- seq(alpha_lower, alpha_upper, length.out = nr_intervals)
          if(nr_intervals==1){
            alphas <- alpha_upper
          }
          polygon(x = c((lags+1):nr_insample, rev((lags+1):nr_insample)),
                  y = c(quant_store[k,1:(nr_insample-lags),i],rev(quant_store[k+1,1:(nr_insample-lags),i])),
                  col = scales::alpha("red", alphas[k] ),
                  border = NA)
          if(length(quantiles)>2){
            polygon(x = c((lags+1):nr_insample, rev((lags+1):nr_insample)),
                    y = c(quant_store[length(quantiles)-k,1:(nr_insample-lags),i],rev(quant_store[length(quantiles)+1-k,1:(nr_insample-lags),i])),
                    col = scales::alpha("red", alphas[k]),
                    border = NA)
          }
        }
        if(!is.null(observed)){
          lines(1:nr_insample, observed[,i], col=1, lwd=1, lty = 4)
        }
        if(!even){
          lines((lags + 1):nr_insample, quant_store[ceiling(length(quantiles)/2),1:(nr_insample-lags),i] , col="red", lwd=2, lty = 1)
        }
      }
    }

    if(nr_outsample>0){
      for(j in seq.int(quantiles)){
        lines((nr_insample):(nr_insample + nr_outsample), quant_store[j,(nr_insample-lags):(nr_insample - lags + nr_outsample),i])
      }
    }

  }
}

#' Visualization of in-sample fit of an estimated VAR.
#'
#' @param x A `bayesianVARs_fitted` object.
#' @param dates optional vector of dates for labelling the x-axis. The default
#'   values is `NULL`; in this case, the axis will be labeled with numbers.
#' @param vars character vector containing the names of the variables to be
#'   visualized. The default is `"all"` indicating that the fit of all variables
#'   is visualized.
#' @param quantiles numeric vector indicating which quantiles to plot.
#' @param n_col integer indicating the number of columns to use for plotting.
#' @param ... Currently ignored.
#'
#' @return returns `x` invisibly
#' @seealso
#' * fitted method for class 'bayesianVARs_bvar': [fitted.bayesianVARs_bvar()].
#' * Other plotting [`plot.bayesianVARs_bvar()`],
#' [`plot.bayesianVARs_fitted()`], [`plot.bayesianVARs_predict()`],
#' [`pairs.bayesianVARs_predict()`], [`posterior_heatmap()`].
#' @export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Estimate a model
#' mod <- bvar(data, sv_keep = "all", quiet = TRUE)
#'
#' # Simulate predicted historical values including the error term.
#' pred <- fitted(mod, error_term = TRUE)
#'
#' # Visualize
#' plot(pred)
plot.bayesianVARs_fitted <- function(x,
                                     dates = NULL,
                                     vars = "all",
                                     quantiles = c(0.05,0.5,0.95),
                                     n_col = 1L, ...){

  if(length(vars)==1L & any(vars == "all")){
    vars <- 1:ncol(x$Yraw)
  }else if(any(base::isFALSE(vars %in% colnames(x$Yraw)))){
    stop("Elements of 'vars' must coincide with 'colnames(x$Yraw)'!")
  }else{
    vars <- which(colnames(x$Yraw) %in% vars)
  }

  if(is.null(dates)){
    dates <- 1:nrow(x$Yraw)
    if(!is.null(rownames(x$Yraw))){
      dates <- tryCatch(as.Date(rownames(x$Yraw)), error = function(e) dates)
    }
  }else{
    if(length(dates) != nrow(x$Yraw)){
      stop("Length of argument 'dates' differs from nrow(x$Yraw)!")
    }
  }
  dates <- as.character(dates)

  plot_predvals(x$fitted[,vars,], quantiles, observed = x$Yraw, colnames(x$Yraw[,vars]),
                dates, n_col, nrow(x$Yraw), 0)

  invisible(x)
}

#' Plot method for bayesianVARs_bvar
#'
#' Visualization of in-sample fit. Can also be used to display prediction
#' intervals of future values.
#'
#' @param x An object of class `bayesianVARs_bvar` obtained via [`bvar()`].
#' @param predictions Optional array of out of sample predictions, e.g. obtained
#'   via [`predict.bayesianVARs_bvar()`].
#' @param quantiles numeric vector indicating which quantiles to plot.
#' @param dates optional vector of dates for labelling the x-axis. The default
#'   values is `NULL`; in this case, the axis will be labeled with numbers.
#' @param n_col integer indicating the number of columns to use for plotting.
#' @param ... Currently ignored!
#' @seealso Other plotting [`plot.bayesianVARs_fitted()`],
#' [`plot.bayesianVARs_predict()`], [`pairs.bayesianVARs_predict()`],
#' [`posterior_heatmap()`].
#' @return Returns `x` invisibly.
#' @export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Estimate a model
#' mod <- bvar(data, sv_keep = "all", quiet = TRUE)
#'
#' # Simulate from posterior predictive
#' predictions <- predict(mod, ahead = 1:3)
#'
#' # Visualize
#' plot(mod, predictions = predictions)
plot.bayesianVARs_bvar <- function(x, predictions = NULL,
                                   quantiles = c(0.05,0.5,0.95),
                                   dates = NULL,
                                   n_col = 1,...){
  fit <- fitted(x, error_term = TRUE)
  nr_observed <- nrow(x$Yraw)
  if(!is.null(predictions)){
    if(!inherits(predictions, "bayesianVARs_predict")){
      stop("Argument 'predictions' must be of class 'bayesianVARs_predict' (?predict.bayesianVARs_bvar)!")
    }
    predictions <- predictions$predictions
    if(any(dim(predictions)==0L)){
      stop("Object supplied to argument 'predictions' does not contain draws of the predictive distribution. Make sure that 'simulate_predictive=TRUE' when calling predict.bayesianVARs_bvar()!")
    }
  }
  nr_outsample <- if(!is.null(predictions)) nrow(predictions) else 0L

  if(!is.null(dates)){
    if(length(dates) != nr_observed+nr_outsample){
      stop("'length(dates)' does not coincide with number of total observations (plus horizon of out-of-sample predictions)!")
    }
    dates <- as.character(dates)
  }else{
    dates <- as.character(1:(nr_observed+nr_outsample))
  }

  fitted_plus_predictions <- array(as.numeric(NA), c(nrow(fit$fitted)+nr_outsample, ncol(x$Yraw), dim(fit$fitted)[3]))
  fitted_plus_predictions[1:nrow(fit$fitted),,] <- fit$fitted
  if(!is.null(predictions)){
    fitted_plus_predictions[(nrow(fit$fitted)+1):(nrow(fit$fitted)+nr_outsample),,1:(dim(predictions)[3])] <- predictions
  }

  plot_predvals(fitted_plus_predictions, quantiles, observed = x$Yraw, colnames(x$Yraw),
                dates = dates, n_col= n_col, nr_insample = nr_observed, nr_outsample)
  #plot(fit)
  invisible(x)
}

#' @name pairs_predict
#' @title Pairwise visualization of out-of-sample posterior predictive
#'   densities.
#'
#' @param x An object of class `bayesianVARs_predict` obtained via
#'   [`predict.bayesianVARs_bvar()`].
#' @param vars Integer vector (or coercible to such) indicating which variables
#'   to plot.
#' @param ahead Integer vector (or coercible to such) indicating which step
#'   ahead to plot. `max(ahead)` must be smaller equal to
#'   `dim(x$predictions)[1]`.
#' @param ... Currently ignored!
#' @note Note that that `bayesianVARs_predict` can also be used withing [`plot.bayesianVARs_bvar()`].
#'
#' @seealso Other plotting [`plot.bayesianVARs_bvar()`],
#' [`plot.bayesianVARs_fitted()`], [`plot.bayesianVARs_predict()`]
#' [`posterior_heatmap()`].
#' @return Returns `x` invisibly.
#' @export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Estimate a model
#' mod <- bvar(data, sv_keep = "all", quiet = TRUE)
#'
#' # Simulate from posterior predictive
#' predictions <- predict(mod, ahead = 1:3)
#'
#' # Visualize
#' pairs(predictions, vars = 1:3, ahead = 1:3)
pairs.bayesianVARs_predict <- function(x, vars, ahead, ...){

  if(is.null(x$predictions)){
    stop("For plotting, set 'simulate_predictive = TRUE' when calling 'predict.bayesianVARs_bvar()'!")
  }

  vars <- as.integer(vars)
  ahead <- as.integer(ahead)
  panel.cor <- function(xx, y, digits = 2, prefix = "", cex.cor){
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(xx, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
  }

  panel.lm <- function (xx, y, col = par("col"), bg = NA, pch = par("pch"),
                        cex = 1, col.lm = 2) #, Y_obs = NA
  {
    points(xx, y, pch = pch, col = col, bg = bg, cex = cex)
    # if(!any(is.na(Y_obs))){
    #   abline(v = Y_obs[1], h = Y_obs[2])
    # }
    ok <- is.finite(xx) & is.finite(y)
    if (any(ok))
      abline(stats::lm(y ~ xx), col = col.lm)
  }

  for(i in seq_along(ahead)){
    # Y_obs <- NA
    # if(!any(is.na(x$Y_obs))){
    #   Y_obs <- x$Y_obs[i,]
    # }
    pairs(x=t(x$predictions[ahead[i],vars,]),
                            upper.panel = panel.lm,
                            lower.panel = panel.cor,
                            gap=0, row1attop=TRUE)
    mtext(paste0("t+", ahead[i]), side = 3, line = 3)
  }

  invisible(x)
}


#' Fan chart
#'
#' Visualization of (out-of-sample) predictive distribution.
#'
#' @param x An object of type `bayesianVARs_predict` obtained via
#'   [`predict.bayesianVARs_bvar()`].
#' @param dates optional vector of dates for labeling the x-axis. The default
#'   values is `NULL`; in this case, the axis will be labeled with numbers.
#' @param vars character vector containing the names of the variables to be
#'   visualized. The default is `"all"` indicating that all variables are
#'   visualized.
#' @param ahead Integer vector (or coercible to such) indicating which step
#'   ahead to plot. `max(ahead)` must be smaller equal to
#'   `dim(x$predictions)[1]`.
#' @param quantiles numeric vector indicating which quantiles to plot.
#' @param n_col integer indicating the number of columns to use for plotting.
#' @param first_obs integer indicating the first observation to be used for plotting.
#' @param ... Currently ignored!
#'
#' @return Returns `x` invisibly!
#' @seealso Other plotting [`plot.bayesianVARs_bvar()`],
#' [`plot.bayesianVARs_fitted()`], [`pairs.bayesianVARs_predict()`]
#' [`posterior_heatmap()`].
#' @export
#'
#' @examples
#' # Access a subset of the usmacro_growth dataset
#' data <- usmacro_growth[,c("GDPC1", "CPIAUCSL", "FEDFUNDS")]
#'
#' # Estimate a model
#' mod <- bvar(data, sv_keep = "all", quiet = TRUE)
#'
#' # Simulate from posterior predictive
#' predictions <- predict(mod, ahead = 1:3)
#'
#' # Visualize
#' plot(predictions, vars = 1:3, ahead = 1:3)
plot.bayesianVARs_predict <- function(x, dates = NULL, vars = "all", ahead = NULL,
                                      quantiles = c(0.05,0.25,0.5,0.75,0.95),
                                      n_col = 1L,
                                      first_obs = 1L,
                                      ...){

  Tobs <- nrow(x$Yraw)
  if(is.null(x$predictions)){
    stop("For plotting, set 'simulate_predictive = TRUE' when calling 'predict.bayesianVARs_bvar()'!")
  }

  n_ahead <- nrow(x$predictions)

  if(length(vars)==1L & any(vars == "all")){
    vars <- 1:ncol(x$Yraw)
  }else if(is.character(vars)){
    if(any(base::isFALSE(vars %in% colnames(x$Yraw)))){
      stop("Elements of 'vars' must coincide with 'colnames(x$Yraw)'!")
    }else{
      vars <- which(colnames(x$Yraw) %in% vars)
    }
  }else if(is.numeric(vars)){
    vars <- as.integer(vars)
    if(any(vars > ncol(x$Yraw))){
      stop("'max(vars)' must be at most 'ncol(x$Yraw)'!")
    }
  }
  var_names <- colnames(x$Yraw)

  if(is.null(dates)){
    dates <- first_obs:(nrow(x$Yraw) + n_ahead)
  }
  dates <- as.character(dates)

  quantiles <- sort(quantiles)
  nr_quantiles <- length(quantiles)
  nr_intervals <- floor(nr_quantiles/2)
  even <- nr_quantiles%%2 == 0

  ahead <- as.integer(ahead)
  pred_quants <- apply(x$predictions, 1:2, quantile, quantiles)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  M <- length(vars)
  par(mfrow=c(min(5,ceiling(M/n_col)),n_col), mar = c(2,2,2,1), mgp = c(2,.5,0))
  for(j in seq_along(vars)){
    plot(first_obs:(nrow(x$Yraw) + n_ahead), rep(0, length(dates)), type = "n", xlab="", ylab = "",
         xaxt="n", ylim = range(c(x$Yraw[first_obs:Tobs,vars[j]], pred_quants[,,vars[j]])))
    lines(first_obs:nrow(x$Yraw), x$Yraw[first_obs:Tobs,vars[j]])

    myaxis <- axis(side = 1, labels = FALSE, tick = FALSE) + 1
    equidist <- mean(diff(myaxis))
    myaxis <- seq(floor(equidist/2), length(dates), by = equidist)
    axis(side=1, at = myaxis + first_obs -1, labels = dates[myaxis])
    mtext(var_names[j], side = 3)

    if(nr_intervals>0){
      for(r in seq.int(nr_intervals)){
        alpha_upper <- 0.4
        alpha_lower <- 0.2
        alphas <- seq(alpha_lower, alpha_upper, length.out = nr_intervals)
        if(nr_intervals==1){
          alphas <- alpha_upper
        }
        polygon(x = c(nrow(x$Yraw):(nrow(x$Yraw)+n_ahead),
                      rev(nrow(x$Yraw):(nrow(x$Yraw)+n_ahead))),
                y = c(x$Yraw[Tobs,vars[j]], pred_quants[r,,vars[j]],
                      rev(c(x$Yraw[Tobs,vars[j]],pred_quants[r+1,,vars[j]]))),
                col = scales::alpha("red", alphas[r]),
                border = NA)
        if(length(quantiles)>2){
          polygon(x = c(nrow(x$Yraw):(nrow(x$Yraw)+n_ahead),
                        rev(nrow(x$Yraw):(nrow(x$Yraw)+n_ahead))),
                  y = c(x$Yraw[Tobs,vars[j]],
                        pred_quants[nrow(pred_quants)+1-r,,vars[j]],
                        rev(c(x$Yraw[Tobs,vars[j]],
                              pred_quants[nrow(pred_quants)-r,,vars[j]]))),
                  col = scales::alpha("red", alphas[r]),
                  border = NA)
        }
      }
      if(!even){
        lines(nrow(x$Yraw):(nrow(x$Yraw)+n_ahead),
              c(x$Yraw[Tobs,vars[j]],
                pred_quants[ceiling(length(quantiles)/2),,vars[j]]),
              col = "red", lwd = 2)
      }
    }
  }

  invisible(x)
}
