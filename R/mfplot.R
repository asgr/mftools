#' Display fitted galaxy mass function
#'
#' This function displays the galaxy mass function (MF) fitted using \code{\link{mffit}}.
#'
#' @param mf List produced by \code{\link{mffit}}
#' @param nbins Number of bins to be plotted. This is purely for illustrative purposes. The fitting does not use bins. Choose \code{nbins=NULL} (default) to determine the number of bins automatically or \code{nbins=0} to suppress the bins.
#' @param bin.xmin Left edge of first bin
#' @param bin.xmax Right edge of last bin
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param xlim 2-element vector with x-axis plotting limits
#' @param ylim 2-element vector with y-axis plotting limits
#' @param show.quantile.50.90 If TRUE, 50 (25 to 70) and 90 (5 to 95) percentile regions are plotted, otherwise only the 68 (16 to 84) percentile region is shown. Percentiles are only available if the MF has been fittet using \code{resampling = TRUE} in \code{\link{mffit}}. Otherwise the MF will be displayed with the 68 percent confidence region determined from direct forward-propagation of the parameter covariance matrix.
#'
#' @return If used as \code{mf = mfplot(mf)} (with \code{nbins=NULL} or \code{nbins>0}), the list \code{mf} is appended an additional entry \code{bin}, containing the binned galaxy data points.
#'
#' @examples
#' data = mfdata()
#' mf = mffit(data$x, data$veff, data$sigma, write.fit = FALSE)
#' mfplot(mf, nbins = 12, bin.xmin = 6.5, bin.xmax = 9.5, xlim=c(2e6,5e10), ylim=c(1e-3,2))
#'
#' @seealso \code{\link{mffit}}
#'
#' @author Danail Obreschkow
#'
#' @export

mfplot <- function(mf,
                   nbins = NULL,
                   bin.xmin = NULL,
                   bin.xmax = NULL,
                   xlab = expression('M [M'['sun']*']'),
                   ylab = expression(phi~'[Mpc'^-3~'dex'^-1~']'),
                   xlim = NULL,
                   ylim = NULL,
                   show.quantile.50.90 = FALSE) {

  # make binned data
  make.bins = is.null(nbins) || (!is.na(nbins) && (nbins>0))
  if (make.bins) {
    bin = list()
    if (is.null(nbins)) {
      bin$n = min(100,round(sqrt(length(mf$input$x))))
    } else {
      bin$n = nbins
    }
    if (is.null(bin.xmin)) {
      bin$xmin = min(mf$input$x)-(max(mf$input$x)-min(mf$input$x))/bin$n*0.25
    } else {
      bin$xmin = bin.xmin
    }
    if (is.null(bin.xmax)) {
      bin$xmax = max(mf$input$x)+(max(mf$input$x)-min(mf$input$x))/bin$n*0.25
    } else {
      bin$xmax = bin.xmax
    }
    bin$wx = bin$xmax-bin$xmin
    bin$dx = bin$wx/bin$n
    bin$xcenter = bin$xmin+(seq(bin$n)-0.5)*bin$dx

    # fill data into bins
    if (is.function(mf$input$veff)) {
      v = mf$input$veff(mf$input$x)
    } else {
      v = mf$input$veff
    }
    bin$phi = bin$count = bin$xmean = array(0,bin$n)
    for (i in seq(length(mf$input$x))) {
      k = floor((mf$input$x[i]-bin$xmin)/bin$wx*0.99999999*bin$n)+1
      bin$phi[k] = bin$phi[k]+1/bin$dx/v[i]
      bin$count[k] = bin$count[k]+1
      bin$xmean[k] = bin$xmean[k]+mf$input$x[i]
    }
    bin$xmean = bin$xmean/pmax(1,bin$count)
    mf = append(mf,list(bin=bin))
  }

  # define plot limits
  if (is.null(xlim)) {
    xlim = 10^range(mf$fit$fn$x)
  }
  if (is.null(ylim)) {
    ylim = c(1e-3*max(mf$fit$fn$y),2*max(max(mf$fit$fn$y),bin$phi))
  }

  # open plot
  plot(1,1,type='n',log='xy',xaxs='i',yaxs='i',xaxt='n',yaxt='n',
       xlim = xlim, ylim = ylim, xlab = '', ylab = '',bty='n')

  # plot polygons
  poly.x = 10^c(mf$fit$fn$x,rev(mf$fit$fn$x))
  show.quantile.50.90 = show.quantile.50.90 & length(mf$fit$fn$y.quantile.05)>0
  if (show.quantile.50.90) {
    poly.y.50 = pmax(ylim[1],c(mf$fit$fn$y.quantile.25,rev(mf$fit$fn$y.quantile.75)))
    poly.y.90 = pmax(ylim[1],c(mf$fit$fn$y.quantile.05,rev(mf$fit$fn$y.quantile.95)))
    polygon(poly.x,poly.y.90,col='#ffdddd',border=NA)
    polygon(poly.x,poly.y.50,col='#ffaaaa',border=NA)
  } else {
    if (length(mf$fit$fn$y.quantile.16)>0) {
      poly.y.68 = pmax(ylim[1],c(mf$fit$fn$y.quantile.16,rev(mf$fit$fn$y.quantile.84)))
    } else {
      poly.y.68 = pmax(ylim[1],c(mf$fit$fn$y-mf$fit$fn$y.error.neg,
                                 rev(mf$fit$fn$y+mf$fit$fn$y.error.pos)))
    }
    polygon(poly.x,poly.y.68,col='#ffbbbb',border=NA)
  }

  # plot central fit
  lines(10^mf$fit$fn$x[mf$fit$fn$y>0],
        mf$fit$fn$y[mf$fit$fn$y>0],col='red',lwd=1.5)

  # plot binned data
  if (make.bins) {
    list = bin$phi>0
    points(10^bin$xmean[list],bin$phi[list],pch=20)
    f.16 = pmax(1e-3,qpois(0.16,bin$count[list])/bin$count[list])
    f.84 = qpois(0.84,bin$count[list])/bin$count[list]
    segments(10^bin$xmean[list],bin$phi[list]*f.16,10^bin$xmean[list],bin$phi[list]*f.84)
    segments(10^(bin$xmin+seq(0,bin$n-1)[list]*bin$dx),bin$phi[list],
             10^(bin$xmin+seq(1,bin$n)[list]*bin$dx),bin$phi[list])
  }

  # axes
  if (requireNamespace("magicaxis", quietly = TRUE)) {
    suppressWarnings(suppressMessages(library(magicaxis, quietly = TRUE, verbose = FALSE)))
    suppressWarnings(suppressMessages(library(magicaxis, quietly = TRUE, verbose = FALSE)))
    magaxis(side=1,xlab=xlab,lwd=NA,lwd.ticks=1)
    magaxis(side=2,ylab=ylab,lwd=NA,lwd.ticks=1)
    magaxis(side=3,labels=FALSE,lwd=NA,lwd.ticks=1)
    magaxis(side=4,labels=FALSE,lwd=NA,lwd.ticks=1)
  } else {
    axis(side=1,xlab=xlab,lwd=NA,lwd.ticks=1)
    axis(side=2,ylab=ylab,lwd=NA,lwd.ticks=1)
  }

  box(which = "plot", lty = "solid", lwd = 1)

  invisible(mf)
}

if (requireNamespace("magicaxis", quietly = TRUE)) {
  # to avoid conflicts if that some packages take a while to load
  suppressWarnings(suppressMessages(library(magicaxis, quietly = TRUE, verbose = FALSE)))
}
