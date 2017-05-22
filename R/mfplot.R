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
#' @param show.uncertainties If \code{TRUE}, uncertainties are displayed around the best fit model.
#' @param uncertainty.type \code{1}: plot Gaussian 1-sigma uncertanties propagated from the Hessian matrix of the likelihood. \code{2}: plot 68 percentile region (from 16 to 84 percent). \code{3} plot 50 (25 to 70) and 90 (5 to 95) percentile regions.
#' @param add If \code{TRUE}, the lines are overplotted on the currently open plot.
#' @param col Color of ML fit and uncertainty regions.
#' @param lwd Line width of ML fit.
#' @param lty Line type of ML fit.
#' @param col.bias.correction Line properties of bias corrected MF fit.
#' @param lwd.bias.correction Line properties of bias corrected MF fit.
#' @param lty.bias.correction Line properties of bias corrected MF fit.
#'
#' @return If used as \code{mf = mfplot(mf)} (with \code{nbins=NULL} or \code{nbins>0}), the list \code{mf} is appended an additional entry \code{bin}, containing the binned galaxy data points.
#'
#' @seealso See examples in \code{\link{mffit}}.
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
                   show.uncertainties = TRUE,
                   uncertainty.type = NULL,
                   show.bias.correction = FALSE,
                   add = FALSE,
                   col = 'red',
                   lwd = 1.5,
                   lty = 1,
                   col.bias.correction = col,
                   lwd.bias.correction = lwd,
                   lty.bias.correction = 2) {

  r = col2rgb(col)[1]/255
  g = col2rgb(col)[2]/255
  b = col2rgb(col)[3]/255

  # mass handling
  if (mf$input$log) {
    x = mf$input$mass
  } else {
    x = log10(mf$input$mass)
  }

  # make binned data
  make.bins = is.null(nbins) || (!is.na(nbins) && (nbins>0))
  if (make.bins) {
    bin = list()
    if (is.null(nbins)) {
      bin$n = min(100,round(sqrt(length(x))))
    } else {
      bin$n = nbins
    }
    if (is.null(bin.xmin)) {
      bin$xmin = min(x)-(max(x)-min(x))/bin$n*0.25
    } else {
      bin$xmin = bin.xmin
    }
    if (is.null(bin.xmax)) {
      bin$xmax = max(x)+(max(x)-min(x))/bin$n*0.25
    } else {
      bin$xmax = bin.xmax
    }
    bin$wx = bin$xmax-bin$xmin
    bin$dx = bin$wx/bin$n
    bin$xcenter = bin$xmin+(seq(bin$n)-0.5)*bin$dx

    # fill data into bins
    if (is.function(mf$input$veff)) {
      v = mf$input$veff(x)
    } else {
      v = mf$input$veff
    }
    bin$phi = bin$count = bin$xmean = array(0,bin$n)
    for (i in seq(length(x))) {
      k = floor((x[i]-bin$xmin)/bin$wx*0.99999999*bin$n)+1
      bin$phi[k] = bin$phi[k]+1/bin$dx/v[i]
      bin$count[k] = bin$count[k]+1
      bin$xmean[k] = bin$xmean[k]+x[i]
    }
    bin$xmean = bin$xmean/pmax(1,bin$count)
    mf = append(mf,list(bin=bin))
  }

  # define plot limits
  if (is.null(xlim)) {
    xlim = 10^range(mf$fit$vectors$x)
  }
  if (is.null(ylim)) {
    if (make.bins) {
      ylim = c(1e-3*max(mf$fit$vectors$y),2*max(max(mf$fit$vectors$y),bin$phi))
    } else {
      ylim = c(1e-3*max(mf$fit$vectors$y),2*max(mf$fit$vectors$y))
    }
  }

  # open plot
  if (add==FALSE) {
    plot(1,1,type='n',log='xy',xaxs='i',yaxs='i',xaxt='n',yaxt='n',
         xlim = xlim, ylim = ylim, xlab = '', ylab = '',bty='n')
  }

  # plot polygons
  if (show.uncertainties) {
    poly.x = 10^c(mf$fit$vectors$x,rev(mf$fit$vectors$x))
    if (is.null(uncertainty.type)) {
      if (length(mf$fit$vectors$y.quantile.05)>0) {
        uncertainty.type = 2
      } else {
        uncertainty.type = 1
      }
    }
    if ((uncertainty.type>1) & (!length(mf$fit$vectors$y.quantile.05)>0)) stop('Quantiles not available. Use resampling in mffit.')
    if (uncertainty.type == 3) {
      poly.y.50 = pmax(ylim[1],c(mf$fit$vectors$y.quantile.25,rev(mf$fit$vectors$y.quantile.75)))
      poly.y.90 = pmax(ylim[1],c(mf$fit$vectors$y.quantile.05,rev(mf$fit$vectors$y.quantile.95)))
      polygon(poly.x,poly.y.90,col=rgb(r,g,b,0.2),border=NA)
      polygon(poly.x,poly.y.50,col=rgb(r,g,b,0.3),border=NA)
    } else {
      if (uncertainty.type == 2) {
        poly.y.68 = pmax(ylim[1],c(mf$fit$vectors$y.quantile.16,rev(mf$fit$vectors$y.quantile.84)))
      } else if (uncertainty.type == 1) {
        poly.y.68 = pmax(ylim[1],c(mf$fit$vectors$y-mf$fit$vectors$y.error.neg,
                                   rev(mf$fit$vectors$y+mf$fit$vectors$y.error.pos)))
      }
      polygon(poly.x,poly.y.68,col=rgb(r,g,b,0.3),border=NA)
    }
  }

  # plot central fit
  x = mf$fit$vectors$x[mf$fit$vectors$y>0]
  lines(10^x,mf$fit$vectors$y[mf$fit$vectors$y>0],col=col,lwd=lwd,lty=lty)
  if (show.bias.correction) {
    if (length(mf$fit$parameters$p.optimal.bias.corrected)==0) stop('Bias corrected MLE parameters not available. Use bias.correction in mffit.')
    lines(10^x,mf$model$mass.function(x,mf$fit$parameters$p.optimal.bias.corrected),col=col.bias.correction,lwd=lwd.bias.correction,lty=lty.bias.correction)
  }

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
  if (!add) {
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
  }

  box(which = "plot", lty = "solid", lwd = 1)

  invisible(mf)
}

if (requireNamespace("magicaxis", quietly = TRUE)) {
  # to avoid conflicts if that some packages take a while to load
  suppressWarnings(suppressMessages(library(magicaxis, quietly = TRUE, verbose = FALSE)))
}
