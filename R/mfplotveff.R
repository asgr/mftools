#' Plot effective volume
#'
#' This function plots the function \code{Veff(x)} used for the mass function fit, as stored in the sublist \code{fit} of the output produced by \code{\link{mffit}}.
#'
#' @param mf List produced by \code{\link{mffit}}
#'
#' @seealso See examples in \code{\link{mffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

mfplotveff <- function(mf) {
  if (mf$input$log) {
    x = mf$input$mass
  } else {
    x = log10(mf$input$mass)
  }
  plot(0,0,type='n',xlim=range(x),ylim=c(0,max(mf$model$veff.function(x))),
       pch=20,xlab='log-mass x',ylab='Veff')
  if (!is.function(mf$input$veff)) {
    points(x,mf$input$veff,pch=20)
  }
  lines(mf$fit$vectors$x,mf$model$veff.function(mf$fit$vectors$x))
}
