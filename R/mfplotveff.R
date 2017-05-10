#' Plot effective volume
#'
#' This function plots the function \code{Veff(x)} used for the mass function fit, as stored in the sublist \code{fit} of the output produced by \code{\link{mffit}}.
#'
#' @param mf List produced by \code{\link{mffit}}
#'
#' @examples
#' data = mfdata()
#' mf = mffit(data$x, data$veff, data$sigma, write.fit = FALSE)
#' mfplotveff(mf)
#'
#' @seealso \code{\link{mffit}}
#'
#' @author Danail Obreschkow
#'
#' @export

mfplotveff <- function(mf) {
  plot(0,0,type='n',xlim=range(mf$input$x),ylim=c(0,max(mf$fit$veff.fn(mf$input$x))),
       pch=20,xlab='log-mass x',ylab='Veff')
  if (!is.function(mf$input$veff)) {
    points(mf$input$x,mf$input$veff,pch=20)
    if (!is.null(mf$input$sigma)) {
      if (mf$input$log) {
        d = mf$input$sigma
      } else {
        d = mf$input$sigma/log(10)/10^mf$input$x
      }
      segments(mf$input$x-d,mf$input$veff,mf$input$x+d)
    }
  }
  lines(mf$fit$fn$x,mf$fit$veff.fn(mf$fit$fn$x))
}
