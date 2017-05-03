#' Plot effective volume
#'
#' This function plots the function \code{Vmax(x)} used for the mass function fit, as stored in the sublist \code{fit} of the list produced by \code{\link{mffit}}.
#' 
#' @param mf List produced by \code{\link{mffit}}
#' 
#' @examples
#' data = mfdata()
#' mf = mffit(data$x, data$vmax, data$xerr, write.fit = FALSE)
#' mfplotvmax(mf)
#' 
#' @seealso \code{\link{mffit}}
#' 
#' @author Danail Obreschkow, 2017
#' 
#' @export

mfplotvmax <- function(mf) {
  plot(0,0,type='n',xlim=range(mf$input$x),ylim=range(mf$input$vmax),
       pch=20,xlab='Mass scale x',ylab='Vmax')
  if (!is.function(mf$input$vmax)) {
    points(mf$input$x,mf$input$vmax,pch=20)
    if (!is.null(mf$input$xerr)) {
      segments(mf$input$x-mf$input$xerr,mf$input$vmax,mf$input$x+mf$input$xerr)
    }
  }
  lines(mf$fit$plot$x,mf$fit$vmax.fn(mf$fit$plot$x))
}