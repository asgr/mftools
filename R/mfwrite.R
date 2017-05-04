#' Write best fitting parameters
#'
#' This function displays the best-fitting model parameters of a mass function previously fitted using \code{\link{mffit}}.
#'
#' @param mf List produced by \code{\link{mffit}}
#'
#' @examples
#' data = mfdata()
#' mf = mffit(data$x, data$vmax, data$xerr, write.fit = FALSE)
#' mfwrite(mf)
#'
#' @seealso \code{\link{mffit}}
#'
#' @author Danail Obreschkow, 2017
#'
#' @export

mfwrite <- function(mf) {
  if (length(mf$fit$parameters$p.quantile.16)>0) {
    mfmodel(p = mf$fit$parameters$p.optimal,
            sigma = rbind(mf$fit$parameters$p.quantile.84-mf$fit$parameters$p.optimal,
                          mf$fit$parameters$p.optimal-mf$fit$parameters$p.quantile.16),
            type = mf$input$type, output='parameters')
  } else {
    mfmodel(p = mf$fit$parameters$p.optimal,
            sigma = mf$fit$parameters$p.sigma,
            type = mf$input$type, output='parameters')
  }
}
