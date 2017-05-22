#' Write best fitting parameters
#'
#' This function displays the best-fitting model parameters of a mass function previously fitted using \code{\link{mffit}}.
#'
#' @param mf List produced by \code{\link{mffit}}
#'
#' @examples
#' data = mfdata()
#' mf = mffit(data$mass, data$veff, write.fit = FALSE)
#' mfwrite(mf)
#'
#' @seealso See examples in \code{\link{mffit}}.
#'
#' @author Danail Obreschkow
#'
#' @export

mfwrite <- function(mf) {

  p = mf$fit$parameters$p.optimal

  if (!is.null(mf$model$mass.function.equation)) {
    cat(sprintf('%s\n',mf$model$mass.function.equation))
  }

  if (length(mf$fit$parameters$p.quantile.16)>0) {
    sigma.84 = mf$fit$parameters$p.quantile.84-mf$fit$parameters$p.optimal
    sigma.16 = mf$fit$parameters$p.optimal-mf$fit$parameters$p.quantile.16
    for (i in seq(length(p))) {
      cat(sprintf('p[%d] = %7.2f (+%3.2f -%3.2f)\n',i,p[i],sigma.84[i],sigma.16[i]))
    }
  } else {
    sigma = mf$fit$parameters$p.sigma
    for (i in seq(length(p))) {
      cat(sprintf('p[%d] = %7.2f (+-%3.2f)\n',i,p[i],sigma[i]))
    }
  }

}
