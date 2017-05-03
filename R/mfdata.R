#' Example data
#'
#' This function produces an example of log-masses \code{x} with uncertainties \code{xerr} and effective volumes \code{vmax}. The data is taken from a publication by Westmeier et al. (2017) and represents HI-Masses in the Sculptor structure.
#' 
#' @examples
#' data <- mfdata()
#' print(data$x)
#' print(data$vmax)
#' print(data$xerr)
#' 
#' # These data can then be used to fit a MF, e.g.
#' mf = mffit(data$x,data$vmax,data$xerr)
#' 
#' @seealso \code{\link{mffit}}
#' 
#' @author Danail Obreschkow, 2017
#' 
#' @export

mfdata <- function() {
  
  # data (Westmeier 2017, rounded)
  data = list(x = c(9.3,9.3,9.2,7.0,7.0,9.3,6.6,8.1,9.1,8.5,8.4,7.4,7.0,6.7,9.3,
                    8.0,8.1,8.9,9.4,7.4,8.4,7.5,8.5,9.1,8.1,7.4,8.2,7.5,7.7,8.0,8.0),
              xerr = c(0.1,0.1,0.2,0.2,0.1,0.1,0.5,0.1,0.2,0.1,0.1,0.1,0.1,0.4,0.1,0.2,
                       0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.1,0.2,0.2,0.3,0.2,0.1,0.1,0.1),
              vmax = c(98,98,98,36,33,98,7,98,98,98,98,70,25,11,98,98,
                       98,98,98,56,98,83,98,98,98,73,98,71,89,98,98))
}


