#' Analytical galaxy mass functions
#'
#' This function handles and evaluates different models of the galaxy mass function (MF), such as the Schechter function.
#'
#' @param x Vector of log-masses, typically \code{x = log10(M/Msun)}.
#' @param p Parameters of the analytical function. See argument \code{type}.
#' @param output Specifies what the function is doing. \code{'density'} evaluates the MF at \code{x}, \code{'npara'} returns the number of parameters of the analytical function, \code{'initial'} returns a vector of typical parameters \code{p} that are used as default initial values when fitting the MF, \code{'equation'} returns a string with the equation of the MF.
#' @param type Kind of MF: \code{'Schechter'} for Schechter function or \code{'PL'} for a simple power law.
#'
#' @examples
#' # Evaluate and plot a Schechter function
#' x = seq(7,11,length=100)
#' mass = 10^x
#' parameters = c(-2,10,-1.5)
#' phi = mfmodel(x, parameters, type = 'Schechter')
#' plot(mass, phi, type='l', log='xy')
#'
#' @seealso \code{\link{mffit}}
#'
#' @author Danail Obreschkow
#'
#' @export

mfmodel <- function(x = NULL, p = NULL, output = 'density', type = 'Schechter') {

  if (output=='npara') {
    return(length(mfmodel(output = 'initial')))
  }

  # call model function
  if (type == 'Schechter') {
    .mfmodel.Schechter(x = x, p = p, output = output)
  } else if (type == 'PL') {
    .mfmodel.PL(x = x, p = p, output = output)
  } else {
    stop('Unknown MF model.')
  }
}

.mfmodel.Schechter = function(x = NULL, p = NULL, output = 'density') {
  if (output == 'density') {
    mu = 10^(x-p[2])
    return(log(10)*10^p[1]*mu^(p[3]+1)*exp(-mu))
  } else if (output == 'initial') {
    return(c(-2,10,-1))
  } else if (output == 'equation') {
    return('dN/(dVdx) = log(10)*10^p[1]*mu^(p[3]+1)*exp(-mu), where mu=10^(x-p[2])')
  }
}

.mfmodel.PL = function(x = NULL, p = NULL, sigma = NULL, output = 'density') {
  if (output == 'density') {
    return(10^p[1]*(10^x)^p[2])
  } else if (output == 'initial') {
    return(c(2,-1))
  } else if (output == 'equation') {
    return('dN/(dVdx) = 10^p[1]*(10^x)^p[2]')
  }
}
