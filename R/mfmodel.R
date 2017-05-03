#' Analytical galaxy mass functions
#'
#' This function handles and evaluates different models of the galaxy mass function (MF), such as the Schechter function.
#' 
#' @param x Vector of log-masses, typically \code{x = log10(M/Msun)}.
#' @param p Parameters of the analytical function. See argument \code{type}.
#' @param sigma Uncertainties of the parameters \code{p}. If sigma is a vector of same length as \code{p}, these uncertainties are considered to be 1-sigma Gaussian uncertainties. If sigma is a 2-row array, the first row contains the upper uncertanties to the 84-percentile value, whereas the second row contains the lower uncertainties down to the 16-percentile value.  
#' @param output Specifies what the function is doing. \code{'density'} evaluate the MF at \code{x}, \code{'npara'} returns the number of parameters of the analytical function, \code{'default'} returns a vector of typical parameters \code{p} that can be used as initial values for fitting the MF, \code{'equation'} returns a string with the equation of the MF, \code{'parameters'} displays the best fitting parameters with uncertainties.
#' @param type Kind of MF: \code{'Schechter'} for Schechter function or \code{'PL'} for a simple power law.
#' 
#' @examples
#' # Evaluate and plot a Schechter function
#' x = seq(6,10,length=100)
#' mass = 10^x
#' parameters = c(-2,9.6,-1.5)
#' phi = mfmodel(x, parameters, type = 'Schechter')
#' plot(mass, phi, type='l', log='xy')
#' 
#' # Display the equation and best-fitting parameters wit uncertainties
#' uncertainties = c(0.3,0.23,0.13)
#' mfmodel(p = parameters, sigma = uncertainties, output = 'parameters')
#' 
#' @seealso \code{\link{mffit}}
#' 
#' @author Danail Obreschkow, 2017
#' 
#' @export

mfmodel <- function(x = NULL, p = NULL, sigma = NULL, output = 'density', type = 'Schechter') {
  if (output!='npara' & output!='default' & output!='equation') {
    np = mfmodel(output='npara', type = type)
    if (length(p)!=np) stop('Wrong number of model parameters.')
    if (!is.null(sigma)) {
      if (!is.array(sigma) & length(sigma)!=np) stop('Wrong number of model parameter uncertainties sigma.')
      if (is.array(sigma) && dim(sigma)[2]!=np) stop('Wrong number of model parameter uncertainties sigma.')
      if (min(sigma)<=0) stop('Model parameter uncertainties sigma must be non-negative.')
    }
  }
  if (output=='parameters') {
    cat(sprintf('MODEL: %s\n',mfmodel(type = type, output = 'equation')))
  }
  if (type == 'Schechter') {
    .mfmodel.Schechter(x = x, p = p, sigma = sigma, output = output)
  } else if (type == 'PL') {
    .mfmodel.PL(x = x, p = p, sigma = sigma, output = output)
  } else {
    stop('Unknown mf model.')
  }
}

.mfmodel.Schechter = function(x = NULL, p = NULL, sigma = NULL, output = 'density') {
  log10phi_star = p[1]
  log10m_star   = p[2]
  alpha         = p[3]
  if (output == 'density') {
    mu = 10^(x-log10m_star)
    density = log(10)*10^log10phi_star*mu^(alpha+1)*exp(-mu)
    return(density)
  } else if (output == 'npara') {
    return(3)
  } else if (output == 'default') {      
    return(c(-2,10,-1))
  } else if (output == 'equation') {
    return('phi(M) [Mpc^-3 dex^-1] = log(10)*phi_star*(M/Mstar)^(alpha+1)*exp(-M/M_star)')
  } else if (output == 'parameters') {
    if (is.null(sigma)) {
      cat(sprintf('log10(phi_star) = %9.2f\n',p[1]))
      cat(sprintf('log10(M_star)   = %9.2f\n',p[2]))
      cat(sprintf('alpha           = %9.2f\n',p[3]))
    } else if (is.array(sigma)) {
      cat(sprintf('log10(phi_star) = %9.2f (+%3.2f -%3.2f)\n',p[1],sigma[1,1],sigma[2,1]))
      cat(sprintf('log10(M_star)   = %9.2f (+%3.2f -%3.2f)\n',p[2],sigma[1,2],sigma[2,2]))
      cat(sprintf('alpha           = %9.2f (+%3.2f -%3.2f)\n',p[3],sigma[1,3],sigma[2,3]))
    } else {
      cat(sprintf('log10(phi_star) = %9.2f (+-%3.2f)\n',p[1],sigma[1]))
      cat(sprintf('log10(M_star)   = %9.2f (+-%3.2f)\n',p[2],sigma[2]))
      cat(sprintf('alpha           = %9.2f (+-%3.2f)\n',p[3],sigma[3]))
    }
  } else {
    stop('Unknown output.')
  }
}

.mfmodel.PL = function(x = NULL, p = NULL, sigma = NULL, output = 'density') {
  log10scale = p[1]
  alpha      = p[2]
  if (output == 'density') {
    density = 10^log10scale*(10^x)^alpha
    return(density)
  } else if (output == 'npara') {
    return(2)
  } else if (output == 'default') {      
    return(c(2,-1))
  } else if (output == 'equation') {
    return('phi(M) [Mpc^-3 dex^-1] = scale*(M/Msun)^alpha')
  } else if (output == 'parameters') {
    if (is.null(sigma)) {
      cat(sprintf('phi_star  = %9.2e\n',p[1]))
      cat(sprintf('alpha     = %9.2f\n',p[2]))
    } else if (is.array(sigma)) {
      cat(sprintf('phi_star  = %9.2f (+%3.2f -%3.2f)\n',p[1],sigma[1,1],sigma[2,1]))
      cat(sprintf('alpha     = %9.2f (+%3.2f -%3.2f)\n',p[2],sigma[1,2],sigma[2,2]))
    } else {
      cat(sprintf('phi_star  = %9.2f (+-%3.2f)\n',p[1],sigma[1]))
      cat(sprintf('alpha     = %9.2f (+-%3.2f)\n',p[2],sigma[2]))
    }
  } else {
    stop('Unknown output.')
  }
}