#' Fit a galaxy mass function
#'
#' This function fits an analytical model of a galaxy mass function (MF) to a discrete set of galaxies with observed masses and effective volumes.
#'
#' @param mass A vector of masses (\code{M/Msun}) or log-masses (\code{x = log10(M/Msun)}). It can also represent (log-)luminosities or absolute magnitudes, if used in conjunction with sensible initial parameters \code{p.initial}.
#' @param veff This is either a function or a vector. If veff is a function, it will be interpreted as the analytical function \code{Veff(x)} representing the effective cosmic volume, in which the log-masses \code{x = log10(M/Msun)} can be detected; normally this is an increasing function of \code{x}. If \code{veff} is a vector of same length as \code{mass}, it will be interpreted as the values of the maximum cosmic volumes associated with the different values of \code{mass}, and a function \code{Veff(x)} will be fitted automatically. In both cases the continuous function \code{Veff(x)} will be returned as \code{veff.function} in the output list \code{model}. To plot this function, use \code{\link{mfplotveff}}.
#' @param mass.sigma Optional vector of same length as \code{mass}, containing the 1-sigma uncertainties of \code{mass}. If \code{mass} are log-masses the \code{mass.sigma} are symmetric Gaussian uncertanties in log-mass; otherwise \code{mass.sigma} are symmetric Gaussian uncertainties in linear mass.
#' @param log If \code{TRUE} \code{mass} and \code{mass.sigma} refer to log-masses, otherwise to linear masses.
#' @param mass.function Eighter a string of a function. A string specifies the type of preset mass function. Currently available are \code{'Schechter'} for Schechter function or \code{'PL'} for a simple power law. Alternatively, \code{mass.function} can be a function of log-mass \code{x = log10(M/Msun)} and a list of parameters \code{p}, defined via \code{mass.function <- function(x,p) {...}}. In the latter case, the argument \code{p.initial} is mandatory.
#' @param n.iterations Number of iterations in the repeated fit-and-debias algorithm to evaluate the maximum likelihood.
#' @param resampling If \code{TRUE}, the MF is resampled \code{resampling.iterations} times to determine quantiles of the parameter distributions and fitted MF.
#' @param resampling.iterations Integer (>0) number of iterations for the resampling. Ignored if \code{resampling = FALSE}.
#' @param bias.correction If true, bias-corrected parameters are estimated using jackknifing.
#' @param write.fit Logical argument specifying whether the best-fitting parameters are displayed in the console.
#' @param p.initial Initial model parameters for fitting the MF.
#' @param integration.range Range of log-mass values on which the likelihood function is evaluate and maximized. Note that this range always is in log10(mass) units, even if \code{log = FALSE}.
#'
#' @return Returns a structured list. The sublist \code{input} contains all relevant input arguments. The sublist \code{model} specified the model function to be fitted. The sublist \code{fit} contains all the output arguments of the MLE algorithm. The output can be visualized using \code{\link{mfplot}} and \code{\link{mfwrite}}.
#'
#' @keywords schechter function
#' @keywords mass function
#' @keywords fit
#'
#' @examples
#' # basic example
#' data = mfdata()
#' mf = mffit(data$mass, data$veff)
#' mfplot(mf, xlim=c(2e6,5e10))
#'
#' # include measurement errors in fit
#' mf = mffit(data$mass, data$veff, data$mass.sigma)
#' mfplot(mf, add = TRUE, show.uncertainties = FALSE, lty = 2, col = 'purple')
#'
#' # evaluate bias correction and parameter quantiles
#' mf = mffit(data$mass, data$veff, data$mass.sigma, resampling = TRUE, resampling.iterations = 1e2, bias.correction = TRUE)
#' mfplot(mf, bin.xmin = 6.5, bin.xmax = 9.5, xlim=c(2e6,5e10), ylim=c(1e-3,1), uncertainty.type = 3, show.bias.correction = TRUE)
#'
#' # show fitted Veff function
#' mfplotveff(mf)
#'
#' @author Danail Obreschkow
#'
#' @export

mffit <- function(mass, # Mass or log-mass
                  veff, # either a vector of Vmax for all points x or a function Veff(x)
                  mass.sigma = NULL, # Gaussian uncertaintiy on mass, or log-mass if log = T
                  log = FALSE, # if true, mass-values are taken as log10(mass)
                  mass.function = 'Schechter',
                  n.iterations = 5,
                  resampling = FALSE,
                  resampling.iterations = 100,
                  bias.correction = FALSE,
                  write.fit = TRUE,
                  p.initial = NULL,
                  log.integration.range = seq(4,12,by=0.01)) {

  # Input handling
  if (length(mass)<1) stop('Give at least one data point.')
  n.data = length(mass)
  if (!is.null(mass.sigma)) {
    if (length(mass.sigma)!=n.data) stop('Length of mass.sigma must be the same as the length of mass.')
  }
  if (length(log.integration.range)<2) stop('Invalid integration range.')
  dx = log.integration.range[2]-log.integration.range[1]
  if (any(abs(log.integration.range[3:length(log.integration.range)]-log.integration.range[2:(length(log.integration.range)-1)]-dx)>dx*1e-10)) {
    stop('log.integration.range must be an equally spaced vector.')
  }
  if (!resampling) {resamling.iterations = NULL}

  # Make mass function phi(x,p)
  if (is.function(mass.function)) {
    phi = mass.function
    phi.eq = NULL
    if (is.null(p.initial)) stop('For user-defined mass function initial parameters must be given.')
  } else {
    phi <- function(x,p) {mfmodel(x,p,type=mass.function)}
    phi.eq = mfmodel(output = 'equation', type = mass.function)
    p.initial = mfmodel(output = 'initial', type = mass.function)
  }

  # Make function veff(x)
  if (is.function(veff)) {
    veff.function = veff
  } else {
    veff.function = .fit.veff(mass,veff,log)
  }

  # Fitting
  best.fit = .corefit(mass, mass.sigma, log, veff.function,
                      phi, p.initial, log.integration.range,
                      n.iterations = n.iterations)

  # Initialize output parameters
  input = list(mass = mass, veff = veff, mass.sigma = mass.sigma, log = log,
               mass.function = mass.function, n.iterations = n.iterations,
               resampling = resampling, resampling.iterations = resampling.iterations,
               p.initial = p.initial,
               log.integration.range = log.integration.range)
  model = list(mass.function = phi,
               mass.function.equation = phi.eq,
               veff.function = veff.function)
  fit = list(parameters = list(p.optimal = best.fit$p.optimal,
                               p.covariance = best.fit$covariance),
             vectors = list(x = log.integration.range,
                            y = phi(log.integration.range,best.fit$p.optimal)))
  mf = list(input = input, model = model, fit = fit)

  # Bias correction
  if (bias.correction) {mf = .bias.correction(mf)}

  # Determine uncertainties
  mf = .add.Gaussian.errors(mf)
  if (resampling) {mf = .resample(mf)}

  # Write fit on screen
  if (write.fit) mfwrite(mf)

  # Return results
  invisible(mf)

}

.add.Gaussian.errors <- function(mf) {

  # make Gaussian uncertainties of parameters
  cov = mf$fit$parameters$p.covariance
  mf$fit$parameters$p.sigma = sqrt(diag(cov))

  # make Gaussian uncertainties of MF
  eig = eigen(cov)
  np = length(mf$fit$parameters$p.optimal)
  index = 0
  nsteps = 6 # a larger number of steps leads to a more accurate sampling of the covariance ellipsoid
  y.new = array(NA,c(length(mf$fit$vectors$x),nsteps^np))
  k = seq(-1,1,length=nsteps)
  step = array(1,np)
  step[1] = 0
  for (i in seq(nsteps^np)) {

    # next step
    overflow = TRUE
    j = 1
    while (overflow) {
      step[j] = step[j]+1
      if (step[j]>nsteps) {
        step[j] = 1
        j = j+1
      } else {
        overflow = FALSE
      }
    }

    # evaluate MF
    if (any(k[step]!=0)) {
      e = k[step]/sqrt(sum(k[step]^2))
      v = array(0,np)
      for (j in seq(np)) {
        v = v+e[j]*sqrt(eig$values[j])*eig$vectors[,j]
      }
      p.new = mf$fit$parameters$p.optimal+v
      y.new[,i] = mf$model$mass.function(mf$fit$vectors$x,p.new)
      if (any(!is.finite(y.new[,i]))) {
        y.new[,i] = NA
      } else if (any(y.new[,i]<0)) {
        y.new[,i] = NA
      }

    } else {
      y.new[,i] = NA
    }

  }

  mf$fit$vectors$y.error.neg = array(NA,length(mf$fit$vectors$x))
  mf$fit$vectors$y.error.pos = array(NA,length(mf$fit$vectors$x))
  for (i in seq(length(mf$fit$vectors$x))) {
    mf$fit$vectors$y.error.neg[i] = mf$fit$vectors$y[i]-min(c(y.new[i,],Inf),na.rm=TRUE)
    mf$fit$vectors$y.error.pos[i] = max(c(y.new[i,],-Inf),na.rm=TRUE)-mf$fit$vectors$y[i]
  }

  return(mf)
}

#' @export
.corefit <- function(mass, mass.sigma, log, veff.function,
                     mass.function, p.initial, log.integration.range,
                     n.iterations = 1, debias.first.iteration = T) {

  # Input handling
  if (log) {
    x = mass
  } else {
    x = log10(mass)
    linear.integration.range = 10^log.integration.range
  }
  n.data = length(x)
  dx = log.integration.range[2]-log.integration.range[1]
  veff = pmax(0,veff.function(log.integration.range))
  if (is.null(mass.sigma)) n.iterations = 1
  if (min(x)<min(log.integration.range)) stop('Smallest mass below integration range.')
  if (max(x)>max(log.integration.range)) stop('Largest mass above integration range.')
  if (n.iterations<=0) stop('Number of iterations must be a positive integer.')

  # iterative algorithm
  for (k in seq(max(1,n.iterations))) {

    # make unbiased source density function
    rho.unbiased = array(0,length(log.integration.range))
    if (is.null(mass.sigma)) {
      for (i in seq(n.data)) {
        index = which.min(abs(x[i]-log.integration.range))[1]
        rho.unbiased[index] = rho.unbiased[index]+1/dx
      }
    } else {
      if ((n.iterations==1)&(!debias.first.iteration)) {
        prior = array(1,length(log.integration.range))
      } else {
        # weighing function (= predicted source counts)
        prior = mass.function(log.integration.range, p.initial)*veff*dx
        prior[!is.finite(prior)] = 0
      }
      for (i in seq(n.data)) {
        if (log) {
          rho.observed = exp(-(x[i]-log.integration.range)^2/2/mass.sigma[i]^2)
        } else {
          # dn/dm = exp(-m^2/2/s^2)
          # m -> x = log10(m) = log(m)/log(10) => dx/dm = 1/m/log(10)
          # dn/dx = dn/dm*dm/dx = exp(-m^2/2/s^2)*m*log(10)
          rho.observed = exp(-(mass[i]-linear.integration.range)^2/2/mass.sigma[i]^2)*linear.integration.range
        }
        rho.corrected = rho.observed*prior
        rho.corrected = rho.corrected/(sum(rho.corrected)*dx)
        rho.unbiased = rho.unbiased+rho.corrected
      }
    }

    # make -ln(L)
    neglogL = function(p) {
      phi = mass.function(log.integration.range, p)
      phi[!is.finite(phi)] = 0
      phi = pmax(1e-323,phi)
      return(sum(phi*veff-log(phi)*rho.unbiased)*dx)
    }

    # maximize ln(L)
    opt = optim(p.initial,neglogL,hessian=TRUE,control=(parscale=abs(p.initial)))
    p.initial = opt$par
  }

  # make output
  if (det(opt$hessian)<1e-12) {
    cov = array(NA,dim(opt$hessian))
  } else {
    cov = solve(opt$hessian)
  }
  return(list(p.optimal = opt$par, covariance = cov))

}

.bias.correction <- function(mf) {
  n = length(mf$input$mass)
  np = length(mf$fit$parameters$p.optimal)
  if (n>=1e3) {
    cat('Hint: bias correction normally not relevant for more than 1000 objects.\n')
  }
  veff.fn <- function(x) {mf$model$veff.function(x)*(n-1)/n}
  p.new = array(NA,c(np,n))
  for (i in seq(n)) {
    list = setdiff(seq(n),i)
    p.new[,i] = .corefit(mf$input$mass[list], mf$input$mass.sigma[list], mf$input$log,
                         veff.fn, mf$model$mass.function,
                         mf$fit$parameters$p.optimal, mf$input$log.integration.range,
                         n.iterations = 1)$p.optimal
  }
  p.reduced = apply(p.new, 1, mean, na.rm = T)
  mf$fit$parameters$p.optimal.bias.corrected = n*mf$fit$parameters$p.optimal-(n-1)*p.reduced
  return(mf)
}

.fit.veff <- function(mass,veff,log) {
  if (length(veff)!=length(mass)) stop('Length of veff must be the same as the length of mass.')
  if (min(veff)<=0) stop('All values of veff must be positive.')
  if (log) {
    x = mass
  } else {
    x = log10(mass)
  }

  # fit ln(veff-a[3]) = a[2]*(x-a[1])
  list = veff<max(veff)
  minfct <- function(a) sum((exp(a[2]*(x[list]-a[1]))+a[3]-veff[list])^2)
  a = optim(c(9,10,0),minfct)$par
  minfct <- function(a) sum((1/(exp(a[2]*(x[list]-a[1]))+a[3])-1/veff[list])^2)
  a = optim(a,minfct)$par
  veff.function <- function(xval) {
    v = pmax(0,pmin(max(veff),exp(a[2]*(xval-a[1]))+a[3]))
    return(v)
  }
  return(veff.function)
}

.resample <- function(mf, seed = 1) {

  # randomly resample and refit the MF
  set.seed(seed)
  np = length(mf$fit$parameters$p.optimal)
  x = mf$input$log.integration.range
  density = pmax(0,mf$fit$vectors$y*mf$model$veff.function(x))
  cum = cumsum(density/sum(density))
  n.data = length(mf$input$mass)
  p.new = array(NA,c(mf$input$resampling.iterations,np))
  for (iteration in seq(mf$input$resampling.iterations)) {
    n.new = max(2,rpois(1,n.data))
    x.obs = array(NA,n.new)
    r = runif(n.new)
    for (i in seq(n.new)) {
      index = which.min(abs(cum-r[i]))
      x.obs[i] = x[index]
    }
    p.new[iteration,] = .corefit(mass = x.obs, mass.sigma = NULL, log = T,
                                 veff.function = mf$model$veff.function,
                                 mass.function = mf$model$mass.function,
                                 p.initial = mf$fit$parameters$p.optimal,
                                 log.integration.range = mf$input$log.integration.range,
                                 n.iterations = 1)$p.optimal
  }

  # make parameter quantiles
  q = c(0.05,0.16,0.25,0.75,0.84,0.95)
  p.quant = array(NA,c(length(q),np))
  for (i in seq(np)) {
    p.quant[,i] = quantile(p.new[,i],q,names=FALSE)
  }
  p.quantile = list(p.quantile.05 = p.quant[1,], p.quantile.16 = p.quant[2,],
                    p.quantile.25 = p.quant[3,], p.quantile.75 = p.quant[4,],
                    p.quantile.84 = p.quant[5,], p.quantile.95 = p.quant[6,])

  # make MF quantiles
  s = array(NA,c(mf$input$resampling.iterations,length(x)))
  for (iteration in seq(mf$input$resampling.iterations)) {
    s[iteration,] = mf$model$mass.function(x,p.new[iteration,])
  }
  y.quant = array(NA,c(6,length(x)))
  for (i in seq(length(x))) {
    list = !is.na(s[,i]) & is.finite(s[,i]) & (s[,i]>0)
    y.quant[,i] = quantile(s[list,i],q,names=FALSE)
  }
  y.quantile = list(y.quantile.05 = y.quant[1,], y.quantile.16 = y.quant[2,],
                    y.quantile.25 = y.quant[3,], y.quantile.75 = y.quant[4,],
                    y.quantile.84 = y.quant[5,], y.quantile.95 = y.quant[6,])

  # output parameters
  mf$fit$parameters = append(mf$fit$parameters, p.quantile)
  mf$fit$vectors = append(mf$fit$vectors, y.quantile)
  return(mf)

}

.derivative <- function(f,p,d,delta=array(0.1,length(p))) {
  dp = array(0,length(p))
  dp[d[1]] = delta[d[1]]/2
  if (length(d)>1) {
    df = derivative(f,p+dp,d[2:length(d)],delta)-derivative(f,p-dp,d[2:length(d)],delta)
  } else {
    df = f(p+dp)-f(p-dp)
  }
  return(df/delta[d[1]])
}

.third_derivative <- function(f,p) {
  library(pracma)
  n = length(p)
  T = array(NA,c(n,n,n))
  for (i in seq(n)) {
    for (j in seq(n)) {
      for (k in seq(n)) {
        H = function(x) {
          s = array(0,n)
          s[i] = x
          hessian(f,p+s)[j,k]
        }
        T[i,j,k] = numderiv(H,0)$df
      }
    }
  }
  return(T)
}
