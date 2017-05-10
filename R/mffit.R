#' Fit a galaxy mass function
#'
#' This function fits an analytical model of a galaxy mass function (MF) to a discrete set of galaxies with observed masses and effective volumes.
#'
#' @param x A vector of masses or log-masses (\code{x = log10(M/Msun)}). It can also represent (log-)luminosities or absolute magnitudes, if used in conjunction with sensible initial parameters \code{p.initial}.
#' @param veff This is either a function or a vector. If veff is a function, it will be interpreted as the analytical function \code{Veff(x)} representing the effective cosmic volume, in which the log-masses \code{x} can be detected; normally this is an increasing function of \code{x}. If \code{veff} is a vector of same length as \code{x}, it will be interpreted as the values of the maximum cosmic volumes associated with the different values of \code{x}, and a function \code{Veff(x)} will be fitted automatically. In both cases the continuous function \code{Veff(x)} will be returned as \code{veff.fn} in the output list. To plot this function, use \code{\link{mfplotveff}}.
#' @param sigma Optional vector of same length as \code{x}, containing the 1-sigma uncertainties of \code{x}. If \code{x} are log-masses the \code{sigma} are symmetric Gaussian uncertanties in log-mass; otherwise \code{sigma} are symmetric Gaussian uncertainties in linear mass.
#' @param log If \code{TRUE} \code{x} and \code{sigma} refer to log-masses, otherwise to linear masses.
#' @param type Kind of MF: \code{'Schechter'} for Schechter function or \code{'PL'} for a simple power law.
#' @param niterations Number of iterations in the repeated fit-and-debias algorithm to evaluate the maximum likelihood.
#' @param resampling If \code{TRUE}, the MF is resampled \code{resampling.iterations} times to determine quantiles of the parameter distributions and fitted MF.
#' @param resampling.iterations Integer (>0) number of iterations for the resampling. Ignored if \code{resampling = FALSE}.
#' @param write.fit Logical argument specifying whether the best-fitting parameters are displayed in the console.
#' @param p.initial Initial model parameters for fitting the MF.
#' @param integration.range Range of log-mass values on which the likelihood function is evaluate and maximized. Note that this range always is in log10(mass) units, even if \code{log = FALSE}.
#'
#' @return Returns a structured list. The entry \code{input} contains all relevant input arguments. The entry \code{fit} contains all the output arguments. Most importantly, the equation of the fitted MF and the best-fitting parameters are contained in the sub-list \code{fit$parameters}.
#'
#' @keywords schechter function
#' @keywords mass function
#' @keywords fit
#'
#' @examples
#' data = mfdata()
#' mf = mffit(data$x, data$veff, data$sigma, resampling = TRUE)
#' mfplot(mf, bin.xmin = 6.5, bin.xmax = 9.5, xlim=c(2e6,5e10), ylim=c(1e-3,1))
#'
#' @author Danail Obreschkow
#'
#' @export

mffit <- function(x, # Mass or log-mass
                  veff, # either a vector of Vmax for all points x or a function Veff(x)
                  sigma = NULL, # Gaussian uncertaintiy on mass, or log-mass if log = T
                  log = TRUE, # if true, mass-values are taken as log10(mass)
                  type = 'Schechter',
                  niterations = 5,
                  resampling = FALSE,
                  resampling.iterations = 100,
                  write.fit = TRUE,
                  p.initial = NULL,
                  integration.range = seq(4,12,by=0.01)) {

  # Input handling
  if (length(x)<=2) stop('Give at least 2 data points.')
  npoints = length(x)
  if (!is.null(sigma)) {
    if (length(sigma)!=npoints) stop('Length of x.sigma must be the same as the length of x.')
  }
  if (length(integration.range)<2) stop('Invalid integration range.')
  dx = integration.range[2]-integration.range[1]
  if (any(abs(integration.range[3:length(integration.range)]-integration.range[2:(length(integration.range)-1)]-dx)>dx*1e-10)) {
    stop('integration.range must be an equally spaced vector.')
  }
  if (!log) x = log10(x)

  # Make function veff(x)
  if (is.function(veff)) {
    veff.fn = veff
  } else {
    veff.fn = .fit.veff(x,veff,sigma,log)
  }

  # Fitting
  if (is.null(p.initial)) {p.initial = mfmodel(output = 'initial', type = type)}
  best.fit = .corefit(x = x, sigma = sigma, log = log, p.initial = p.initial,
                    integration.range = integration.range, veff.fn = veff.fn,
                    type = type, niterations = niterations)

  # Initialize output parameters
  input = list(x = x, veff = veff, sigma = sigma, log = log, type = type, niterations = niterations, resampling = resampling, resampling.iterations = resampling.iterations)
  fit = list(parameters = list(equation = mfmodel(output='equation', type=type), p.optimal = best.fit$p.optimal, p.covariance = best.fit$covariance),
             fn = list(x = integration.range, y = mfmodel(x = integration.range, p = best.fit$p.optimal, type = type)),
             veff.fn = veff.fn)
  mf = list(input = input, fit = fit)

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
  nsteps = 3 # a larger number of steps leads to a more accurate sampling of the covariance ellipsoid
  y.new = array(NA,c(length(mf$fit$fn$x),nsteps^np))
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
      y.new[,i] = mfmodel(x = mf$fit$fn$x, p = p.new, type = mf$input$type)
    } else {
      y.new[,i] = NA
    }

  }

  mf$fit$fn$y.error.neg = array(NA,length(mf$fit$fn$x))
  mf$fit$fn$y.error.pos = array(NA,length(mf$fit$fn$x))
  for (i in seq(length(mf$fit$fn$x))) {
    mf$fit$fn$y.error.neg[i] = mf$fit$fn$y[i]-min(y.new[i,],na.rm=TRUE)
    mf$fit$fn$y.error.pos[i] = max(y.new[i,],na.rm=TRUE)-mf$fit$fn$y[i]
  }

  return(mf)
}

#' @export
.corefit <- function(x, sigma, log, p.initial, integration.range, veff.fn, type, niterations = 1) {

  # Input handling
  npoints = length(x)
  dx = integration.range[2]-integration.range[1]
  veff = veff.fn(integration.range)
  if (is.null(sigma)) niterations = 1
  if (!log) {unlog.range = 10^integration.range}

  # iterative algorithm
  for (k in seq(max(1,niterations))) {

    # make unbiased source density function
    rho.unbiased = array(0,length(integration.range))
    if (is.null(sigma)) {
      for (i in seq(npoints)) {
        index = which.min(abs(x[i]-integration.range))[1]
        rho.unbiased[index] = rho.unbiased[index]+1/dx
      }
    } else {
      if (niterations==0) {
        prior = array(1,length(integration.range))
      } else {
        # weighing function (= predicted source counts)
        prior = mfmodel(integration.range, p.initial, type = type)*veff*dx
      }
      for (i in seq(npoints)) {
        if (log) {
          rho.observed = exp(-(x[i]-integration.range)^2/2/sigma[i]^2)
        } else {
          # dn/dm = exp(-m^2/2/s^2)
          # m -> x = log10(m) = log(m)/log(10) => dx/dm = 1/m/log(10)
          # dn/dx = dn/dm*dm/dx = exp(-m^2/2/s^2)*m*log(10)
          rho.observed = exp(-(10^x[i]-unlog.range)^2/2/sigma[i]^2)*unlog.range
        }
        rho.corrected = rho.observed*prior
        rho.corrected = rho.corrected/(sum(rho.corrected)*dx)
        rho.unbiased = rho.unbiased+rho.corrected
      }
    }

    # make -ln(L)
    neglogL = function(p) {
      phi = mfmodel(integration.range, p, type=type)
      list = phi>0 & veff>0
      return(sum(phi[list]*veff[list]-log(phi[list])*rho.unbiased[list])*dx)
    }

    # maximize ln(L)
    opt = optim(p.initial,neglogL,hessian=TRUE,control=(parscale=abs(p.initial)))
    p.initial = opt$par
  }

  # make output
  return(list(p.optimal = opt$par, covariance = solve(opt$hessian)))

}

.fit.veff <- function(x,veff,sigma,log) {
  if (length(veff)!=length(x)) stop('Length of veff must be the same as the length of x.')
  if (min(veff)<=0) stop('All values of veff must be positive.')
  # fit ln(veff-a[3]) = a[2]*(x-a[1])
  list = veff<max(veff)
  if (is.null(sigma)) {
    weight = 1
  } else {
    if (log) {
      weight = 1/sigma[list]^2
    } else {
      delta = sigma[list]/log(10)/10^x[list]
      weight = 1/delta^2
    }
  }
  minfct <- function(a) sum((exp(a[2]*(x[list]-a[1]))+a[3]-veff[list])^2*weight)
  a = optim(c(9,10,0),minfct)$par
  minfct <- function(a) sum((1/(exp(a[2]*(x[list]-a[1]))+a[3])-1/veff[list])^2*weight)
  a = optim(a,minfct)$par
  veff.fn <- function(xval) {
    v = pmax(0,pmin(max(veff),exp(a[2]*(xval-a[1]))+a[3]))
    return(v)
  }
}

.resample <- function(mf, seed = 1) {

  # randomly resample and refit the MF
  set.seed(seed)
  np = length(mf$fit$parameters$p.optimal)
  x = mf$fit$fn$x
  density = pmax(0,mf$fit$fn$y*mf$fit$veff.fn(x))
  cum = cumsum(density/sum(density))
  npoints = length(mf$input$x)
  p.new = array(NA,c(mf$input$resampling.iterations,np))
  x.obs = array(NA,npoints)
  for (iteration in seq(mf$input$resampling.iterations)) {
    r = runif(npoints)
    for (i in seq(npoints)) {
      index = which.min(abs(cum-r[i]))
      x.obs[i] = x[index]
    }
    p.new[iteration,] = .corefit(x.obs, sigma = NULL, log = mf$input$log,
                                 p.initial = mf$fit$parameters$p.optimal,
                                 integration.range = x, veff.fn = mf$fit$veff.fn,
                                 type = mf$input$type)$p.optimal
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
    s[iteration,] = mfmodel(x,p.new[iteration,],type=mf$input$type)
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
  mf$fit$fn = append(mf$fit$fn, y.quantile)
  return(mf)

}
