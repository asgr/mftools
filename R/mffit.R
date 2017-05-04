#' Fit a galaxy mass function
#'
#' This function fits an analytical model of a galaxy mass function (MF) to a discrete set of galaxies with observed masses and effective volumes.
#'
#' @param x Normally a vector of log-masses \code{x = log10(M/Msun)}, but also be used with log-luminosities or absolute magnitudes.
#' @param vmax This is either a function or a vector. If vmax is a function, it will be interpreted as the analytical function \code{Veff(x)} representing the effective cosmic volume, in which the log-masses \code{x} can be detected; normally this is an increasing function of \code{x}. If \code{vmax} is a vector of same length as \code{x}, it will be interpreted as the values of the maximum cosmic volumes associated with the different values of \code{x}, and a function \code{Veff(x)} will be fitted automatically. In both cases the continuous function \code{Vmax(x)} will be returned as \code{veff.fn} in the output list.
#' @param x.sigma Optional vector of same length as \code{x}, containing the 1-sigma uncertainties of \code{x}.
#' @param type Kind of MF: \code{'Schechter'} for Schechter function or \code{'PL'} for a simple power law.
#' @para niterations Number of iterations in the repeated fit-and-debias algorithm to evaluate the maximum likelihood.
#' @param resampling If true, the MF is resampled resampling.iterations times to determine quantiles of the parameter distributions and fitted functions.
#' @param resampling.iterations Integer (>0) number of iterations for the resampling. Ignored if \code{resampling = FALSE}.
#' @param write.fit Logical argument specifying whether the best-fitting parameters are displayed in the console.
#'
#' @return Returns a structured list. The entry \code{input} contains all the input arguments. The entry \code{fit} contains all the output arguments. Most importantly, the equation of the fitted MF and the best-fitting parameters are contained in the sub-list \code{fit$parameters}.
#'
#' @keywords schechter function
#' @keywords mass function
#' @keywords fit
#'
#' @examples
#' data = mfdata()
#' mf = mffit(data$x, data$vmax, data$x.sigma, resampling = TRUE)
#' mfplot(mf, bin.xmin = 6.5, bin.xmax = 9.5, xlim=c(2e6,5e10), ylim=c(1e-3,1))
#'
#' @author Danail Obreschkow, 2017
#'
#' @export

mffit <- function(x,
                  vmax, # either a vector of Vmax for all points x or a function Vmax(x)
                  x.sigma = NULL,
                  type = 'Schechter',
                  niterations = 5,
                  resampling = FALSE,
                  resampling.iterations = 100,
                  write.fit = TRUE) {

  # Input handling
  if (length(x)<=2) stop('give at least 2 data points')
  npoints = length(x)
  wx = max(x)-min(x)

  # make x-range
  if (!is.null(x.sigma)) {
    if (length(x.sigma)!=npoints) stop('Length of x.sigma must be the same as the length of x.')
  }
  if (is.null(x.sigma)) {
    xmin = min(x)-wx
    xmax = max(x)+wx
  } else {
    xmin = min(x-x.sigma*2)-wx
    xmax = max(x+x.sigma*2)+wx
  }
  dx = min(0.01,(xmax-xmin)/500)
  xrange = seq(xmin+dx/2,xmax,by=dx)

  # Make Vmax function
  if (is.function(vmax)) {
    veff.fn = vmax
  } else {
    if (length(vmax)!=npoints) stop('Length of vmax must be the same as the length of x.')
    if (min(vmax)<=0) stop('All values of vmax must be positive.')
    # fit ln(vmax-a[3]) = a[2]*(x-a[1])
    list = vmax<max(vmax)
    if (is.null(x.sigma)) {
      minfct <- function(a) {
        return(sum((exp(a[2]*(x[list]-a[1]))+a[3]-vmax[list])^2))
      }
    } else {
      minfct <- function(a) {
        return(sum((exp(a[2]*(x[list]-a[1]))+a[3]-vmax[list])^2/x.sigma[list]))
      }
    }
    a = optim(c(9,10,0),minfct)$par
    veff.fn <- function(xval) {
      v = pmax(0,pmin(max(vmax),exp(a[2]*(xval-a[1]))+a[3]))
      v[xval<xmin] = 0
      return(v)
    }
  }

  # Make source density function
  best.fit = .corefit(x.obs = x, x.sigma = x.sigma,
                      p.initial = mfmodel(output = 'default', type = type),
                      x = xrange, v.eff = veff.fn(xrange), type = type, nit = niterations)

  # Determine 1-sigma range of model function
  eig = eigen(best.fit$covariance)
  np = length(best.fit$p.optimal)
  index = 0
  phi.new = array(NA,c(length(xrange),np*2))
  for (i in seq(np)) {
    for (k in c(-1,1)) {
      index = index+1
      v = k*sqrt(eig$values[i])*eig$vectors[,i]
      p.new = best.fit$p.optimal+v
      phi.new[,index] = mfmodel(x = xrange, p = p.new, type = type)
    }
  }
  phi.sigma = array(NA,length(xrange))
  for (i in seq(length(xrange))) {
    phi.sigma[i] = (max(phi.new[i,])-min(phi.new[i,]))/2
  }

  # Initialize output parameters
  input = list(x = x, x.sigma = x.sigma, vmax = vmax, type = type, niterations = niterations,
               resampling = resampling, resampling.iterations = resampling.iterations)
  fit = list(parameters = list(equation = mfmodel(output='equation', type=type), p.optimal = best.fit$p.optimal, p.sigma = best.fit$sigma, p.cov = best.fit$covariance),
             plot = list(x = xrange, phi.optimal = mfmodel(x = xrange, p = best.fit$p.optimal, type = type),
                         phi.sigma = phi.sigma),
             veff.fn = veff.fn)
  mf = list(input = input, fit = fit)

  # Resample model
  if (resampling) {
    output.resample = .mfresample(mf)
    mf$fit$parameters = append(mf$fit$parameters,output.resample$p.quantile)
    mf$fit$plot = append(mf$fit$plot,output.resample$phi.quantile)
  }

  # Write fit on screen
  if (write.fit) mfwrite(mf)

  # Return results
  invisible(mf)

}

#' @export
.corefit <- function(x.obs,x.sigma,p.initial,x,v.eff,type,nit) {

  # checks
  if (!is.null(x.sigma)) {
    if (length(x.sigma)!=length(x.obs)) {stop('x.obs and x.sigma must have the same length.')}
    if (min(x.sigma)<=0) {stop('All values of x.sigma must be larger than 0.')}
  }

  # some intro stuff
  dx = x[2]-x[1]
  if (is.null(x.sigma)) {nit = 1}

  # iterative algorithm
  for (k in seq(max(1,nit))) {

    # make unbiased source density function
    rho.unbiased = x*0
    n.sources = length(x.obs)
    if (is.null(x.sigma)) {
      for (i in seq(n.sources)) {
        index = which.min(abs(x.obs[i]-x))[1]
        rho.unbiased[index] = rho.unbiased[index]+1/dx
      }
    } else {

      # weighing function
      if (nit==0) {
        rho.initial = x*0+1
      } else {
        rho.initial = mfmodel(x, p.initial, type = type)*v.eff*dx
      }

      for (i in seq(n.sources)) {
        rho.observed = 1/sqrt(2*pi)/x.sigma[i]*exp(-(x.obs[i]-x)^2/2/x.sigma[i]^2)
        rho.corrected = rho.observed*rho.initial
        rho.corrected = rho.corrected/(sum(rho.corrected)*dx)
        rho.unbiased = rho.unbiased+rho.corrected
      }
    }

    # make -ln(L)
    neglogL = function(p) {
      phi = mfmodel(x,p,type=type)
      list = phi>0 & v.eff>0
      return(sum(phi[list]*v.eff[list]-log(phi[list])*rho.unbiased[list])*dx)
    }

    # maximize ln(L)
    opt = optim(p.initial,neglogL,hessian=TRUE,control=(parscale=abs(p.initial)))
    p.initial = opt$par
  }

  # determine uncertainties
  covariance = solve(opt$hessian)
  sigma = sqrt(diag(covariance))

  # make output
  output = list(p.optimal = opt$par,
                sigma = sigma,
                covariance = covariance)
  return(output)

}

.mfresample <- function(mf, seed = 1) {

  set.seed(seed)
  niterations = mf$input$resampling.iterations
  np = mfmodel(type = mf$input$type, output = 'npara')
  p.initial = mf$fit$parameters$p.optimal

  dx = min(0.1,(max(mf$fit$plot$x)-min(mf$fit$plot$x))/100)
  x = seq(min(mf$fit$plot$x),max(mf$fit$plot$x),by=dx)
  v = mf$fit$veff.fn(x)
  density = pmax(0,mfmodel(x, p.initial, type = mf$input$type)*v)
  density = density/sum(density)
  cum = cumsum(density)
  npoints = length(mf$input$x)
  p.new = array(NA,c(niterations,np))
  x.obs = array(NA,npoints)

  for (iteration in seq(niterations)) {
    r = runif(npoints)
    for (i in seq(npoints)) {
      index = which.min(abs(cum-r[i]))
      x.obs[i] = x[index]
    }
    p.new[iteration,] = .corefit(x.obs, x.sigma = NULL,
                                 p.initial = p.initial,
                                 x = x, v.eff = v,
                                 type = mf$input$type,
                                 nit = 1)$p.optimal
  }

  p.quant.05 = array(NA,np)
  p.quant.16 = array(NA,np)
  p.quant.25 = array(NA,np)
  p.quant.75 = array(NA,np)
  p.quant.84 = array(NA,np)
  p.quant.95 = array(NA,np)
  for (i in seq(np)) {
    p.quant.05[i] = quantile(p.new[,i],0.05,names=FALSE)
    p.quant.16[i] = quantile(p.new[,i],0.16,names=FALSE)
    p.quant.25[i] = quantile(p.new[,i],0.25,names=FALSE)
    p.quant.75[i] = quantile(p.new[,i],0.75,names=FALSE)
    p.quant.84[i] = quantile(p.new[,i],0.84,names=FALSE)
    p.quant.95[i] = quantile(p.new[,i],0.95,names=FALSE)
  }
  p.quantile = list(p.quantile.05 = p.quant.05,
                    p.quantile.16 = p.quant.16,
                    p.quantile.25 = p.quant.25,
                    p.quantile.75 = p.quant.75,
                    p.quantile.84 = p.quant.84,
                    p.quantile.95 = p.quant.95)

  x = mf$fit$plot$x
  s = array(NA,c(niterations,length(x)))
  for (iteration in seq(niterations)) {
    s[iteration,] = mfmodel(x,p.new[iteration,],type=mf$input$type)
  }
  phi.quant = array(NA,c(6,length(x)))
  for (i in seq(length(x))) {
    list = !is.na(s[,i]) & is.finite(s[,i]) & (s[,i]>0)
    phi.quant[,i] = quantile(s[list,i],c(0.05,0.16,0.25,0.75,0.84,0.95),names=FALSE)
  }

  phi.quantile = list(phi.quantile.05 = phi.quant[1,],
                      phi.quantile.16 = phi.quant[2,],
                      phi.quantile.25 = phi.quant[3,],
                      phi.quantile.75 = phi.quant[4,],
                      phi.quantile.84 = phi.quant[5,],
                      phi.quantile.95 = phi.quant[6,])

  invisible(list(phi.quantile = phi.quantile, p.quantile = p.quantile))

}
