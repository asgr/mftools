#' Fit a galaxy mass function
#'
#' This function fits an analytical model of a galaxy mass function (MF) to a discrete set of galaxies with observed masses and effective volumes.
#'
#' @param x Normally a vector of log-masses \code{x = log10(M/Msun)}, but also be used with log-luminosities or absolute magnitudes.
#' @param vmax This is either a function or a vector. If vmax is a function, it will be interpreted as the analytical function \code{Vmax(x)} representing the effective cosmic volume, in which the log-masses \code{x} can be detected; normally this is an increasing function of \code{x}. If \code{vmax} is a vector of same length as \code{x}, it will be interpreted as the values of the effective cosmic volumes associated with the different values of \code{x}, and a function \code{Vmax(x)} will be fitted automatically. In both cases the continuous function \code{Vmax(x)} will be returned as \code{vmax.fn} in the output list.
#' @param xerr Optional vector of same length as \code{x}, containing the 1-sigma uncertainties of \code{x}.
#' @param type Kind of MF: \code{'Schechter'} for Schechter function or \code{'PL'} for a simple power law.
#' @param resampling If true, the MF is resampled resampling.iterations times to determine quantiles of the parameter distributions and fitted functions.
#' @param resampling.iterations Integer (>0) number of iterations for the resampling. Ignored if \code{resampling.method = 'none'}.
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
#' mf = mffit(data$x, data$vmax, data$xerr, resampling = TRUE)
#' mfplot(mf, bin.xmin = 6.5, bin.xmax = 9.5, xlim=c(2e6,5e10), ylim=c(1e-3,1))
#'
#' @author Danail Obreschkow, 2017
#'
#' @export

mffit <- function(x,
                  vmax, # either a vector of Vmax for all points x or a function Vmax(x)
                  xerr = NULL,
                  type = 'Schechter',
                  resampling = FALSE,
                  resampling.iterations = 100,
                  write.fit = TRUE) {

  # Input handling
  if (length(x)<=2) stop('give at least 2 data points')
  npoints = length(x)
  wx = max(x)-min(x)

  # make x-range
  if (!is.null(xerr)) {
    if (length(xerr)!=npoints) stop('Length of xerr must be the same as the length of x.')
  }
  if (is.null(xerr)) {
    xmin = min(x)
    xmax = max(x)+wx
  } else {
    xmin = min(x-xerr*2)
    xmax = max(x+xerr*2)+wx
  }
  dx = min(0.01,(xmax-xmin)/500)
  xrange = seq(xmin+dx/2,xmax,by=dx)

  # Make Vmax function
  if (is.function(vmax)) {
    vmax.fn = vmax
  } else {
    if (length(vmax)!=npoints) stop('Length of vmax must be the same as the length of x.')
    if (min(vmax)<=0) stop('All values of vmax must be positive.')
    # fit ln(vmax-a[3]) = a[2]*(x-a[1])
    list = vmax<max(vmax)
    if (is.null(xerr)) {
      minfct <- function(a) {
        return(sum((exp(a[2]*(x[list]-a[1]))+a[3]-vmax[list])^2))
      }
    } else {
      minfct <- function(a) {
        return(sum((exp(a[2]*(x[list]-a[1]))+a[3]-vmax[list])^2/xerr[list]))
      }
    }
    a = optim(c(9,10,0),minfct)$par
    vmax.fn <- function(xval) {
      v = pmax(0,pmin(max(vmax),exp(a[2]*(xval-a[1]))+a[3]))
      v[xval<xmin] = 0
      return(v)
    }
  }

  # Make source density function
  phi_data = xrange*0
  if (is.null(xerr)) {
    for (i in seq(npoints)) {
      index = which.min(abs(x[i]-xrange))[1]
      phi_data[index] = phi_data[index]+1/dx
    }
  } else {
    list = vmax.fn(xrange)<=0
    for (i in seq(npoints)) {
      phi_source = 1/sqrt(2*pi)/xerr[i]*exp(-(x[i]-xrange)^2/2/xerr[i]^2)
      phi_source[list] = 0
      phi_source = phi_source/(sum(phi_source)*dx)
      phi_data = phi_data+phi_source
    }
  }

  # Find most likely model
  v = vmax.fn(xrange)
  neglogL = function(p) {
    m = mfmodel(xrange,p,type=type)
    list = m>0 & v>0
    return(sum(m[list]*v[list]-log(m[list])*phi_data[list]))
  }
  p.initial = mfmodel(output = 'default', type = type)
  opt = optim(p.initial,neglogL,hessian=TRUE,control=(parscale=abs(p.initial)))
  hessian = opt$hessian*dx
  fisher.matrix = solve(hessian)
  sigma = sqrt(diag(fisher.matrix))

  # Determine 1-sigma range of model function
  eig = eigen(fisher.matrix)
  np = length(p.initial)
  index = 0
  phi.new = array(NA,c(length(xrange),np*2))
  for (i in seq(np)) {
    for (k in c(-1,1)) {
      index = index+1
      v = k*sqrt(eig$values[i])*eig$vectors[,i]
      p.new = opt$par+v
      phi.new[,index] = mfmodel(x = xrange, p = p.new, type = type)
    }
  }
  phi.sigma = array(NA,length(xrange))
  for (i in seq(length(xrange))) {
    phi.sigma[i] = (max(phi.new[i,])-min(phi.new[i,]))/2
  }

  # Initialize output parameters
  input = list(x = x, xerr = xerr, vmax = vmax, type = type,
               resampling = resampling, resampling.iterations = resampling.iterations)
  fit = list(parameters = list(equation = mfmodel(output='equation', type=type),p.optimal = opt$par, p.sigma = sigma, p.cov = fisher.matrix),
             plot = list(x = xrange, phi.optimal = mfmodel(x = xrange, p = opt$par, type = type),
                         phi.sigma = phi.sigma),
             vmax.fn = vmax.fn)
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

.mfresample = function(mf, seed = 2) {

  set.seed(seed)
  niterations = mf$input$resampling.iterations
  np = mfmodel(type = mf$input$type, output = 'npara')

  dx = min(0.1,(max(mf$fit$plot$x)-min(mf$fit$plot$x))/100)
  x = seq(min(mf$fit$plot$x),max(mf$fit$plot$x),by=dx)
  density = pmax(0,mfmodel(x, mf$fit$parameters$p.optimal, type = mf$input$type)*mf$fit$vmax.fn(x))
  density = density/sum(density)
  cum = cumsum(density)
  npoints = length(mf$input$x)
  badlist = density==0
  p.new = array(NA,c(niterations,np))

  for (iteration in seq(niterations)) {
    phi_data = x*0
    r = runif(npoints)
    for (i in seq(npoints)) {
      index = which.min(abs(cum-r[i]))
      phi_data[index] = phi_data[index]+1/dx
    }
    phi_data[badlist] = 0 # just to be sure
    v = mf$fit$vmax.fn(x)
    neglogL = function(p) {
      m = mfmodel(x,p,type=mf$input$type)
      list = m>0 & v>0
      return(sum(m[list]*v[list]-log(m[list])*phi_data[list]))
    }
    p.new[iteration,] = optim(mf$fit$parameters$p.optimal,neglogL,control=list(parscale=abs(mf$fit$parameters$p.optimal)))$par
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
