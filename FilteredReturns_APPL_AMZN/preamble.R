#This file contains a range of functions and packages required for fitting the modelling framework 

library(evgam)

packages = c("SphericalCubature","qgam","circular","rugarch","tseries") 
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

Laplace_inverse = function(u){ 
  x = c()
  x[u<=0.5] = log(2*u[u<=0.5])
  x[u>0.5] = -log(2*(1-u[u>0.5]))
  return(x)
}

Laplace_cdf = function(x){ 
  u = c()
  u[x<0] = exp(x[x<0])/2
  u[x>=0] = 1-exp(-x[x>=0])/2
  return(u)
}

dLaplace = function(x, location = 0, logscale = log(1),log=F) {
  scale = exp(logscale)
  if(log){
    return(  log((1 / (2 * scale)) * exp(-abs(x - location) / scale)) )
  } else {
    return(  (1 / (2 * scale)) * exp(-abs(x - location) / scale) )
  }

}

Laplace_density<-function(z){ #Standard Laplace density function
  nz<-length(z)
  result<-numeric(nz)
  result[z<0]<- 0.5*exp(z[z<0])
  result[z>=0]<- 0.5*exp(-z[z>=0])
  return(result)
}

Laplace_nll = function(par,dat){
  location = par[1]
  logscale = par[2]
  
  nll = -sum(dLaplace(dat,location=location,logscale = logscale,log=T))
  
  return(nll)
  
}

inverse_semi_empirical_cdf = function(u,alpha,gpd_paras,marg_vec){ 
  
  nvec = c()
  
  nvec[u<=1-alpha & alpha <= u] = quantile(marg_vec,u[u<=1-alpha & alpha <= u]) #for probabilities below gpd quantile, we use empirical distribution to do inverse transformation
  
  if(sum(u>1-alpha)>0){ #probabilities above some high threshold. Had to add if statement as qgpd doesn't like it when you input an empty vector
    nvec[u>1-alpha] = evd::qgpd((u[u>1-alpha]-(1-alpha))/(alpha),loc=gpd_paras[1,1], scale = gpd_paras[2,1], shape = gpd_paras[3,1]) # had issue here, since all value simulated were below threshold
  }
  
  if(sum(u<alpha)>0){ #probabilities above some high threshold. Had to add if statement as qgpd doesn't like it when you input an empty vector
    nvec[u<alpha] = -evd::qgpd(1-(u[u<alpha]/alpha),loc=-gpd_paras[1,2], scale = gpd_paras[2,2], shape = gpd_paras[3,2]) # had issue here, since all value simulated were below threshold
  }
  
  return(nvec)
}

l1_norm = function(x){
  return(sum(abs(x)))
}

l2_norm = function(x){
  return(sqrt(sum(x^2)))
}

linf_norm = function(x){
  return(max(abs(x)))
}

gauge_normal<-function(wpts,rho){
  x<- wpts[1];y<-wpts[2]
  if(x < 0 & y >= 0 | x >= 0 & y < 0){
    return((abs(x)+abs(y)+2*rho*sqrt(abs(x*y)))/(1-rho^2) )
  } else {
    return((abs(x)+abs(y)-2*rho*sqrt(abs(x*y)))/(1-rho^2) )
  }
}

gauge_inv_logistic = function(x,dep){
  d = 2
  if(sum(sign(x)) == -length(x)){
    return( (1/dep)*sum(-x) + (1-(d/dep))*min(abs(x)) )
  } else if(sum(sign(x)) == length(x)){
    return( ( sum( (x)^(1/dep)  ) )^(dep) )
  } else {
    return( (1/dep)*sum(-x[sign(x)==-1]) +( sum( (x[sign(x)==1])^(1/dep)  ) )^(dep) )
  }
  
}

gauge_mv_t = function(x,nu){
  d = 2
  return( -(1/nu)*sum(abs(x)) + (1 + (d/nu))*max(abs(x)) )
}

gauge_frank <- function(w){
  x<- w[1];y<-w[2]
  return(abs(x) + abs(y))
}

gauge_joe <- function(wpts,alpha){
  x<- wpts[1];y<-wpts[2]
  w1<- x/(abs(x)+abs(y));w2<- y/(abs(x)+abs(y))
  if(w1 < 0 & w2 < 0){
    return(abs(x)+abs(y))
  } else if(w1 >= 0 & w2 < 0){
    return( (abs(x)+abs(y)) *(1 + (alpha-1)*abs(w1) ))
  } else if(w1>=0 & w2>=0){
    return((abs(x)+abs(y)) *(1 - alpha + (2*alpha - 1)*max(w1,w2)))
  } else {
    return((abs(x)+abs(y)) *(1 + (alpha-1)*abs(w2)))
  }
}

gauge_hw<-function(wpts,par,ll=-10,ul=10)
{
  x<-wpts[1];y<-wpts[2]
  gamma<-par[1]
  dummy<-function(q)
  {
    gauge_normal(wpts=(wpts-gamma*q),rho=par[2])+abs(q)
  }
  opt<-optimize(dummy,lower = ll,upper = ul)
  out<-opt$obj
  
  # Numerical optimization doesn't always work correctly
  # These two lines add to the time for calculating the gauge, but help if there are optimization issues
  #===========================================================
  check<-sapply(seq(ll,ul,len=500),dummy)
  if(any(check<opt$obj)){out<-min(check)}
  #===========================================================
  
  # Rescale if gamma>1:
  if(gamma>1){out<-out*gamma}
  return(out)
}

sim_data_normal = function(rho){
  d = 2
  normc = normalCopula(param = rho, dim = d)
  return(apply(rCopula(1, copula = normc),2,Laplace_inverse))
}

sim_data_invlog = function(alpha){
  return(Laplace_inverse(pexp(1/rbvevd(n=1,dep=alpha,model="log",mar1 = c(1,1,1)))))
}

sim_data_t = function(nu,rho){
  tc = tCopula(param = rho,df=nu,dim = 2)
  return(apply(rCopula(1, copula = tc),2,Laplace_inverse))
}

sim_data_frank = function(theta){
  d = 2
  fc = frankCopula(theta, dim = d)
  return(apply(rCopula(1, copula = fc),2,Laplace_inverse))
}

sim_data_joe = function(alpha){
  d = 2
  jc = joeCopula(alpha, dim = d)
  return(apply(rCopula(1, copula = jc),2,Laplace_inverse))
}

FX<-function(x,gamma)
{
  if(gamma<=1)
  {
    to.int<-function(s)
    {
      Laplace_cdf(x-gamma*s)*Laplace_density(s)
    }
    cdf<-integrate(to.int,lower=-Inf,upper=Inf)$value
  }
  else if(gamma>1)
  {
    to.int<-function(v)
    {
      Laplace_cdf(x-v/gamma)*Laplace_density(v)
    }
    cdf<-integrate(to.int,lower=-Inf,upper=Inf)$value
  }
  return(cdf)
}

sim_data_hw = function(gamma,rho){
  
  V = sim_data_normal(rho=rho)
  S = Laplace_inverse(runif(1))
  
  if(gamma>1){
    X<-S+V/gamma
  } else {
    X<-gamma*S+V
  }

  FX.vec<-Vectorize(FX,vectorize.args = "x")
  
  Xlap <- Laplace_inverse(FX.vec(x=X,gamma=gamma))
  
  return(Xlap)
}

truncated_gamma_nll = function(log_gauge,r,r_thresh){
  alpha = 2
  gauge = exp(log_gauge)
  ll = alpha*log(gauge) + (alpha-1)*log(r) - gauge*r - lgamma(alpha) - log(pgamma(r_thresh,shape = alpha,rate=gauge,lower.tail = FALSE))
  return(-ll)
}


#=====================================================================
# Functions for Generalised Pareto Distribution.
# Added option to use nu parameterisation
# checks that param values are valid
#=====================================================================
# pgpd
# qgpd
# dgpd
# rgpd

# GPD_LL
#=====================================================================
#' Generalised Pareto Distribution
#'
#' Cumulative density function of the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0 or p<mu.
#'
#' @author Zak Varty
#'
#' @param q vector of quantiles.
#' @param shape shape parameter (xi)
#' @param scale scale parameter (sigma)
#' @param nu  alternative scale parameter: nu = sigma/(1+xi)
#' @param mu  location parameter
#' @param skip_checks logical. Speed up evaluation by skipping checks on inputs? (Beware!)
#' @return Probability of the GPD X<=q
#' @importFrom stats pexp
#' @examples
#' pgpd(q = c(-1,1.5,3), shape = 1, scale = 1)
#' pgpd(q = 1.5, shape = c(0,-1), scale = c(0.1,1))
#' @export

pgpd <- function(q, shape, scale = NULL, nu = NULL, mu = 0, skip_checks = FALSE){
  
  if (!skip_checks) {
    # one and only one of {nu, scale} may be specified
    if (is.null(scale) & is.null(nu)) {
      stop('Define one of the parameters nu or scale.')
    }
    if (!is.null(scale) & !is.null(nu)) {
      stop('Define only one of the parameters nu and scale.')
    }
    # Calculate scale from nu if required
    if (!is.null(nu) & is.null(scale)) {
      scale <- nu / (1 + shape)
      if (any(scale <= 0)) {
        stop('Implied scale parameter(s) must be positive.')
      }
      
    }
    # Check that scale value(s) are positive
    if (any(scale <= 0)) {
      stop('Scale parameter(s) must be positive.')
    }
    
    # Ensure q, scale, shape and mu are of same length.
    if (length(scale) == 1 & length(q) > 1) {
      scale <- rep(scale, length(q))
    }
    if (length(shape) == 1 & length(q) > 1) {
      shape <- rep(shape, length(q))
    }
    if (length(mu) == 1 & length(q) > 1) {
      mu <- rep(mu, length(q))
    }
  } else {
    if (!is.null(nu) & is.null(scale)) {
      scale <- nu / (1 + shape)
    }
  }
  #calculate probabilities
  p <- (1 - (1 + (shape * (q - mu))/scale)^(-1/shape))
  #correct probabilities below mu or above upper end point
  p[q < mu] <- 0
  p[(shape < 0) & (q >= (mu - scale/shape))] <- 1
  
  #correct probabilities where xi = 0
  if (any(abs(shape) < 1e-10)) {
    #ex <- which(shape ==0)
    ex <- which(abs(shape) < 1e-10)
    p[ex] <- pexp(q = q[ex] - mu[ex], rate = 1 / scale[ex])
  }
  
  return(p)
}

#' Generalised Pareto Distribution
#'
#' Cumulative density function of the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0 or p is not a valid
#' probability.
#'
#' @author Zak Varty
#'
#' @param p vector of quantiles.
#' @param shape shape parameter (xi)
#' @param scale scale parameter (sigma)
#' @param nu  alternative scale parameter: nu = sigma/(1+xi)
#' @param mu  location parameter
#' @return Probability of the GPD X<=x
#' @examples
#' qgpd(p = 0.5, shape = 0.5, scale = 0.5)
#' \dontrun{ qgpd(p = -0.1, shape = 0, scale = 1, mu = 0.1) }
#' @export
qgpd <- function(p, shape, scale = NULL, nu = NULL, mu = 0){
  # one and only one of {nu, scale} may be specified
  if (is.null(scale) & is.null(nu)) {
    stop('Define one of the parameters nu or scale.')
  }
  if (!is.null(scale) & !is.null(nu)) {
    stop('Define only one of the parameters nu and scale.')
  }
  
  # Probabilities must all be positive
  if (!all((p >= 0) & (p <= 1))) {
    stop('Probabilities p must be in the range [0,1].')
  }
  
  # Calculate scale from nu if required
  if (!is.null(nu) & is.null(scale)) {
    scale <- nu / (1 + shape)
    if (any(scale <= 0)) {
      stop('Implied scale parameter(s) must be positive.')
    }
    
  }
  
  # Check that scale value(s) are positive
  if (any(scale <= 0)) {
    stop('Scale parameter(s) must be positive.')
  }
  # Ensure p, scale, shape and mu are of same length.
  if (length(scale) == 1 & length(p) > 1) {
    scale <- rep(scale, length(p))
  }
  if (length(shape) == 1 & length(p) > 1) {
    shape <- rep(shape, length(p))
  }
  if (length(mu) == 1 & length(p) > 1) {
    mu <- rep(mu, length(p))
  }
  
  #calculate quantiles
  q <- mu + (scale/shape) * ((1 - p)^(-shape) - 1)
  
  #correct quantiles where xi = 0
  #ex <- which(shape ==0)
  if (any(abs(shape) < 1e-10)) {
    ex <- which(abs(shape) < 1e-10)
    q[ex] <- mu[ex] + stats::qexp(p = p[ex],rate = 1/scale[ex])
  }
  return(q)
}

# The following functions are used for fitting the GPD distribution to the margins. They are obtained from the code published at https://github.com/conor-murphy4/automated_threshold_selection 

#' Generalised Pareto Distribution
#'
#' Density function of the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0 or x is
#' outside of the domain of the given distribution.
#'
#' @author Zak Varty
#'
#' @param x vector of values as which to evaluate density.
#' @param shape shape parameter (xi)
#' @param scale scale parameter (sigma)
#' @param nu  alternative scale parameter
#' @param mu  location parameter
#' @param log  locical. Return log
#' @return density of the GPD at x
#' @examples
#' dgpd(x = c(-1,0.5,1,1.9,5),shape = -0.5, scale = 1)
#' @export
#'
dgpd <- function(x, shape, scale = NULL, nu = NULL, mu = 0, log = FALSE){
  # one and only one of {nu, scale} may be specified
  if (is.null(scale) & is.null(nu)) {
    stop('Define one of the parameters nu or scale.')
  }
  if (!is.null(scale) & !is.null(nu)) {
    stop('Define only one of the parameters nu and scale.')
  }
  
  # Calculate scale from nu if required
  if (!is.null(nu) & is.null(scale)) {
    scale <- nu / (1 + shape)
    if (any(scale <= 0)) {
      stop('Implied scale parameter(s) must be positive.')
    }
  }
  
  # Check that scale value(s) are positive
  if (any(scale <= 0)) {
    stop('Scale parameter(s) must be positive.')
  }
  # Ensure x, scale, shape and mu are of same length.
  if (length(scale) == 1 & length(x) > 1) {
    scale <- rep(scale, length(x))
  }
  if (length(shape) == 1 & length(x) > 1) {
    shape <- rep(shape, length(x))
  }
  if (length(mu) == 1 & length(x) > 1) {
    mu <- rep(mu, length(x))
  }
  
  if (log == FALSE) {
    out <- (scale^(-1)) * pmax((1 + shape * (x - mu)/scale),0)^((-1/shape) - 1)
    # amend values below threshold
    out[which(x < mu)] <- 0
    # amend values above upper endpoint (if it exists)
    out[which((shape < 0) & (x >= (mu - scale/shape)))] <- 0
    # amend values where xi = 0 (if they exist)
    if (any(abs(shape < 1e-10))) {
      ex <- which(abs(shape) < 1e-10)
      out[ex] <- stats::dexp(x = x[ex] - mu[ex], rate = 1/scale[ex])
    }
  } else {
    out <-  -log(scale) + ((-1/shape) - 1)*log(pmax((1 + shape * (x - mu)/scale),0))
    # amend values below threshold
    out[which(x < mu)] <- -Inf
    # amend values above upper endpoint (if it exists)
    out[which((shape < 0) & (x >= (mu - scale/shape)))] <- -Inf
    # amend values where xi = 0 (if they exist)
    if (any(abs(shape) < 1e-10)) {
      ex <- which(abs(shape) < 1e-10)
      out[ex] <- stats::dexp(x = x[ex] - mu[ex], rate = 1 / scale[ex],log = TRUE)
    }
  }
  return(out)
}

#' Generalised Pareto Distribution
#'
#' Sample the GPD specified in terms of (sigma,xi) or (nu,xi).
#' Improvement over evir function as it returns an error if the (implied) shape
#' parameter is non-positive. Also properly handles cases where and xi=0.
#'
#' @author Zak Varty
#'
#' @param n sample size.
#' @param shape shape parameter (xi).
#' @param scale scale parameter (sigma).
#' @param nu  alternative scale parameter.
#' @param mu  location parameter.
#' @return Random sample from generalised pareto distirbution.
#'
#' @examples
#' rgpd(n = 100, shape = 0, scale = 1:100)
#' @export
rgpd <- function(n, shape, scale = NULL, nu = NULL, mu = 0){
  ## Input checks
  # one and only one of {nu, scale} may be specified
  if (is.null(scale) & is.null(nu)) {
    stop('Define one of the parameters nu or scale.')
  }
  if (!is.null(scale) & !is.null(nu)) {
    stop('Define only one of the parameters nu and scale.')
  }
  # Calculate scale from nu if required
  if (!is.null(nu) & is.null(scale)) {
    scale <- nu / (1 + shape)
    if (any(scale <= 0)) {
      stop('Implied scale parameter(s) must be positive.')
    }
  }
  # Check that scale value(s) are positive
  if (any(scale <= 0)) {
    stop('Scale parameter(s) must be positive.')
  }
  # Ensure q, scale, shape and mu are of same length.
  if ((length(scale) == 1) & (n > 1)) {
    scale <- rep(scale, n)
  }
  if ((length(shape) == 1) & (n > 1)) {
    shape <- rep(shape, n)
  }
  if ((length(mu) == 1) & (n > 1)) {
    mu <- rep(mu, n)
  }
  
  #simulate sample
  sample <- mu + (scale/shape) * ((1 - stats::runif(n))^(-shape) - 1)
  #correct sample values where xi = 0
  #ex <- which(shape ==0)
  if (any(abs(shape) < 1e-10)) {
    ex <- which(abs(shape) < 1e-10)
    sample[ex] <- mu[ex] +
      stats::rexp(n = length(ex),rate = 1/scale[ex])
  }
  return(sample)
}

#' evalue probability mass function of rounded generalised Pareto distribution
#'
#' @author Zak Varty
#'
#' @param x Vector values at which to evaluate mass function
#' @param u Vector of latent threshold values
#' @param sig_u Vector of latent scale parameters (for exceedances of u)
#' @param xi Latent shape parameter
#' @param to_nearest Level of rounding
#'
#' @return pmf evaluated at x. NOTE: does not check validity of x values.
#'
#' @examples
#' dgpd_rd(x = seq(0.1, 1.5, by = 0.1), to_nearest = 0.1, u = 0, sig_u = 1, xi = 0)
#' dgpd_rd(x = seq(0.1, 1.5, by = 0.1), to_nearest = 0.1, u = seq(0, 1.4, by = 0.1), sig_u = 1, xi = 0)
#' # CAUTION:
#' gpd_rd(x = 0.15, to_nearest = 0.1,  u = 0, sig_u = 1, xi = 0)
dgpd_rd <- function(x, u, sig_u, xi, to_nearest){
  
  # If (Y_i - u_i | Y_i > u_i) ~ GPD(sig_i, xi)
  # then Z_i  = [(Y_i - u_i)/sig_i | Y_i - u_i > 0] ~ GPD(1, xi)
  
  # range of z values that lead to observing x
  x_low <- pmax(x - to_nearest/2, u)
  z_low <- (x_low - u) / sig_u
  z_high <- (x - u + to_nearest/2)/sig_u
  
  # calculate probability of z in that range
  p_high <- pgpd(q = z_high, scale = 1, shape = xi, mu = 0)
  p_low  <- pgpd(q = z_low,  scale = 1, shape = xi, mu = 0)
  p <- p_high - p_low
  
  return(p)
}


#' Generalised Pareto log-likelihood
#'
#' @author Conor Murphy
#'
#' @param par A numeric vector of parameter values of length 2.
#' @param z A numeric vector of excesses of some threshold.
#'
#' @returns A numeric value of the log-likeihood.
#'
#' @examples
#' test1 <- rgpd(1000, shape = 0.1, scale=0.5, mu=1)
#' excess <- test1[test1>1.5] - 1.5
#' GPD_LL(par=c(1,0.4), z=excess)


GPD_LL <- function(par, z){
  sigma <- par[1]
  xi <- par[2]
  if (sigma > 0) {
    if (abs(xi) < 1e-10) {
      return(-length(z) * log(sigma) - ((1 / sigma) * sum(z)))
    }
    else {
      if (all(1 + (xi * z) / sigma > 0)) {
        return(-(length(z) * log(sigma)) - ((1 + 1 / xi)*(sum(log(1 + (xi * z) / sigma)))))
      }
      else{
        return(-1e6)
      }
    }
  }
  else{
    return(-1e7)
  }
}


transform_to_exp <- function (y,sig, xi){
  std_exp <- (1 / xi) * log( 1 + xi * (y/sig))  
  return(std_exp)
}



