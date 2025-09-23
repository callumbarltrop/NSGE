
rm(list=ls())
library(SphericalCubature)
library(copula)
library(evgam)
library(qgam)
library(evd)
library(numDeriv)

Laplace_inverse = function(u){ #standard Laplace quantile function
  x = c()
  x[u<=0.5] = log(2*u[u<=0.5])
  x[u>0.5] = -log(2*(1-u[u>0.5]))
  return(x)
}

Laplace_cdf = function(x){ #Standard Laplace cumulative distribution function
  u = c()
  u[x<0] = exp(x[x<0])/2
  u[x>=0] = 1-exp(-x[x>=0])/2
  return(u)
}

Laplace_density<-function(z){ #Standard Laplace density function
  nz<-length(z)
  result<-numeric(nz)
  result[z<0]<- 0.5*exp(z)
  result[z>=0]<- 0.5*exp(-z)
  return(result)
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

adjustment_func = function(x){
  upper = x[1]
  lower = x[2]
  
  x = x[3:length(x)]
  
  indices = x > 0
  
  adjusted_est = c()
  adjusted_est[which(indices == T)] = x[which(indices == T)]/upper
  adjusted_est[which(indices == F)] = x[which(indices == F)]/lower
  return(adjusted_est)
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

trun_gamma_d0 <- function(pars, likdata) {

  # this is needed, and should be fine
  pars <- split(pars, likdata$idpars)

  # likdata$X should be set up by evgam
  log_gauge <- likdata$X[[1]] %*% pars[[1]]

  nllh = sum(mapply(truncated_gamma_nll,log_gauge,polar_df_exc$r,polar_df_exc$r_thresh))

  return(nllh)
}

trun_gamma_d12 <- function(pars, likdata, sandwich = TRUE) {

  # this is needed, and should be fine
  pars <- split(pars, likdata$idpars)

  # likdata$X should be set up by evgam
  log_gauge <- likdata$X[[1]] %*% pars[[1]]

  deriv_matrix = matrix(NA,ncol=2,nrow=length(log_gauge))

  for(i in 1:length(log_gauge)){
    point = log_gauge[i]

    grad_result <- grad(truncated_gamma_nll,point,r=polar_df_exc$r[i],r_thresh = polar_df_exc$r_thresh[i])

    hessian_result <- hessian(truncated_gamma_nll, point,r=polar_df_exc$r[i],r_thresh = polar_df_exc$r_thresh[i])

    deriv_matrix[i,] = c(grad_result,c(hessian_result))
  }

  return(deriv_matrix)
}

trun_gamma_fns <- list(d0 = trun_gamma_d0, d120 = trun_gamma_d12, d340 = NULL)

unlink <- list(function(x){exp(x)})
attr(unlink[[1]], "deriv") <- unlink[[1]]
#tgamma_quantile <- function(p){ -log(1-p) } #ignore this part

#trun_gamma_fns$q <- tgamma_quantile
trun_gamma_fns$unlink <- unlink
