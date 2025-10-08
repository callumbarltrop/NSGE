rm(list=ls())

source("preamble.R")

symbols = sort(c("AAPL","AMZN"))

symbols2 = c("Apple","Amazon")[order(c("AAPL","AMZN"))]

# Data accessed on 25.08.2025

returns_datasets = readRDS(file="data/ret_data.rds")

dates = readRDS(file="data/ret_dates.rds")

# Remove any zero log-returns entries 

dates = dates[rowSums(returns_datasets != 0) == ncol(returns_datasets)]

returns_datasets = returns_datasets[rowSums(returns_datasets != 0) == ncol(returns_datasets), ]

# Compute length of data frame 

n_data = as.numeric(lengths(returns_datasets)[1])

# Plot ACF for returns and squared returns 

pdf(file="figures/ret_acfs.pdf",width=10,height=5)

par(mfrow=c(1,2),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

for(s in symbols){
  dat = returns_datasets[[s]]

  acf(dat,main=symbols2[which(symbols==s)],
       cex.lab=1.3, cex.axis=1.2,cex.main=1.8, cex=0.5)

}

dev.off()

pdf(file="figures/vol_acfs.pdf",width=10,height=5)

par(mfrow=c(1,2),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

for(s in symbols){
  dat = returns_datasets[[s]]

  acf(dat^2,main=symbols2[which(symbols==s)],
      cex.lab=1.3, cex.axis=1.2,cex.main=1.8, cex=0.5)

}

dev.off()


# Filtering the returns and testing for stationarity  ---------------------

# Fitting a GARCH(1,1) model to filter the returns

filtered_returns_datasets = list()

garch_paras = list()

garch_spec = ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0)), 
  distribution.model = "norm"             
)

for(s in symbols){
  dat = returns_datasets[[s]]

  garch_fit = ugarchfit(spec=garch_spec,data=dat,xts_check_TZ = FALSE)

  filtered_returns_datasets[[s]] = as.vector(residuals(garch_fit, standardize = TRUE))

  garch_paras[[s]] = data.frame(mu = rep(coef(garch_fit)[1],length(dat)),sigma = as.vector(sigma(garch_fit)),t = 1:length(dat))
}

saveRDS(garch_paras,file = "data/garch_paras.rds")

saveRDS(filtered_returns_datasets,file="data/fil_ret_data.rds")

filtered_returns_datasets = readRDS(file="data/fil_ret_data.rds")

# ACF plots for filtered and squared filtered returns  

pdf(file="figures/fil_ret_acfs.pdf",width=10,height=5)

par(mfrow=c(1,2),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

for(s in symbols){
  dat = filtered_returns_datasets[[s]]
  
  acf(dat,main=symbols2[which(symbols==s)],
      cex.lab=1.3, cex.axis=1.2,cex.main=1.8, cex=0.5)
  
}

dev.off()

pdf(file="figures/fil_vol_acfs.pdf",width=10,height=5)

par(mfrow=c(1,2),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

for(s in symbols){
  dat = filtered_returns_datasets[[s]]
  
  acf(dat^2,main=symbols2[which(symbols==s)],
      cex.lab=1.3, cex.axis=1.2,cex.main=1.8, cex=0.5)
  
}

dev.off()

# Plotting time series for filtered returns 

pdf(file="figures/fil_ret_ts.pdf",width=10,height=5)

par(mfrow=c(1,2),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

for(s in symbols){
  plot(dates,filtered_returns_datasets[[s]],main=paste0(symbols2[which(symbols==s)]),type="l", lwd=2, ylab = "Filtered log-returns", xlab = "Date",cex.lab=1.5, cex.axis=1.5,cex.main=1.7)
}

dev.off()

# Testing for stationarity in the filtered time series' 

alpha = 0.05

adf_tests = list()

for(s in symbols){
  dat = filtered_returns_datasets[[s]]

  adf_test = adf.test(dat,alternative = "stationary")

  if(adf_test$p.value < alpha){
    adf_tests[[s]] = "Indicates Stationarity"
  } else {
    adf_tests[[s]] = "Indicates Non-stationarity"
  }
}

adf_tests = as.data.frame(adf_tests)

names(adf_tests) = symbols2

print(adf_tests)

kpss_tests = list()

for(s in symbols){
  dat = filtered_returns_datasets[[s]]

  kpss_test = kpss.test(dat, null = "Level")

  if(kpss_test$p.value < alpha){
    kpss_tests[[s]] = "Indicates Non-stationarity"
  } else {
    kpss_tests[[s]] = "Indicates Stationarity"
  }
}

kpss_tests = as.data.frame(kpss_tests)

names(kpss_tests) = symbols2

print(kpss_tests)

# Marginal modelling ------------------------------------------------------

alpha = 0.03 # This corresponds to the threshold level in the GPD fits 

nboot = 200 # Number of bootstraps for parametric bootstrap 

xis = matrix(NA,nrow=2,ncol=length(symbols))

sigmas = matrix(NA,nrow=2,ncol=length(symbols))

us = matrix(NA,nrow=2,ncol=length(symbols))

# Fitting GPD model to upper tails 

pdf(file="figures/gpd_qq_plots_ret_upper.pdf",width=12,height=4)

par(mfrow=c(1,3),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

for(s in symbols){
  dat = filtered_returns_datasets[[s]]

  thresh = quantile(dat,1-alpha)

  us[1,which(s == symbols)] = thresh

  data_over_u = dat[dat > thresh]-thresh

  m = length(data_over_u)

  observed_quants = sort(data_over_u) + thresh

  mle0 = mean(data_over_u)
  init.fit = optim(GPD_LL, z = data_over_u, par = c(mle0,0.1), control = list(fnscale = -1))
  xi = init.fit$par[[2]]
  sigma = init.fit$par[[1]]

  xis[1,which(s == symbols)] = xi
  sigmas[1,which(s == symbols)] = sigma

  theoretical_quants = qgpd((1:m)/(m+1),shape=xi,scale = sigma) + thresh

  boot_theoretical_quants = matrix(0,nrow=nboot,ncol=m)

  for(nb in 1:nboot){

    set.seed(nb)

    boot_data_over_u = data_over_u[sample(1:m,m,replace = T)]

    mle0 = mean(boot_data_over_u)
    boot.init.fit = optim(GPD_LL, z = boot_data_over_u, par = c(mle0,0.1), control = list(fnscale = -1))

    boot_theoretical_quants[nb,] = qgpd((1:m)/(m+1),shape=boot.init.fit$par[[2]],scale = boot.init.fit$par[[1]]) + thresh

  }

  plot(theoretical_quants,observed_quants,xlim=range(theoretical_quants,apply(boot_theoretical_quants,2,quantile,probs = c(0.025,0.975)),observed_quants),
       ylim=range(theoretical_quants,apply(boot_theoretical_quants,2,quantile,probs = c(0.025,0.975)),observed_quants),pch=16,col=1,ylab="Observed",xlab="Theoretical",main=paste0(symbols2[which(symbols==s)]),
       cex.lab=1.3, cex.axis=1.2,cex.main=1.8, cex=0.5)
  polygon(x=c(apply(boot_theoretical_quants,2,quantile,probs = 0.025), rev(apply(boot_theoretical_quants,2,quantile,probs = 0.975))  ),y = c(observed_quants, rev(observed_quants)),col = 'grey80', border = NA)
  abline(a=0,b=1,lwd=3,col=2)
  points(theoretical_quants,observed_quants,pch=16,col="black", cex=1.3)

}

plot(1, 1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", axes = FALSE)

text(0.5, 0.5, "Upper tail", cex = 3)

dev.off()

pdf(file="figures/gpd_qq_plots_ret_lower.pdf",width=12,height=4)

par(mfrow=c(1,3),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

for(s in symbols){
  dat = filtered_returns_datasets[[s]]

  thresh = quantile(dat,alpha)

  us[2,which(s == symbols)] = thresh

  data_over_u = - (dat[dat < thresh]-thresh)

  m = length(data_over_u)

  observed_quants = thresh - sort(data_over_u)

  mle0 = mean(data_over_u)
  init.fit = optim(GPD_LL, z = data_over_u, par = c(mle0,0.1), control = list(fnscale = -1))
  xi = init.fit$par[[2]]
  sigma = init.fit$par[[1]]

  xis[2,which(s == symbols)] = xi
  sigmas[2,which(s == symbols)] = sigma

  theoretical_quants = thresh - qgpd((1:m)/(m+1),shape=xi,scale = sigma)

  boot_theoretical_quants = matrix(0,nrow=nboot,ncol=m)

  for(nb in 1:nboot){

    set.seed(nb)

    boot_data_over_u = data_over_u[sample(1:m,m,replace = T)]

    mle0 = mean(boot_data_over_u)
    boot.init.fit = optim(GPD_LL, z = boot_data_over_u, par = c(mle0,0.1), control = list(fnscale = -1))

    boot_theoretical_quants[nb,] = thresh - qgpd((1:m)/(m+1),shape=boot.init.fit$par[[2]],scale = boot.init.fit$par[[1]])

  }

  plot(theoretical_quants,observed_quants,xlim=range(theoretical_quants,apply(boot_theoretical_quants,2,quantile,probs = c(0.025,0.975)),observed_quants),
       ylim=range(theoretical_quants,apply(boot_theoretical_quants,2,quantile,probs = c(0.025,0.975)),observed_quants),pch=16,col=1,ylab="Observed",xlab="Theoretical",main=paste0(symbols2[which(symbols==s)]),
       cex.lab=1.3, cex.axis=1.2,cex.main=1.8, cex=0.5)
  polygon(x=c(apply(boot_theoretical_quants,2,quantile,probs = 0.025), rev(apply(boot_theoretical_quants,2,quantile,probs = 0.975))  ),y = c(observed_quants, rev(observed_quants)),col = 'grey80', border = NA)
  abline(a=0,b=1,lwd=3,col=2)
  points(theoretical_quants,observed_quants,pch=16,col="black", cex=1.3)

}

plot(1, 1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", axes = FALSE)

text(0.5, 0.5, "Lower tail", cex = 3)

dev.off()

# Save the estimated GPD parameters for the upper and lower tails 

saveRDS(list(thresh = us,scale = sigmas,shape=xis),file="data/gpd_marg_paras.rds")

# Transform the data to the standard Laplace scale 

laplace_datasets = matrix(NA,ncol=length(symbols),nrow=n_data)

for(s in symbols){
  unif = rep(NA,n_data)

  dat = filtered_returns_datasets[[s]]

  u_upper = us[1,which(s == symbols)]

  xi_upper = xis[1,which(s == symbols)]
  sigma_upper = sigmas[1,which(s == symbols)]

  unif[dat > u_upper] = 1-(alpha)*(1-pgpd(dat[dat>u_upper],mu=u_upper,shape=xi_upper,scale = sigma_upper))

  u_lower = us[2,which(s == symbols)]

  xi_lower = xis[2,which(s == symbols)]
  sigma_lower = sigmas[2,which(s == symbols)]

  unif[dat < u_lower] = (alpha)*(1-pgpd(-dat[dat<u_lower],mu=-u_lower,shape=xi_lower,scale = sigma_lower))

  unif[is.na(unif)] = ((rank(dat,ties.method = "random"))/(n_data + 1))[is.na(unif)]

  lap_data = Laplace_inverse(unif)

  laplace_datasets[,which(symbols == s)] = lap_data

}

laplace_datasets = as.data.frame(laplace_datasets)

names(laplace_datasets) = symbols

saveRDS(laplace_datasets,file="data/ret_data_laplace.rds")

# We check the margins have been sucessfully standardised by computing MLE Laplace parameter estimates across rolling windows 

n_window = 20 # Number of rolling windows 

window_lap_para_func_boot = function(vec,dat,nboot = 100){ # Wrapper for computing bootstrapped parameter estimates for each rolling window

  dat = dat[vec]

  par_mat = matrix(NA,ncol=2,nrow=nboot+1)

  est_par = optim(par=c(0,0),fn=Laplace_nll,dat=dat,method="BFGS")

  par_mat[1,] = c(est_par$par[1],exp(est_par$par[2]))

  par_mat[-1,] = t(sapply(1:nboot,function(i){
    set.seed(i)
    bdat = dat[sample(1:length(dat),length(dat),replace = T)]
    est_par = optim(par=c(0,0),fn=Laplace_nll,dat=bdat,method="BFGS")
    return(c(est_par$par[1],exp(est_par$par[2])))
  }))

  est_mu = par_mat[1,1]

  est_sig = par_mat[1,2]

  ci_mu = quantile(par_mat[,1],probs=c(0.025,0.975))
  ci_sig = quantile(par_mat[,2],probs=c(0.025,0.975))

  return(c(est_mu,ci_mu,est_sig,ci_sig))

}

pdf(file="figures/ret_laplace_para_estimates_boot.pdf",width=12,height=6)

par(mfrow=c(1,2),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

for(s in symbols){

  dat = laplace_datasets[[s]]

  vecs = sapply(1:n_window,function(i){
    if(i == n_window){
      return( ((i-1)*round(n_data/n_window)+1):n_data )
    } else {
      return( ((i-1)*round(n_data/n_window)+1):(i*round(n_data/n_window)) )
    }
  })

  mid_points = lapply(vecs,function(x){return(round(median(x)))})

  if(!file.exists(paste0("data/laplace_marg_paras_",s,".rds"))){

    window_paras = lapply(vecs,window_lap_para_func_boot,dat=dat)

    saveRDS(window_paras,file=paste0("data/laplace_marg_paras_",s,".rds"))
  }

  window_paras = readRDS(file=paste0("data/laplace_marg_paras_",s,".rds"))

  mu = c(unlist(lapply(window_paras,function(x){return(x[1])})))

  mu_l = c(unlist(lapply(window_paras,function(x){return(x[2])})))

  mu_u = c(unlist(lapply(window_paras,function(x){return(x[3])})))

  sigma = c(unlist(lapply(window_paras,function(x){return(x[4])})))

  sigma_l = c(unlist(lapply(window_paras,function(x){return(x[5])})))

  sigma_u = c(unlist(lapply(window_paras,function(x){return(x[6])})))

  plot(dates, unlist(sapply(1:n_window,function(i){return(rep(mu[i],length(vecs[[i]])))})) ,ylim=c(-0.8,1.6),lwd=3,col="grey"
       ,ylab="Parameter estimates",xlab="Date",type="l",
       cex.lab=1.3, cex.axis=1.2,cex.main=1.6, cex=0.5,main=paste0("Laplace parameter estimates ",symbols2[which(symbols==s)]))
  lines(dates,unlist(sapply(1:n_window,function(i){return(rep(mu_l[i],length(vecs[[i]])))})),lwd=3,lty=2,col="grey")
  lines(dates,unlist(sapply(1:n_window,function(i){return(rep(mu_u[i],length(vecs[[i]])))})),lwd=3,lty=2,col="grey")
  abline(h = 0,lwd=3,col=1)
  lines(dates,unlist(sapply(1:n_window,function(i){return(rep(sigma[i],length(vecs[[i]])))})),lwd=3,lty=1,col="red4")
  lines(dates,unlist(sapply(1:n_window,function(i){return(rep(sigma_l[i],length(vecs[[i]])))})),lwd=3,lty=2,col="red4")
  lines(dates,unlist(sapply(1:n_window,function(i){return(rep(sigma_u[i],length(vecs[[i]])))})),lwd=3,lty=2,col="red4")
  abline(h=1,lwd=3,col="red")
  if(s == symbols[1]){
    legend("bottomright", legend = c("True location","Estimated location", "True scale","Estimated scale"),
           col = c(1, "grey","red","red4"), lwd=3,cex=1,bg="white")
  }

}

dev.off()

# Fitting the geometric extremes model -----------------------------------------------------------

pred_phis = seq(0,2*pi,length.out=201) # Fine grid of polar angles upon which to evaluate model 

s = symbols[1]

s2 = symbols[2]

data_lap = cbind(laplace_datasets[[s]],laplace_datasets[[s2]]) # Laplace data as a matrix 

d = 2

unit_circle = cbind(cos(pred_phis),sin(pred_phis)) # Unit ball for L2 norm 

knot_angle = 17 # Number of knots for the polar angle

knot_time = 10 # Number of knots for time 

tau = 0.8 # Quantile level for threshold regression 

time = 1:n_data # Discrete time vector 

time_points = round(seq(1,n_data,length.out = 10)) # 10 space time points over the observation period 

num_win = length(time_points)

colfunc = colorRampPalette(c("blue", "orange")) # Defining a colour gradient vector 

grad_cols = colfunc(length(time_points))

polar = rect2polar(t(data_lap)) # Transforming to polar coordinates 

polar_df = data.frame(r=polar$r,phi = polar$phi[1,],logr=log(polar$r),t = time) # Forming a data frame of polar coordinates

fmla_qgam = paste0("logr ~ te(phi,t,bs=c('cc','cr'),k=c(",knot_angle,",",knot_time,"))") # Quantile regression formula for qgam 

phi_knots = seq(0,2*pi,length.out=knot_angle) # Specifying the angular knots 

time_knots = seq(min(time),max(time),length.out=knot_time) # Specifying the time knots 

knots = list(phi = phi_knots,t = time_knots) 

if(!file.exists(paste0("data/qr_model_coeff_",s,"_",s2,".rds"))){ # Fitting the quantile regression model to estimate the threshold function 

  qr_model = qgam(as.formula(fmla_qgam), data=polar_df,  qu=tau,argGam = list(knots = knots) )

  saveRDS(list(coeff = qr_model$coefficients, learn_rate = qr_model$calibr$lsig),file=paste0("data/qr_model_coeff_",s,"_",s2,".rds"))

}

qr_coeff_info = readRDS(file=paste0("data/qr_model_coeff_",s,"_",s2,".rds"))

# Specify a dummy/redundant model. This avoids us having to save the whole qgam object - instead we just save spline coefficients
  
dummy_model = gam(as.formula(paste0("r ~ te(phi,t,bs=c('cc','cr'),k=c(",knot_angle,",",knot_time,"))")),data = polar_df,knots = knots,sp=c(0.1,0.1))

dummy_design = predict(dummy_model, newdata=polar_df,type="lpmatrix") # Compute design matri from dummy model 

thresh_function = exp(c(dummy_design%*%as.matrix(qr_coeff_info$coeff))) # Computing the estimated threshold function

positive_exceedances_indicator = which(polar_df$r - thresh_function>0) # Indices of threshold exceedances  

polar_df_exc = data.frame(r=polar_df$r[positive_exceedances_indicator],phi = (polar$phi[1,])[positive_exceedances_indicator],r_thresh=thresh_function[positive_exceedances_indicator],t=time[positive_exceedances_indicator]) # Data frame of threshold exceeding data 

fmla_trun_gamma = as.formula(paste0("r ~ te(phi,t,bs=c('cc','cr'),k=c(",knot_angle,",",knot_time,"))")) # evgam formula for truncated gamma

if(!file.exists(paste0("data/tg_model_coeff_",s,"_",s2,".rds"))){
  tg_model = evgam(fmla_trun_gamma, data = polar_df_exc, family = 'ltgammab', args = list(alpha = 2,left = polar_df_exc$r_thresh), trace = 2, knots = knots)

  saveRDS(list(coeff = tg_model$coefficients,smooth_penalties = tg_model$sp ),file = paste0("data/tg_model_coeff_",s,"_",s2,".rds"))
}

tg_coeff_info = readRDS(file = paste0("data/tg_model_coeff_",s,"_",s2,".rds"))

# Specify a second dummy/redundant model (for the same reasons)  

dummy_model2 = gam(as.formula(paste0("r ~ te(phi,t,bs=c('cc','cr'),k=c(",knot_angle,",",knot_time,"))")),data = polar_df_exc,knots = knots,sp=c(0.1,0.1))

# Plot threshold quantile sets over time from model fit

pdf(file=paste0("figures/all_quantile_sets_",s,"_",s2,".pdf"),width=5,height=5)

par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

for(i in 1:length(time_points)){
  
  dummy_design = predict(dummy_model, newdata=list(phi=pred_phis,t = rep(time_points[i],length(pred_phis))),type="lpmatrix")
  
  pred_thresh = exp(c(dummy_design%*%as.matrix(qr_coeff_info$coeff)))
  
  quantile_set = as.vector(pred_thresh)*unit_circle
  
  if(i == 1){
    plot(quantile_set,pch=16,col=grad_cols[i],type="l",lwd=3,ylab = paste0(symbols2[which(symbols==s2)]), xlab = paste0(symbols2[which(symbols==s)]),main=paste0("Quantile sets, ",symbols2[which(symbols==s)],"-",symbols2[which(symbols==s2)]),cex.lab=1.2, cex.axis=1.2,cex.main=1.5,xlim=c(-5,5),ylim=c(-5,5))
  } else {
    lines(quantile_set,type="l",lwd=3,col=grad_cols[i])
  }
  
  legend("topleft", legend = c("Start", "End"),
         col = c("blue", "orange"), lwd=3,cex=1.1,bg="white")

  
}

dev.off()

# Plot estimated limit sets over time from model fit 

pdf(file=paste0("figures/all_boundary_sets_",s,"_",s2,".pdf"),width=5,height=5)

par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

for(i in 1:length(time_points)){
  
  dummy_design = predict(dummy_model2, newdata=list(phi=pred_phis,t = rep(time_points[i],length(pred_phis))),type="lpmatrix")
  
  pred_gauge = exp(c(dummy_design%*%as.matrix(tg_coeff_info$coeff)))
  
  est_limit_set = unit_circle/pred_gauge
  
  if(i == 1){
    plot(est_limit_set,pch=16,col=grad_cols[i],type="l",lwd=3,ylab = paste0(symbols2[which(symbols==s2)]), xlab = paste0(symbols2[which(symbols==s)]),main=paste0("Limit sets, ",symbols2[which(symbols==s)],"-",symbols2[which(symbols==s2)]),cex.lab=1.2, cex.axis=1.2,cex.main=1.5,xlim=c(-1.3,1.3),ylim=c(-1.3,1.3))
  } else {
    lines(est_limit_set,type="l",lwd=3,col=grad_cols[i])
  }
  
  rect(-1,-1,1,1,lwd=4,lty=2,col=NULL)
  
  legend("topleft", legend = c("Start", "End"), 
          col = c("blue", "orange"), lwd=3,cex=1.1,bg="white")
  
  
}

dev.off()

# Compute return level sets over time from fitted model 

rl_set_prob = 0.999

trunc_prob = (rl_set_prob - tau)/(1-tau)

pdf(file=paste0("figures/all_return_level_sets_",s,"_",s2,".pdf"),width=5,height=5)

par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

for(i in 1:length(time_points)){
  
  dummy_design = predict(dummy_model, newdata=list(phi=pred_phis,t = rep(time_points[i],length(pred_phis))),type="lpmatrix")
  
  pred_thresh = exp(c(dummy_design%*%as.matrix(qr_coeff_info$coeff)))
  
  dummy_design2 = predict(dummy_model2, newdata=list(phi=pred_phis,t = rep(time_points[i],length(pred_phis))),type="lpmatrix")
  
  pred_gauge = exp(c(dummy_design2%*%as.matrix(tg_coeff_info$coeff)))
  
  radial_quants = qgamma(p = pgamma(pred_thresh,shape = 2,rate = pred_gauge,lower.tail = F)*trunc_prob + pgamma(pred_thresh,shape = 2,rate = pred_gauge,lower.tail = T) , shape = 2, rate = pred_gauge )
  
  rl_set = radial_quants*unit_circle
  
  if(i == 1){
    plot(rl_set,pch=16,col=grad_cols[i],type="l",lwd=3,ylab = paste0(symbols2[which(symbols==s2)]), xlab = paste0(symbols2[which(symbols==s)]),sub=paste0("p = ",rl_set_prob),main=paste0("Return level sets, ",symbols2[which(symbols==s)],"-",symbols2[which(symbols==s2)]),cex.lab=1.2, cex.axis=1.2,cex.main=1.5,xlim=c(-10,10),ylim=c(-10,10))
  } else {
    lines(rl_set,type="l",lwd=3,col=grad_cols[i])
  }
  
  legend("topleft", legend = c("Start", "End"), 
         col = c("blue", "orange"), lwd=3,cex=1.1,bg="white")
  
  
}

dev.off()


# Computing diagnostics  --------------------------------------------------

# Evaluate truncated gamma QQ plot diagnostic 

pdf(file=paste0("figures/qq_plot_diag_",s,"_",s2,".pdf"),width=5,height=5)

dummy_design = predict(dummy_model2, newdata=polar_df_exc,type="lpmatrix")

pred_gauge = exp(c(dummy_design%*%as.matrix(tg_coeff_info$coeff)))

unif_exceedances =  1 - exp(pgamma(polar_df_exc$r,shape = 2,rate = pred_gauge,lower.tail = F,log.p=T)-pgamma(polar_df_exc$r_thresh,shape = 2,rate = pred_gauge,lower.tail = F,log.p=T))

observed_quants = qexp(unif_exceedances)

m = length(observed_quants)

theoretical_quants = qexp((1:m)/(m+1))

alpha = 0.05

lower_bound = qexp(qbeta(alpha / 2, 1:m, rev(1:m)))
upper_bound = qexp(qbeta(1 - alpha / 2, 1:m, rev(1:m)))

plot(theoretical_quants,sort(observed_quants),xlim=range(theoretical_quants,observed_quants,upper_bound,lower_bound),
     ylim=range(theoretical_quants,observed_quants,upper_bound,lower_bound),pch=16,col=1,ylab="Observed",xlab="Theoretical",
     cex.lab=1.3, cex.axis=1.2,cex.main=1.8, cex=0.5,main=paste0("QQ plot ",symbols2[which(symbols==s)],"-",symbols2[which(symbols==s2)]))
polygon(c(lower_bound,rev(upper_bound)),c(sort(observed_quants),rev(sort(observed_quants))), col = 'grey80', border = NA)
abline(a=0,b=1,lwd=3,col=2)
points(theoretical_quants,sort(observed_quants),pch=16,col=1, cex=0.5)

dev.off()

# Evaluated return level set probabilities diagnostic 

m = 200

rl_set_probs = seq(tau,0.99,length.out = m)

trunc_probs = (rl_set_probs - tau)/(1-tau)

dummy_design = predict(dummy_model, newdata=polar_df,type="lpmatrix")

pred_thresh = exp(c(dummy_design%*%as.matrix(qr_coeff_info$coeff)))

dummy_design2 = predict(dummy_model2, newdata=polar_df,type="lpmatrix")

pred_gauge = exp(c(dummy_design2%*%as.matrix(tg_coeff_info$coeff)))

emp_probs = sapply(trunc_probs,function(tp){
  radial_quants = qgamma(p = pgamma(pred_thresh,shape = 2,rate = pred_gauge,lower.tail = F)*tp + pgamma(pred_thresh,shape = 2,rate = pred_gauge,lower.tail = T) , shape = 2, rate = pred_gauge )
  return(mean(polar_df$r <= radial_quants))
})

pdf(file=paste0("figures/rl_set_diag_",s,"_",s2,".pdf"),width=5,height=5)

par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

plot(-log(1-rl_set_probs),-log(1-emp_probs),xlim=range(-log(1-rl_set_probs),-log(1-emp_probs)),
     ylim=range(-log(1-rl_set_probs),-log(1-emp_probs)),pch=16,col=1,ylab="Observed",xlab="Theoretical",
     cex.lab=1.3, cex.axis=1.2,cex.main=1.5, cex=0.5,main=paste0("RL set diagnostic ",symbols2[which(symbols==s)],"-",symbols2[which(symbols==s2)]),sub=expression(-log(1-p)))

abline(a=0,b=1,lwd=3,col=2)
points(-log(1-rl_set_probs),-log(1-emp_probs),pch=16,col=1, cex=0.5)

dev.off()

# Joint risk assessment -------------------------------------------------

model_eta_func_cwm = function(index,dummy_model2,pred_phis,tg_coeff_info,quad){ # Wrapper for computing eta from an estimated limit set 
  
  dummy_design = predict(dummy_model2, newdata=list(phi=pred_phis,t = rep(index,length(pred_phis))),type="lpmatrix")
  
  pred_gauge = exp(c(dummy_design%*%as.matrix(tg_coeff_info$coeff)))
  
  est_limit_set = unit_circle/pred_gauge
  
  if(quad == 1){
    qs = apply(est_limit_set,2,max)
  } else if(quad == 2){
    qs = c(min(est_limit_set[,1]),max(est_limit_set[,2]))
  } else if(quad == 3){
    qs = apply(est_limit_set,2,min)
  } else if(quad == 4){
    qs = c(max(est_limit_set[,1]),min(est_limit_set[,2]))
  }
  
  r_q = rect2polar(qs)$r
  
  phi_q = c(rect2polar(qs)$phi)
  
  dummy_design2 = predict(dummy_model2, newdata=list(phi=phi_q,t = index),type="lpmatrix")
  
  gauge_q = exp(c(dummy_design2%*%as.matrix(tg_coeff_info$coeff)))
  
  eta_gauge = 1/(r_q*gauge_q)
  
  return(eta_gauge)
  
}

# Computing eta across all quadrants from the model fit 

if(!file.exists(paste0("data/model_etas_quad4_",s,"_",s2,".rds"))){
  
  eta_gauge_cwm_1 = sapply(1:n_data,model_eta_func_cwm,dummy_model2=dummy_model2,pred_phis=pred_phis,tg_coeff_info=tg_coeff_info,quad=1)
  
  saveRDS(eta_gauge_cwm_1,paste0("data/model_etas_quad1_",s,"_",s2,".rds"))
  
  eta_gauge_cwm_2 = sapply(1:n_data,model_eta_func_cwm,dummy_model2=dummy_model2,pred_phis=pred_phis,tg_coeff_info=tg_coeff_info,quad=2)
  
  saveRDS(eta_gauge_cwm_2,paste0("data/model_etas_quad2_",s,"_",s2,".rds"))
  
  eta_gauge_cwm_3 = sapply(1:n_data,model_eta_func_cwm,dummy_model2=dummy_model2,pred_phis=pred_phis,tg_coeff_info=tg_coeff_info,quad=3)
  
  saveRDS(eta_gauge_cwm_3,paste0("data/model_etas_quad3_",s,"_",s2,".rds"))
  
  eta_gauge_cwm_4 = sapply(1:n_data,model_eta_func_cwm,dummy_model2=dummy_model2,pred_phis=pred_phis,tg_coeff_info=tg_coeff_info,quad=4)
  
  saveRDS(eta_gauge_cwm_4,paste0("data/model_etas_quad4_",s,"_",s2,".rds"))
  
}

eta_gauge_cwm_1 = readRDS(paste0("data/model_etas_quad1_",s,"_",s2,".rds"))

eta_gauge_cwm_2 = readRDS(paste0("data/model_etas_quad2_",s,"_",s2,".rds"))

eta_gauge_cwm_3 = readRDS(paste0("data/model_etas_quad3_",s,"_",s2,".rds"))

eta_gauge_cwm_4 = readRDS(paste0("data/model_etas_quad4_",s,"_",s2,".rds"))

# Plotting eta estimates over time 

pdf(file=paste0("figures/eta_ests_quad1_",s,"_",s2,".pdf"),width=5,height=5)

par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

plot(dates,eta_gauge_cwm_1,col="blue3",lwd=3,main=paste0(symbols2[which(symbols==s)],"-",symbols2[which(symbols==s2)]),type="l", ylab = expression(eta[t]), xlab = "Date",cex.lab=1.5, cex.axis=1.5,cex.main=1.7,cex.sub=1.2,ylim=c(0.45,1.05),sub="Quadrant 1")

abline(h = 1,lwd=2,lty=2)

dev.off()

pdf(file=paste0("figures/eta_ests_quad2_",s,"_",s2,".pdf"),width=5,height=5)

par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

plot(dates,eta_gauge_cwm_2,col="blue3",lwd=3,main=paste0(symbols2[which(symbols==s)],"-",symbols2[which(symbols==s2)]),type="l", ylab = expression(eta[t]), xlab = "Date",cex.lab=1.5, cex.axis=1.5,cex.main=1.7,cex.sub=1.2,ylim=c(0.45,1.05),sub="Quadrant 2")

abline(h = 1,lwd=2,lty=2)

dev.off()

pdf(file=paste0("figures/eta_ests_quad3_",s,"_",s2,".pdf"),width=5,height=5)

par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

plot(dates,eta_gauge_cwm_3,col="blue3",lwd=3,main=paste0(symbols2[which(symbols==s)],"-",symbols2[which(symbols==s2)]),type="l", ylab = expression(eta[t]), xlab = "Date",cex.lab=1.5, cex.axis=1.5,cex.main=1.7,cex.sub=1.2,ylim=c(0.45,1.05),sub="Quadrant 3")

abline(h = 1,lwd=2,lty=2)

dev.off()

pdf(file=paste0("figures/eta_ests_quad4_",s,"_",s2,".pdf"),width=5,height=5)

par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

plot(dates,eta_gauge_cwm_4,col="blue3",lwd=3,main=paste0(symbols2[which(symbols==s)],"-",symbols2[which(symbols==s2)]),type="l", ylab = expression(eta[t]), xlab = "Date",cex.lab=1.5, cex.axis=1.5,cex.main=1.7,cex.sub=1.2,ylim=c(0.45,1.05),sub="Quadrant 4")

abline(h = 1,lwd=2,lty=2)

dev.off()

prob = 0.99

alpha = 0.03

all_gpd_paras = readRDS(file="data/gpd_marg_paras.rds")

cycle1 = as.Date(c("2001-03-01","2001-11-30"))

cycle2 = as.Date(c("2007-12-01","2009-06-30"))

cycle3 = as.Date(c("2020-02-01","2020-04-30"))

# We do not compute the CoVaR estimates in this script since this requires a lot of CPU time + power. Please contact us if you would like access to this code

pdf(file=paste0("figures/os_covar_ests_",s,"_",s2,"_both.pdf"),width=6,height=6)

par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

gpd_paras_s = rbind(all_gpd_paras$thresh[,which(symbols==s)],all_gpd_paras$scale[,which(symbols==s)],all_gpd_paras$shape[,which(symbols==s)])

gpd_paras_s2 = rbind(all_gpd_paras$thresh[,which(symbols==s2)],all_gpd_paras$scale[,which(symbols==s2)],all_gpd_paras$shape[,which(symbols==s2)])

# We load in the CoVaR estimates on Laplace margins, then transform back to the original scale of the log-returns by apply the probability integral transform

covar_over_t_1 = readRDS(file = paste0("data/model_covar_",s,"_",s2,"_quad_",1,".rds")) 

os_var_s_1 = inverse_semi_empirical_cdf(u=prob,alpha=alpha,gpd_paras = gpd_paras_s,marg_vec = filtered_returns_datasets[[s]])

covar_s2_1 = inverse_semi_empirical_cdf(u=Laplace_cdf(covar_over_t_1),alpha=alpha,gpd_paras = gpd_paras_s2,marg_vec = filtered_returns_datasets[[s2]])

os_var_s_1 = os_var_s_1*garch_paras[[s]]$sigma + garch_paras[[s]]$mu

covar_s2_1 = covar_s2_1*garch_paras[[s2]]$sigma + garch_paras[[s2]]$mu

covar_over_t_3 = readRDS(file = paste0("data/model_covar_",s,"_",s2,"_quad_",3,".rds"))

os_var_s_3 = inverse_semi_empirical_cdf(u=1-prob,alpha=alpha,gpd_paras = gpd_paras_s,marg_vec = filtered_returns_datasets[[s]])

covar_s2_3 = inverse_semi_empirical_cdf(u=Laplace_cdf(covar_over_t_3),alpha=alpha,gpd_paras = gpd_paras_s2,marg_vec = filtered_returns_datasets[[s2]])

os_var_s_3 = os_var_s_3*garch_paras[[s]]$sigma + garch_paras[[s]]$mu

covar_s2_3 = covar_s2_3*garch_paras[[s2]]$sigma + garch_paras[[s2]]$mu

plot(dates,covar_s2_1,main=paste0(symbols2[which(symbols==s)]," vs ",symbols2[which(symbols==s2)]),type="l", lwd=1,col=3, ylab = "", xlab = "Date",cex.lab=1.5, cex.axis=1.5,cex.main=1.7,ylim=range(covar_s2_1,os_var_s_1,covar_s2_3,os_var_s_3))

lines(dates,os_var_s_1,lwd=1,col=2)

lines(dates,covar_s2_3,lwd=1,lty=2,col=3)

lines(dates,os_var_s_3,lwd=1,col=2,lty=2)

usr = par("usr")

rect(xleft = cycle1[1], 
     ybottom = par("usr")[3],   
     xright = cycle1[2], 
     ytop = par("usr")[4],      
     col = rgb(0, 0, 1, 0.2),   
     border = NA)

text(x = cycle1[1]-150, 
     y = usr[3]*0.9, 
     labels = "P", 
     pos = 3, cex = 0.9,font=2)  

text(x = cycle1[2]+150, 
     y = usr[3]*0.9, 
     labels = "T", 
     pos = 3, cex = 0.9,font = 2)

rect(xleft = cycle2[1], 
     ybottom = par("usr")[3],   
     xright = cycle2[2], 
     ytop = par("usr")[4],      
     col = rgb(0, 0, 1, 0.2),   
     border = NA)

text(x = cycle2[1]-150, 
     y = usr[3]*0.9, 
     labels = "P", 
     pos = 3, cex = 0.9,font=2)  

text(x = cycle2[2]+150, 
     y = usr[3]*0.9, 
     labels = "T", 
     pos = 3, cex = 0.9,font = 2)

rect(xleft = cycle3[1], 
     ybottom = par("usr")[3],   
     xright = cycle3[2], 
     ytop = par("usr")[4],      
     col = rgb(0, 0, 1, 0.2),   
     border = NA)

text(x = cycle3[1]-150, 
     y = usr[3]*0.9, 
     labels = "P", 
     pos = 3, cex = 0.9,font=2)  

text(x = cycle3[2]+150, 
     y = usr[3]*0.9, 
     labels = "T", 
     pos = 3, cex = 0.9,font = 2)

legend("topright", legend = c("Upside VaR","Upside CoVaR","Downside VaR","Downside CoVaR"), 
       col = c(2,3,2,3),lty=c(1,1,2,2), lwd=2,cex=1,bg="white")

dev.off()

data_os = cbind(filtered_returns_datasets[[s]],filtered_returns_datasets[[s2]]) # Matrix of filtered returns data

rl_set_prob = 0.999 # Probability for computing return level sets 

trunc_prob = (rl_set_prob - tau)/(1-tau) # Corresponding truncated probability for truncated Gamma model 

signif_days = as.Date(c("2015-08-24","2020-03-16","2022-02-24")) # Dates of significance 

event_names = c("Flash crash","COVID pandemic","Russian invasion") # Names of significant events 

leg_names = paste0(signif_days,", ",event_names)

time_points = which(dates %in% signif_days)

grad_cols = c("#0072B2", "#E69F00", "#D55E00")

pdf(file=paste0("figures/os_return_level_sets_signif_days_",s,"_",s2,".pdf"),width=5,height=5)

par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

for(i in 1:length(time_points)){
  
  dummy_design = predict(dummy_model, newdata=list(phi=pred_phis,t = rep(time_points[i],length(pred_phis))),type="lpmatrix")
  
  pred_thresh = exp(c(dummy_design%*%as.matrix(qr_coeff_info$coeff)))
  
  dummy_design2 = predict(dummy_model2, newdata=list(phi=pred_phis,t = rep(time_points[i],length(pred_phis))),type="lpmatrix")
  
  pred_gauge = exp(c(dummy_design2%*%as.matrix(tg_coeff_info$coeff)))
  
  radial_quants = qgamma(p = pgamma(pred_thresh,shape = 2,rate = pred_gauge,lower.tail = F)*trunc_prob + pgamma(pred_thresh,shape = 2,rate = pred_gauge,lower.tail = T) , shape = 2, rate = pred_gauge )
  
  rl_set = radial_quants*unit_circle
  
  rl_set_os = rl_set
  
  rl_set_os[,1] = inverse_semi_empirical_cdf(u=Laplace_cdf(rl_set[,1]),alpha=alpha,gpd_paras = rbind(all_gpd_paras$thresh[,which(symbols==s)],all_gpd_paras$scale[,which(symbols==s)],all_gpd_paras$shape[,which(symbols==s)]),marg_vec = filtered_returns_datasets[[s]])
  
  rl_set_os[,2] = inverse_semi_empirical_cdf(u=Laplace_cdf(rl_set[,2]),alpha=alpha,gpd_paras = rbind(all_gpd_paras$thresh[,which(symbols==s2)],all_gpd_paras$scale[,which(symbols==s2)],all_gpd_paras$shape[,which(symbols==s2)]),marg_vec = filtered_returns_datasets[[s2]])
  
  rl_set_os[,1] = garch_paras[[s]]$sigma[time_points[i]]*rl_set_os[,1] + garch_paras[[s]]$mu[time_points[i]]
  
  rl_set_os[,2] = garch_paras[[s2]]$sigma[time_points[i]]*rl_set_os[,2] + garch_paras[[s2]]$mu[time_points[i]]
  
  if(i == 1){
    plot(rl_set_os,pch=16,col=grad_cols[i],type="l",lwd=4,ylab = paste0(symbols2[which(symbols==s2)]), xlab = paste0(symbols2[which(symbols==s)]),sub=paste0("p = ",rl_set_prob),main=paste0("Return level sets, ",symbols2[which(symbols==s)],"-",symbols2[which(symbols==s2)]),cex.lab=1.2, cex.axis=1.2,cex.main=1.5,xlim=c(-0.6,0.6),ylim=c(-0.6,0.6))
  } else {
    lines(rl_set_os,type="l",lwd=4,col=grad_cols[i])
  }
  
  legend("topleft", legend = event_names, 
    col = grad_cols, lwd=4,cex=1.1,bg="white")
  
  
  
}

dev.off()

# Modelling angular distribution over time -------------------------------------------

num_win_ang = 10 # Number of local non-overlapping windows for assessing angular fits

bw_t = 150 # Bandwidth for Gaussian kernel 

kappa = 5 # Bandwidth/concentration parameter for von-Mises kernel

pdf(file=paste0("figures/circ_dens_nonpara_regrssion_",s,"_",s2,".pdf"),width=15,height=6)

par(mfrow=c(2,5),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

# Compute a histogram for each time window, compare to kernel density estimate at the mid-point of the window

for(i in 1:num_win_ang){
  
  inds = (round((i-1)*n_data/num_win_ang) + 1):(round(i*n_data/num_win_ang))
  
  t0 = round(median(inds))
  
  phi_window = polar_df_exc$phi[ which(polar_df_exc$t %in% inds) ]
  
  weights = dnorm(polar_df_exc$t, mean = t0, sd = bw_t)
  
  kde = sapply(pred_phis, function(phi0) {
    sum(weights * dvonmises(circular(polar_df_exc$phi), mu = circular(phi0), kappa = kappa)) / sum(weights)
  })
  
  hist_data = hist(phi_window, plot=F) 
  
  hist(phi_window, freq = FALSE,breaks = 10,xlab=expression(Phi),ylab=expression(f[Phi](phi)),col="lightblue", main=paste0("Angular density, ",format(dates[t0], "%d/%m/%Y")),sub=paste0("Window: ",format(dates[min(inds)], "%d/%m/%Y"),"-",format(dates[max(inds)], "%d/%m/%Y")),cex.lab=1.5, cex.axis=1.5,cex.main=1.7,ylim=c(0,0.45))#,ylim=range(kde,0,hist_data$density)
  
  lines(pred_phis,kde,lwd=3,col=2)
  
}

dev.off()

# Simulating data in the joint tail -----------------------------------------------------

geometric_simulation_function = function(time_index,nsim=1){ # Wrapper for simulating from the fitted geometric model 
  
  # time_index denotes the time point at which you wish to simulate new observations 
  
  # nsim denotes the number of points you wish to simulate at this time point
  
  weights <- dnorm(polar_df_exc$t, mean = time_index, sd = bw_t)
  
  kde <- sapply(pred_phis, function(phi0) {
    sum(weights * dvonmises(circular(polar_df_exc$phi), mu = circular(phi0), kappa = kappa)) / sum(weights)
  })
  
  f_phi = approxfun(x = pred_phis,y = kde,method = "constant") # Approximate the density function 
  
  max_cdf_val = integrate(f_phi,0,2*pi,subdivisions = 2000,rel.tol = 1e-3,abs.tol = 1e-3)$value # Compute integrated area under density curve 
  
  kd_integral = function(u,x,f_phi,max_cdf_val){ # Integral of estaimted density function, normalised to exist in [0,1]
    return((integrate(f_phi,0,x,subdivisions = 2000,rel.tol = 1e-3,abs.tol = 1e-3)$value)/max_cdf_val - u) #to ensure we integrate to 1
  }
  
  kd_root = function(u,f_phi,max_cdf_val){ # Rootfinder function for inverting the estimated distibution function 
    return(uniroot(f=kd_integral,interval = c(0,2*pi),u=u,f_phi = f_phi,max_cdf_val=max_cdf_val)$root)
  }
  
  gen_phi = sapply(runif(nsim),kd_root,f_phi=f_phi,max_cdf_val=max_cdf_val) # Generate new sample of angles 
  
  # The generated angles are inputted into the truncated Gamma model and used to generate a radial value 
  
  dummy_design = predict(dummy_model, newdata=list(phi=gen_phi,t = rep(time_index,nsim)),type="lpmatrix")
  
  pred_thresh = exp(c(dummy_design%*%as.matrix(qr_coeff_info$coeff)))
  
  dummy_design2 = predict(dummy_model2, newdata=list(phi=gen_phi,t = rep(time_index,nsim)),type="lpmatrix")
  
  pred_gauge = exp(c(dummy_design2%*%as.matrix(tg_coeff_info$coeff)))
  
  r_sample = qgamma(p = pgamma(pred_thresh,shape = 2,rate = pred_gauge,lower.tail = F)*runif(nsim) + pgamma(pred_thresh,shape = 2, rate = pred_gauge,lower.tail = T) , shape = 2, rate = pred_gauge )
  
  # Transform from polar to Cartesian scale  
  
  data_lap_sample = r_sample*cbind(cos(gen_phi),sin(gen_phi))
  
  return(cbind(data_lap_sample,rep(time_index,nsim)))
  
}

# For each time point, we simulate one data point in the joint tail. We compare this to the observed data

if(!file.exists("data/simulated_data_over_time.rds")){
  simulated_samples = t(sapply(time,geometric_simulation_function,nsim = 1,simplify = T))
  
  saveRDS(simulated_samples,file="data/simulated_data_over_time.rds")
  
}

simulated_samples = readRDS(file="data/simulated_data_over_time.rds")

pdf(file="figures/model_simulated_data.pdf",width=5,height=5)

num_win = 10

for(i in 1:num_win){
  
  par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)
  
  inds = (round((i-1)*n_data/num_win) + 1):(round(i*n_data/num_win))
  
  window_sim_data = simulated_samples[ which(simulated_samples[,3] %in% inds) ,1:2]
  
  plot(laplace_datasets[[s]][inds],laplace_datasets[[s2]][inds],type="p", col="grey",pch=16,cex=1.5, ylab = paste0(symbols2[which(symbols==s2)]), xlab = paste0(symbols2[which(symbols==s)]),cex.lab=1.5, cex.axis=1.5,cex.main=1.7,main=paste0("Window: ",format(dates[((i-1)*nrow(data_lap)/num_win + 1)], "%d/%m/%Y"),"-",format(dates[(i*nrow(data_lap)/num_win)], "%d/%m/%Y")),xlim=c(-8,8),ylim=c(-8,8))
  
  points(window_sim_data,type="p", col="green4",pch=16)
  
  legend("topleft",legend = c("Observed", "Simulated"),col = c("grey","green4"),pch=c(16,16),cex=1.1,bg = "white")
  
}

dev.off()

if(!file.exists("data/simulated_data_COVID_onset.rds")){
  
  COVID_t = which(dates == "2020-03-16") # Time index of COVID onset 
  
  nsim = 10000 # Number of joint tail simulations for COVID onset 
  
  simulated_data_COVID = geometric_simulation_function(COVID_t,nsim = nsim)
  
  saveRDS(simulated_data_COVID,file="data/simulated_data_COVID_onset.rds")
  
}

simulated_data_COVID = readRDS(file="data/simulated_data_COVID_onset.rds")

pdf(file="figures/COVID_simulated_data.pdf",width=5,height=5)

par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

plot(simulated_data_COVID[,1:2],type="p", col="green4",pch=16,cex=.8, ylab = paste0(symbols2[which(symbols==s2)]), xlab = paste0(symbols2[which(symbols==s)]),cex.lab=1.5, cex.axis=1.5,cex.main=1.7,main=paste0("Onset of COVID, ",format(dates[COVID_t], "%d/%m/%Y")),xlim=c(-10,10),ylim=c(-10,10))

dev.off()