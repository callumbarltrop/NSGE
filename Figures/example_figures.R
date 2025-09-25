# Packages and functions ----------------------------------------------------------------

source("preamble.R")

library(rgl)

sim_data_normal_fn = function(rho){
  d = 2
  normc = normalCopula(param = rho, dim = d)
  return(apply(rCopula(1, copula = normc),2,Laplace_inverse))
}

pred_phis = seq(0,2*pi,length.out=251)

# Theoretical gauges ------------------------------------------------------

set.seed(1)

n = 10000

rho = 0.5

normc = normalCopula(param = rho, dim = 2)

sim_data_normal = apply(rCopula(n, copula = normc),2,Laplace_inverse)

gauge_unit_circle = apply(cbind(cos(pred_phis),sin(pred_phis)),1, gauge_normal,rho=rho)

limit_set_normal = cbind(cos(pred_phis)/gauge_unit_circle,sin(pred_phis)/gauge_unit_circle)

rho = 0.5
nu = 1.5

gauge_unit_circle = apply(cbind(cos(pred_phis),sin(pred_phis)),1, gauge_mv_t,nu=nu)

limit_set_t = cbind(cos(pred_phis)/gauge_unit_circle,sin(pred_phis)/gauge_unit_circle) #by homogeniety

tc = tCopula(param = rho,df=nu, dim = 2)

sim_data_t = apply(rCopula(n, copula = tc),2,Laplace_inverse)

dep = 0.4

gauge_unit_circle = apply(cbind(cos(pred_phis),sin(pred_phis)),1, gauge_inv_logistic,dep=dep)

limit_set_invlog = cbind(cos(pred_phis)/gauge_unit_circle,sin(pred_phis)/gauge_unit_circle) #by homogeniety

sim_data_invlog = apply(pexp(1/rbvevd(n=n,dep=dep,model="log",mar1 = c(1,1,1))),2,Laplace_inverse)

pdf(file="theoretical_gauges.pdf",width=15,height=5)

#Setting plotting parameters
par(mfrow=c(1,3),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

plot(sim_data_normal/(log(n/2)),xlim=c(-1,1),ylim=c(-1,1),xlab=expression(X[1]/r[n]),ylab=expression(X[2]/r[n]),pch=16,col="grey",main="Gaussian copula",cex.lab=1.5, cex.axis=1.5,cex.main=1.7)
rect(-1,-1,1,1,lwd=4,lty=2,col=NULL)
polygon(limit_set_normal[,1],limit_set_normal[,2],col=rgb(0, 0, 255, max = 255, alpha = 15),border = 0,lwd=3)
lines(limit_set_normal,lwd=4,col=4)

plot(sim_data_t/(log(n/2)),xlim=c(-1,1),ylim=c(-1,1),xlab=expression(X[1]/r[n]),ylab=expression(X[2]/r[n]),pch=16,col="grey",main="Student-t copula",cex.lab=1.5, cex.axis=1.5,cex.main=1.7)
rect(-1,-1,1,1,lwd=4,lty=2,col=NULL)
polygon(limit_set_t[,1],limit_set_t[,2],col=rgb(0, 0, 255, max = 255, alpha = 15),border = 0,lwd=3)
lines(limit_set_t,lwd=4,col=4)

plot(sim_data_invlog/(log(n/2)),xlim=c(-1,1),ylim=c(-1,1),xlab=expression(X[1]/r[n]),ylab=expression(X[2]/r[n]),pch=16,col="grey",main="Inverted logistic copula",cex.lab=1.5, cex.axis=1.5,cex.main=1.7)
rect(-1,-1,1,1,lwd=4,lty=2,col=NULL)
polygon(limit_set_invlog[,1],limit_set_invlog[,2],col=rgb(0, 0, 255, max = 255, alpha = 15),border = 0,lwd=3)
lines(limit_set_invlog,lwd=4,col=4)

dev.off()


# Time changing limit sets ------------------------------------------------

pred_phis = seq(0,2*pi,length.out=601)

unit_circle <- cbind(cos(pred_phis),sin(pred_phis))

rhos = seq(0.2,0.8,length.out=n)

time_indices = seq(1,n,length.out=10)

ns_gauges_normal = sapply(rhos[time_indices],function(p){
  return(unit_circle/apply(unit_circle,1, gauge_normal,rho=p))
},simplify = F)
  
alphas = seq(0.3,0.7,length.out=n)

ns_gauges_invlog = sapply(alphas[time_indices],function(p){
  return(unit_circle/apply(unit_circle,1, gauge_inv_logistic,dep=p))
},simplify = F)
  
rho = 0.5
nus = seq(0.5,2,length.out=n)

ns_gauges_t = sapply(nus[time_indices],function(p){
  return(unit_circle/apply(unit_circle,1, gauge_mv_t,nu=p))
},simplify = F)
 
t1 = 0.45*(0:(n-1))/(n-1)
t2 = sin((5*pi/2)*(0:(n-1))/(n-1))/2

rhos = t1 + t2
  
ns_gauges_normal2 = sapply(rhos[time_indices],function(p){
  return(unit_circle/apply(unit_circle,1, gauge_normal,rho=p))
},simplify = F)
  
gammas = seq(0.7,2.5,len=n)
rho = 0.5
  
ns_gauges_hw = sapply(gammas[time_indices],function(p){
  return(unit_circle/apply(unit_circle,1, gauge_hw,par=c(p,rho)))
},simplify = F)

# Create a color function using colorRampPalette
colfunc <- colorRampPalette(c("blue", "orange"))

# Generate a vector of 10 colors from the palette
cols <- colfunc(length(time_indices))

pdf(file="ns_limit_sets.pdf",width=15,height=10)

#Setting plotting parameters
par(mfrow=c(2,3),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

plot(ns_gauges_normal[[1]],xlim=c(-1,1),ylim=c(-1,1),xlab=expression(X[1]),ylab=expression(X[2]),type="l",lwd=2,col=cols[1],main="Copula 1",cex.lab=1.5, cex.axis=1.5,cex.main=1.7)
rect(-1,-1,1,1,lwd=4,lty=2,col=NULL)
for(i in 2:length(time_indices)){
  lines(ns_gauges_normal[[i]],lwd=2,col=cols[i])
}

plot(ns_gauges_normal2[[1]],xlim=c(-1,1),ylim=c(-1,1),xlab=expression(X[1]),ylab=expression(X[2]),type="l",lwd=2,col=cols[1],main="Copula 2",cex.lab=1.5, cex.axis=1.5,cex.main=1.7)
rect(-1,-1,1,1,lwd=4,lty=2,col=NULL)
for(i in 2:length(time_indices)){
  lines(ns_gauges_normal2[[i]],lwd=2,col=cols[i])
}

plot(ns_gauges_invlog[[1]],xlim=c(-1,1),ylim=c(-1,1),xlab=expression(X[1]),ylab=expression(X[2]),type="l",lwd=2,col=cols[1],main="Copula 3",cex.lab=1.5, cex.axis=1.5,cex.main=1.7)
rect(-1,-1,1,1,lwd=4,lty=2,col=NULL)
for(i in 2:length(time_indices)){
  lines(ns_gauges_invlog[[i]],lwd=2,col=cols[i])
}

plot(ns_gauges_t[[1]],xlim=c(-1,1),ylim=c(-1,1),xlab=expression(X[1]),ylab=expression(X[2]),type="l",lwd=2,col=cols[1],main="Copula 4",cex.lab=1.5, cex.axis=1.5,cex.main=1.7)
rect(-1,-1,1,1,lwd=4,lty=2,col=NULL)
for(i in 2:length(time_indices)){
  lines(ns_gauges_t[[i]],lwd=2,col=cols[i])
}

plot(ns_gauges_hw[[1]],xlim=c(-1,1),ylim=c(-1,1),xlab=expression(X[1]),ylab=expression(X[2]),type="l",lwd=2,col=cols[1],main="Copula 5",cex.lab=1.5, cex.axis=1.5,cex.main=1.7)
rect(-1,-1,1,1,lwd=4,lty=2,col=NULL)
for(i in 2:length(time_indices)){
  lines(ns_gauges_hw[[i]],lwd=2,col=cols[i])
}

plot(1, 1, type = "n", xlab = "", ylab = "", axes = FALSE)

# Add a legend
legend("center", legend = c("Start of time frame", "End of time frame"), 
       col = c("blue", "orange"), lwd=4,cex=2.3)

dev.off()

pdf(file="normal_ns_limit_sets.pdf",width=5,height=5)

#Setting plotting parameters
par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

plot(ns_gauges_normal[[1]],xlim=c(-1,1),ylim=c(-1,1),xlab=expression(X[1]),ylab=expression(X[2]),type="l",lwd=2,col=cols[1],main="",cex.lab=1.2, cex.axis=1.5,cex.main=1.7)
rect(-1,-1,1,1,lwd=4,lty=2,col=NULL)
for(i in 2:length(time_indices)){
  lines(ns_gauges_normal[[i]],lwd=2,col=cols[i])
}

# Add a legend
legend("topleft", legend = c("Start of time frame", "End of time frame"), 
       col = c("blue", "orange"), lwd=4,cex=1.2,bg="white")

dev.off()


# Changing data -----------------------------------------------------------

rhos = seq(-.8,0.8,length.out=n)

set.seed(2311)
sim_data = t(sapply(rhos,sim_data_normal_fn))

num_windows = 8

# Create a color function using colorRampPalette
colfunc <- colorRampPalette(c("blue", "orange"))

# Generate a vector of 10 colors from the palette
cols <- colfunc(num_windows)

for(i in 1:num_windows){
  pdf(file=paste0("ns_data",i,".pdf"),width=5,height=5)
  
  plot(sim_data[1:(n/num_windows),],xlim=c(-10,10),ylim=c(-10,10),xlab=expression(X[1]),ylab=expression(X[2]),type="p",pch=16,col=cols[1],main=paste0("Obsevation window ",i),cex.lab=1.2, cex.axis=1.5,cex.main=1.7)
  
  if(i > 1){
    for(j in 2:i){
      points(sim_data[((j-1)*(n/num_windows)+1):(j*(n/num_windows)),],pch=16,col=cols[j])
    }
  }
  
  dev.off()
}

# Rho function --------------------------------------------------------

n = 10000

t1 = 0.45*(0:(n-1))/(n-1)
t2 = sin((5*pi/2)*(0:(n-1))/(n-1))/2

pdf(file="two_cov.pdf",width=5,height=5)

par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

rhos = t1 + t2

plot(1:n,rhos,xlab="t",xaxt="n",ylab=expression(rho(t)),cex.lab=1.2, cex.axis=1.2,cex.main=1.5,type="l",lwd=4,col="red")

axis(side = 1, at = c(0, n/2, n), labels = c(1, expression(T/2), expression(T)))

dev.off()


# Different norms visualisation -------------------------------------------

pred_phis = seq(0,2*pi,length.out=601)

set.seed(1)

n = 5000

rho = 0.5

normc = normalCopula(param = rho, dim = 2)

sim_data_normal = apply(rCopula(n, copula = normc),2,Laplace_inverse)

polar_df = rect2polar(t(sim_data_normal))

polar_df$phi = c(polar_df$phi)

polar_df = as.data.frame(polar_df)

polar_df = polar_df[order(polar_df$phi),]

#non-exceedance prob for quantile regression
tau = 0.8

#Quantile regression formula for qgam 
fmla_qgam = "r ~ s(phi,bs='cc',k=15)"

#Specifying the knots 
phi_knots = seq(0,2*pi,length.out=15)

knots = list(phi = phi_knots)

qr_model = qgam(as.formula(fmla_qgam), data=polar_df,  qu=tau,argGam = list(knots = knots) )

#Computing estimated threshold function
thresh_function = c(predict(qr_model, newdata=data.frame(phi=pred_phis)))

unit_ball_l1 = cbind(cos(pred_phis),sin(pred_phis))
unit_ball_l1 = unit_ball_l1/apply(unit_ball_l1,1,l1_norm)

shape_l1 = thresh_function*unit_ball_l1

unit_ball_l2 = cbind(cos(pred_phis),sin(pred_phis))

shape_l2 = thresh_function*unit_ball_l2

unit_ball_linf = cbind(cos(pred_phis),sin(pred_phis))
unit_ball_linf = unit_ball_linf/apply(unit_ball_linf,1,linf_norm)

shape_linf = thresh_function*unit_ball_linf



pdf(file="different_norms.pdf",width=15,height=5)

#Setting plotting parameters
par(mfrow=c(1,3),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

plot(pred_phis,thresh_function,col="black",type="l",lwd=4,xlab=expression(Phi),ylab=expression(R),main="Polar coordinates",cex.lab=1.5, cex.axis=1.5,cex.main=1.7)

plot(shape_l1,lwd=4,type="l",col="blue",xlab=expression(X[1]),ylab=expression(X[2]),main="Cartesian coordinates",cex.lab=1.5, cex.axis=1.5,cex.main=1.7,ylim=range(shape_l1,shape_l2,shape_linf),xlim=range(shape_l1,shape_l2,shape_linf))

lines(shape_l2,lwd=4,col="orange")

lines(shape_linf,lwd=4,col="purple")

plot(unit_ball_l1,lwd=4,type="l",col="blue",xlab=expression(X[1]),ylab=expression(X[2]),main="Unit balls",cex.lab=1.5, cex.axis=1.5,cex.main=1.7,ylim=range(unit_ball_l1,unit_ball_l2,unit_ball_linf),xlim=range(unit_ball_l1,unit_ball_l2,unit_ball_linf))

lines(unit_ball_l2,lwd=4,col="orange")

lines(unit_ball_linf,lwd=4,col="purple")

# Add a legend
legend("topleft", legend = c(expression(L^1), expression(L^2),expression(L^{infinity})), 
       col = c("blue", "orange","purple"), lwd=4,cex=1.7,bg="white")

dev.off()


# Different knot locations angle ------------------------------------------

knot_ang1 = seq(0,2*pi,len=9)

knot_ang2 = seq(0,2*pi,len=17)

pdf(file="knot_locations_angle.pdf",width=10,height=5)

#Setting plotting parameters
par(mfrow=c(1,2),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

plot(1, 1, type = "n",xlab=expression(Phi),ylab=expression(R),main="Polar coordinates",cex.lab=1.5, cex.axis=1.5,cex.main=1.7,xlim=c(0,2*pi),ylim=c(0,10))

sapply(knot_ang1,function(phi){abline(abline(v = phi,lwd=4,lty=2,col="blue"))})

sapply(knot_ang2,function(phi){abline(abline(v = phi,lwd=3,lty=3,col="orange"))})

# Add a legend
legend("topleft", legend = c(expression(kappa[phi] == 9), expression(kappa[phi] == 17)), 
       col = c("blue", "orange"), lwd=4:3,lty=2:3,cex=1.4,bg="white")

plot(1, 1, type = "n",xlab=expression(X[1]),ylab=expression(X[2]),main="Cartesian coordinates",cex.lab=1.5, cex.axis=1.5,cex.main=1.7,ylim=c(-10,10),xlim=c(-10,10))

sapply(knot_ang1,function(phi){abline(a = 0,b=tan(phi),lwd=4,lty=2,col="blue")})

sapply(knot_ang2,function(phi){abline(a = 0,b=tan(phi),lwd=3,lty=3,col="orange")})

dev.off()


# Tensor product visualiation ---------------------------------------------

kn = 10

n = 100

phi = seq(0,2*pi,len=n)

t = seq(1,n,len = n)

# Define the smooth specification for cubic regression splines
smooth_spec_phi <- s(phi, bs = "cc", k = kn)  # k = number of basis functions

# Construct the basis
smooth_basis_phi <- smooth.construct(smooth_spec_phi, data = data.frame(phi = phi), knots = NULL)

# Define the smooth specification for cubic regression splines
smooth_spec_time <- s(t, bs = "cr", k = kn)  # k = number of basis functions

# Construct the basis
smooth_basis_time <- smooth.construct(smooth_spec_time, data = data.frame(t = t), knots = NULL)

# Extract the design matrix
X_cc <- PredictMat(smooth_basis_phi, data = data.frame(phi = phi))

X_cr <- PredictMat(smooth_basis_time, data = data.frame(t = t))

k1 <- k2 <- 5

X_cc <- PredictMat(smooth_basis_phi, data = data.frame(phi = phi))[,k1]

X_cr <- PredictMat(smooth_basis_time, data = data.frame(t = t))[,k2]

basis_1 = outer(X_cc,X_cr)

k1 <- k2 <- 9

X_cc <- PredictMat(smooth_basis_phi, data = data.frame(phi = phi))[,k1]

X_cr <- PredictMat(smooth_basis_time, data = data.frame(t = t))[,k2]

basis_2 = outer(X_cc,X_cr)

open3d()

mfrow3d(1, 2)

persp3d(phi,t,basis_1,alpha=.6, col="blue",xlab=expression(Phi),ylab="t",zlab="")

view3d(zoom=.7)

persp3d(phi,t,basis_2,alpha=.6, col="red",xlab=expression(Phi),ylab="t",zlab="")

view3d(zoom=.7)

snapshot3d("tensor_product_fig.png")

clear3d()



# Live demonstration ------------------------------------------------------

set.seed(2311)

n = 5000

pred_phis = seq(0,2*pi,length.out=201)

d = 2

knot_angle = 9

knot_time = 5

tau = 0.7

time = 1:n

rhos = seq(0,0.8,length.out=n)

sim_data = t(sapply(rhos,sim_data_normal_fn))

polar = rect2polar(t(sim_data))

#Forming data frame of polar coordinates
polar_df = data.frame(r=polar$r,phi = polar$phi[1,],logr=log(polar$r),t = time)

#Quantile regression formula for qgam 
fmla_qgam = paste0("logr ~ te(phi,t,bs=c('cc','cr'),k=c(",knot_angle,",",knot_time,"))")

#Specifying the knots 
phi_knots = seq(0,2*pi,length.out=knot_angle)

time_knots = seq(min(time),max(time),length.out=knot_time)

knots = list(phi = phi_knots,t = time_knots)

#Fitting qgam
if(!file.exists("example_tg.rds")){
  qr_model = qgam(as.formula(fmla_qgam), data=polar_df,  qu=tau,argGam = list(knots = knots) )
  
  saveRDS(qr_model,file="example_qr.rds")
}

qr_model = readRDS("example_qr.rds")

#Computing estimated threshold function
thresh_function = exp(predict(qr_model, newdata=polar_df))


#Computing exceedances of threshold function
thresh_exceedances = polar_df$r - thresh_function

print("Observed threshold exceedance probability = ")
print(mean(thresh_exceedances>0))

col_vec = c()
col_vec[thresh_exceedances>0] = "green"
col_vec[thresh_exceedances<=0] = "grey"

pred_phis = seq(0,2*pi,length.out=201)
pred_time = seq(1,n,length.out=201)
pred_grid = expand.grid(pred_phis,pred_time)
pred_df = list(phi = pred_grid$Var1,t = pred_grid$Var2)

pred_thresh = predict(qr_model, newdata=pred_df)
thresh_surface = list(phi = pred_grid$Var1,time_cov = pred_grid$Var2,r=pred_thresh)

open3d()
plot3d(x=polar_df$phi,y=polar_df$t,z=polar_df$r,col="grey",pch=16,cex=2,xlab=expression(phi),ylab="t",zlab="r")

clear3d()

plot3d(x=polar_df$phi,y=polar_df$t,z=polar_df$r,col=col_vec,pch=16,cex=2,xlab=expression(phi),ylab="t",zlab="r")
surface3d(x=matrix(thresh_surface$phi, nrow=sqrt(dim(pred_grid)[1]), ncol=sqrt(dim(pred_grid)[1])),
          y=matrix(thresh_surface$time_cov, nrow=sqrt(dim(pred_grid)[1]), ncol=sqrt(dim(pred_grid)[1])),
          z=matrix(exp(thresh_surface$r), nrow=sqrt(dim(pred_grid)[1]), ncol=sqrt(dim(pred_grid)[1])), alpha=.3, col="blue")

# close3d()
clear3d()

#checking which exceedances are positive 
positive_exceedances_indicator = which(polar_df$r - thresh_function>0)

#Obtaining exceedances
polar_df_exc = data.frame(r=polar_df$r[positive_exceedances_indicator],phi = (polar$phi[1,])[positive_exceedances_indicator],r_thresh=thresh_function[positive_exceedances_indicator],t=time[positive_exceedances_indicator])

#Defining formula for truncated gamma
fmla_trun_gamma = list(gauge = as.formula(paste0("r ~ te(phi,t,bs=c('cc','cr'),k=c(",knot_angle,",",knot_time,"))")))

#Fitting evgam truncated gamma model 
if(!file.exists("example_tg.rds")){
  tg_model <- evgam(fmla_trun_gamma, data = polar_df_exc, family = 'ltgammab', trace = 2, args = list(left = polar_df_exc$r_thresh, alpha = 2), knots = knots)
  
  saveRDS(tg_model,file="example_tg.rds")
}

tg_model = readRDS("example_tg.rds")

pred_gauge = as.vector(exp(predict(tg_model, newdata=pred_df)))

gauge_surface = list(phi = pred_grid$Var1,time_cov = pred_grid$Var2,r=as.vector(pred_gauge))

gauge_surface = as.data.frame(gauge_surface)


plot3d(x=polar_df$phi,y=polar_df$t,z=polar_df$r/log(n/2),col="grey",pch=16,cex=2,xlab=expression(phi),ylab="t",zlab="r")
surface3d(x=matrix(gauge_surface$phi, nrow=sqrt(dim(pred_grid)[1]), ncol=sqrt(dim(pred_grid)[1])),
          y=matrix(gauge_surface$time_cov, nrow=sqrt(dim(pred_grid)[1]), ncol=sqrt(dim(pred_grid)[1])),
          z=matrix(1/gauge_surface$gauge, nrow=sqrt(dim(pred_grid)[1]), ncol=sqrt(dim(pred_grid)[1])), alpha=.3, col="blue")

close3d()

time_points = round(seq(1,n,length.out = 10))

unit_circle = cbind(cos(pred_phis),sin(pred_phis))

colfunc <- colorRampPalette(c("blue", "orange"))

grad_cols <- colfunc(length(time_points))

pdf(file="example_ns_limit_sets.pdf",width=5,height=5)

par(mfrow=c(1,1),mgp=c(2.5,1,0),mar=c(5,4,4,2)+0.1)

for(i in 1:length(time_points)){
  
  pred_gauge = exp(predict(tg_model, newdata=list(phi=pred_phis,t = rep(time_points[i],length(pred_phis))))$gauge)
  
  est_limit_set = unit_circle/pred_gauge
  
  if(i == 1){
    plot(est_limit_set,pch=16,col=grad_cols[i],type="l",lwd=3,ylab = expression(X[2]), xlab = expression(X[1]),main="Limit sets",cex.lab=1.2, cex.axis=1.2,cex.main=1.5,xlim=c(-1,1),ylim=c(-1,1))
  } else {
    lines(est_limit_set,type="l",lwd=3,col=grad_cols[i])
  }
  
  rect(-1,-1,1,1,lwd=4,lty=2,col=NULL)
  
  legend("topleft", legend = c("Start", "End"), 
           col = c("blue", "orange"), lwd=3,cex=1.6,bg="white")
  
  
}

dev.off()
