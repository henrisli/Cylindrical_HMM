ggplot()+
  geom_point(aes(x=rep(1,10), y=1:10)) + geom_line(aes(x=rep(1,10), y=1:10))+
  geom_point(aes(x=rep(2,10), y=1:10)) + geom_line(aes(x=rep(2,10), y=1:10))+
  geom_point(aes(x=rep(3,10), y=1:10)) + geom_line(aes(x=rep(3,10), y=1:10))+
  geom_point(aes(x=rep(4,10), y=1:10)) + geom_line(aes(x=rep(4,10), y=1:10))+
  geom_point(aes(x=rep(5,10), y=1:10)) + geom_line(aes(x=rep(5,10), y=1:10))+
  geom_point(aes(x=rep(6,10), y=1:10)) + geom_line(aes(x=rep(6,10), y=1:10))+
  geom_point(aes(x=rep(7,10), y=1:10)) + geom_line(aes(x=rep(7,10), y=1:10))+
  geom_point(aes(x=rep(8,10), y=1:10)) + geom_line(aes(x=rep(8,10), y=1:10))+
  geom_point(aes(x=rep(9,10), y=1:10)) + geom_line(aes(x=rep(9,10), y=1:10))+
  geom_point(aes(x=rep(10,10), y=1:10)) + geom_line(aes(x=rep(10,10), y=1:10))+
  geom_point(aes(x=1:10, y=rep(1,10))) + geom_line(aes(x=1:10, y=rep(1,10)))+
  geom_point(aes(x=1:10, y=rep(2,10))) + geom_line(aes(x=1:10, y=rep(2,10)))+
  geom_point(aes(x=1:10, y=rep(3,10))) + geom_line(aes(x=1:10, y=rep(3,10)))+
  geom_point(aes(x=1:10, y=rep(4,10))) + geom_line(aes(x=1:10, y=rep(4,10)))+
  geom_point(aes(x=1:10, y=rep(5,10))) + geom_line(aes(x=1:10, y=rep(5,10)))+
  geom_point(aes(x=1:10, y=rep(6,10))) + geom_line(aes(x=1:10, y=rep(6,10)))+
  geom_point(aes(x=1:10, y=rep(7,10))) + geom_line(aes(x=1:10, y=rep(7,10)))+
  geom_point(aes(x=1:10, y=rep(8,10))) + geom_line(aes(x=1:10, y=rep(8,10)))+
  geom_point(aes(x=1:10, y=rep(9,10))) + geom_line(aes(x=1:10, y=rep(9,10)))+
  geom_point(aes(x=1:10, y=rep(10,10))) + geom_line(aes(x=1:10, y=rep(10,10)))+
  xlab("x") + ylab("y")

library(ggpubr)
library(grid)
library(reshape2)
library(circular)
library(ggplot2)
library(fields)
n=300
kappa=1
lambda=0
alpha=2
beta=1
mu=0
#u = runif(n,0,2*pi)
#c = 2*tanh(kappa/2)/(1+tanh(kappa/2)^2)
#theta_1 = (acos((cos(u)+c)/(1+c*cos(u))) + mu)%%(2*pi)
theta_1 = rwrappedcauchy(n, mu = circular(mu), rho = tanh(kappa/2))
theta_1 = ifelse(as.numeric(theta_1)>pi, as.numeric(theta_1)-2*pi, as.numeric(theta_1))
#theta_1 = ifelse(theta_1<rep(-pi,n), theta_1+2*pi, theta_1)
#theta_1 = ifelse(theta_1>rep(pi,n), theta_1-2*pi, theta_1)

u = runif(n)
theta = ifelse(u<(1+lambda*sin(theta_1-mu))/2, theta_1, -theta_1)
#for(i in 1:10){
X = rweibull(n, scale = 1/(beta*(1-tanh(kappa)*cos(theta-mu))^(1/alpha)), shape = alpha)

df = data.frame(x=X, y=theta)
ggplot(df) + geom_bin2d(aes(x=x, y=y))


loglik <- function(param, X, theta){
  n = length(X)
  alpha = param[1]
  beta = param[2]
  mu = param[3]
  kappa = param[4]
  lambda = param[5]
  if(alpha<0|beta<0|mu<(-pi)|mu>pi|kappa<0|abs(lambda)>1){return(Inf)}
  return(-((alpha-1)*sum(log(X)) - beta^alpha*sum(X^alpha*(1-tanh(kappa)*cos(theta-mu)))+sum(log(1+lambda*sin(theta-mu)))+n*(alpha*log(beta)+log(alpha)-log(2*pi*cosh(kappa)))))
}

gradient <- function(param){
  alpha = param[1]
  beta = param[2]
  mu = param[3]
  kappa = param[4]
  lambda = param[5]
  dalpha = sum(log(X)) - beta^alpha*sum(log(beta*X)*X^alpha*(1-tanh(kappa)*cos(theta-mu))) + n*(log(beta)+1/alpha)
  dbeta = -alpha*beta^(alpha-1)*sum(X^alpha*(1-tanh(kappa)*cos(theta-mu))) + n*alpha/beta
  dmu = beta^alpha*tanh(kappa)*sum(X^alpha*sin(theta-mu))-lambda*sum(cos(theta-mu)/(1+lambda*sin(theta-mu)))
  dkappa = beta^alpha/cosh(kappa)^2*sum(X^alpha*cos(theta-mu))-n*tanh(kappa)
  dlambda = sum(sin(theta-mu)/(1+lambda*sin(theta-mu)))
  return(c(dalpha,dbeta,dmu,dkappa,dlambda))
}

init_par = rep(0.5,5)
ans = optim(init_par, loglik, X = X, theta = theta, control = list(trace=6, maxit=10000))
ggplot() + geom_point(aes(x=1:5, y =ans$par, col = "Est")) + geom_point(aes(x=1:5, y = c(alpha,beta,mu,kappa,lambda), col = "True"))

abeley = function(alpha, beta, mu, kappa, lambda, X){
  return(alpha*beta^alpha/(2*pi*cosh(kappa))*(1+lambda*sin(X[2]-mu))*X[1]^(alpha-1)*exp(-(beta*X[1])^alpha*(1-tanh(kappa)*cos(X[2]-mu))))}
X_cor = seq(0,5,l=100)
y_cor = seq(-3,3,l=100)
vals = cbind(rep(X_cor,100), rep(y_cor,each=100))
values = apply(X = vals, MARGIN= 1, FUN = abeley, alpha = alpha, beta = beta, mu=mu, kappa=kappa, lambda = lambda)
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))


n_grid = 30
mat_grid <- matrix(seq(n_grid^2), n_grid)

addresses <- expand.grid(x = 1:n_grid, y = 1:n_grid)

# Relative addresses
z <- rbind(c(0,-1,1,0),c(-1,0,0,1))

get.neighbors <- function(rw) {
  # Convert to absolute addresses 
  z2 <- t(z + unlist(rw))
  # Choose those with indices within mat_grid 
  b.good <- rowSums(z2 > 0)==2  &  z2[,1] <= nrow(mat_grid)  &  z2[,2] <=ncol(mat_grid)
  mat_grid[z2[b.good,]]
}

# Not exactly sure if the simulation is working...

neighbor_list = apply(addresses,1, get.neighbors) # Returns a list with neighbors
k=3
ncolor = k

simulate_posterior = function(beta, l_0, iterations, k){
  l = l_0
  n = length(l_0)
  for (j in 1:iterations){
    proposed_site = sample(1:n,1)
    neighbor = neighbor_matrix[[1]][proposed_site]
    for (k_0 in 2:k){neighbor = c(neighbor, neighbor_matrix[[k_0]][proposed_site])}
    #print(neighbor)
    p = exp(beta*neighbor)/sum(exp(beta*neighbor))
    old_i = l[proposed_site]
    new_i = which(rmultinom(1,1,p)==1)
    neighbor_matrix[[old_i]][neighbor_list[[proposed_site]]] = neighbor_matrix[[old_i]][neighbor_list[[proposed_site]]] - 1
    neighbor_matrix[[new_i]][neighbor_list[[proposed_site]]] = neighbor_matrix[[new_i]][neighbor_list[[proposed_site]]] + 1
    l[proposed_site] = new_i
  }
  return(l)
}

mat = matrix(0,ncol = n_grid^2, nrow = 9)
rho = 1.2

for (i in 1:9){
  neighbor_matrix = list(matrix(0,n_grid,n_grid))
  for (mat_count in 2:k){
    neighbor_matrix[[mat_count]] = matrix(0,n_grid,n_grid)}
  init_l = sample(1:k, replace=T, size = n_grid^2)
  for (grid_count in 1:n_grid^2){
    neighbor_matrix[[init_l[grid_count]]][neighbor_list[[grid_count]]] =  neighbor_matrix[[init_l[grid_count]]][neighbor_list[[grid_count]]]+1
  }
  mat[i,] = init_l
  
  mat[i,] = simulate_posterior(rho, mat[i,], n_grid^2*40, k)
}
x = rep(seq(1,n_grid),n_grid)
y = rep(seq(1,n_grid),each = n_grid)
plot_image <- function(l){
  mtrx3d <- data.frame(x = x, y = y, z=as.vector(l))
  mtrx.melt <- melt(mtrx3d, id.vars = c("x","y"), measure.vars = "z")
  mtrx.melt$value = as.factor(mtrx.melt$value)
  #return(ggplot(mtrx.melt, aes(x = x, y = y, fill = value)) +
  #         geom_raster() + coord_fixed(ratio=1)+ scale_fill_manual("Rock type", values=c("dark blue","dark red", "dark green")) +  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_text(size=4), axis.text.x = element_text(size=4)))
  return(ggplot(mtrx.melt, aes(x = x, y = y, fill = value)) +
                    geom_raster() + coord_fixed(ratio=1)+ scale_fill_discrete() +  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_text(size=4), axis.text.x = element_text(size=4)))
}
p1 = plot_image(mat[1,])
p2 = plot_image(mat[2,])
p3 = plot_image(mat[3,])
p4 = plot_image(mat[4,])
p5 = plot_image(mat[5,])
p6 = plot_image(mat[6,])
p7 = plot_image(mat[7,])
p8 = plot_image(mat[8,])
p9 = plot_image(mat[9,])

ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)

##
#set.seed(12345)
rho = 0.4
potts_param <- c(rep(0, ncolor), rho)
x_potts <- matrix(1, nrow = n_grid, ncol = n_grid)
foo <- packPotts(x_potts, ncolor)
out <- potts(foo, potts_param, nbatch = 10)
pdf(file="C:/Users/henri/Documents/GitHub/Master-Thesis/Images/Potts_04.pdf")
image(unpackPotts(out$final), x = 1:24, y = 1:24, xlab = "", ylab = "", col = tim.colors(64))
dev.off()
plot_image(unpackPotts(out$final))
pdf(file="C:/Users/henri/Documents/GitHub/Master-Thesis/Images/Potts_0.pdf")
image(matrix(apply(rmultinom(24*24,1,rep(1/3,3)), 2, function(i) which(i==1)),nrow=24), x = 1:24, y = 1:24, xlab = "", ylab = "", col = tim.colors(64))
dev.off()

spat_pros = as.vector(unpackPotts(out$final))

dabeley <- function(param, x){
  alpha = param[1]
  beta=param[2]
  mu=param[3]
  kappa=param[4]
  lambda=param[5]
  return(alpha*beta^alpha/(2*pi*cosh(kappa))*(1+lambda*sin(x[2]-mu))*x[1]^(alpha-1)*exp(-(beta*x[1])^alpha*(1-tanh(kappa)*cos(x[2]-mu))))}
X_cor = seq(0,5,l=100)
y_cor = seq(-pi,pi,l=100)
vals = cbind(rep(X_cor,100), rep(y_cor,each=100))
parameters = rbind(c(2,1,0,0,1), c(2,1,0,0,-1), c(2,0.6,0,1.5,0))
values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[1,])
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[2,])
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[3,])
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))

simulated_sample = data.frame(x = rep(NA,n_grid*n_grid), theta = rep(NA,n_grid*n_grid))
for(i in 1:(n_grid*n_grid)){
  k_i = spat_pros[i]
  theta_1 = rwrappedcauchy(1, mu = circular(parameters[k_i,3]), rho = tanh(parameters[k_i,4]/2))
  theta_1 = ifelse(as.numeric(theta_1)>pi, as.numeric(theta_1)-2*pi, as.numeric(theta_1))
  
  u = runif(1)
  simulated_sample$theta[i] = ifelse(u<(1+parameters[k_i,5]*sin(theta_1-parameters[k_i,3]))/2, theta_1, -theta_1)
  simulated_sample$x[i] = rweibull(1, scale = 1/(parameters[k_i,2]*(1-tanh(parameters[k_i,4])*cos(simulated_sample$theta[i]-parameters[k_i,3]))^(1/parameters[k_i,1])), shape = parameters[k_i,1])
}
ggplot(simulated_sample) + geom_point(aes(x=x, y = theta, col = as.factor(spat_pros))) + theme_bw() + theme(legend.position = "none")
 
df <- data.frame(x=rep((1:30),30),y=rep((1:30),each=30),dx=simulated_sample[,1]*cos(simulated_sample[,2]),dy=simulated_sample[,1]*sin(simulated_sample[,2]))
ggplot(data=df, aes(x=x, y=y)) + geom_segment(aes(xend=x+dx, yend=y+dy, col = as.factor(spat_pros)), arrow = arrow(length = unit(0.2,"cm"))) + theme(legend.position="none", 
  panel.background = element_rect(fill = "#BFD5E3", colour = "#BFD5E3",
                                  size = 2, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                  colour = "#BFD5E3"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "#BFD5E3")
) + coord_fixed(ratio=1)
plot_image(spat_pros)


#https://stackoverflow.com/questions/32183495/netcdf-files-in-r
nc = nc_open("C:/Users/henri/Documents/20170107090000-GLOBCURRENT-L4-CUReul_15m-ALT_SUM-v03.0-fv01.0.nc")
nc = nc_open("C:/Users/henri/Documents/20180104000000-GLOBCURRENT-L4-CUReul_15m-ALT_SUM_NRT-v03.0-fv01.0.nc")
current = list()
current$x = ncvar_get(nc, "lon")
current$y = ncvar_get(nc, "lat")
current$u = ncvar_get(nc, "eastward_eulerian_current_velocity")
current$v = ncvar_get(nc, "northward_eulerian_current_velocity")
#current$u = ncvar_get(nc, "u")
#current$v = ncvar_get(nc, "v")
y_points = which(current$y>68 & current$y < 76)
x_points = which(current$x < 8 & current$x > 0)
#y_points = which(current$y>72 & current$y < 76)
#x_points = which(current$x < 32 & current$x > 20)
image(x=x_points, y = y_points, current$u[x_points,y_points])
current$theta = atan(current$v/current$u)
current$theta = ifelse(current$u<0&current$v<0, current$theta-pi, current$theta)
current$theta = ifelse(current$u<0&current$v>0, current$theta + pi, current$theta)
current$r = sqrt(current$u^2+current$v^2)
image.plot(x=x_points, y = y_points, current$theta[x_points,y_points])
image.plot(x=x_points, y = y_points, current$r[x_points,y_points])
r_sample = as.vector(current$r[x_points, y_points])
theta_sample = as.vector(current$theta[x_points, y_points])
df <- data.frame(x=rep(current$x[x_points],length(y_points)),y=rep(current$y[y_points],each=length(x_points)),dx=r_sample*cos(theta_sample)*3,dy=r_sample*sin(theta_sample)*3)
ggplot(data=df, aes(x=x, y=y)) + geom_segment(aes(xend=x+dx, yend=y+dy), col = as.factor(apply(xi_probs_i_normal, 1, which.max)), arrow = arrow(length = unit(0.08,"cm"))) + theme(legend.position="none") + coord_fixed(ratio=1)


data_set = data.frame(x = r_sample, theta = theta_sample)
ggplot(data_set) + geom_point(aes(x=x, y = theta), col = as.factor(apply(xi_probs_i_normal, 1, which.max)), size = 0.1) + theme_classic()# + coord_cartesian(xlim = c(0,50))
loglik_new = function(param1, x){
  return(-sum(log(apply(x, 1, dabeley_reparam, param=param1))))
}
#param_test = c(1.1,0.9,0,1.7,-0.8)
#param_test = c(1.5,0.2,0,0.8,0)
param_test = c(1.4,10,0,0.6,0)
param_test = c(1.9,15,0,0,0)
values = apply(vals, MARGIN= 1, FUN = dabeley, param=param_test)
values[which(values==Inf)]=0
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))

for (i in 1:(length(x_points)*length(y_points))){
  theta_1 = rwrappedcauchy(1, mu = circular(param_test[3]), rho = tanh(param_test[4]/2))
  theta_1 = ifelse(as.numeric(theta_1)>pi,as.numeric(theta_1)-2*pi, as.numeric(theta_1))
  u = runif(1)
  data_set$theta[i] = ifelse(u<(1+param_test[5]*sin(theta_1-param_test[3]))/2, theta_1, -theta_1)
  data_set$x[i] = rweibull(1, scale = 1/(param_test[2]*(1-tanh(param_test[4])*cos(data_set$theta[i]-param_test[3]))^(1/param_test[1])), shape = param_test[1])
}
ggplot(data_set) + geom_point(aes(x=x, y = theta)) + theme_classic()# + coord_cartesian(xlim = c(0,1))

param_test_reparam = rep(0,5)
param_test_reparam[c(1,2,4)] = log(param_test[c(1,2,4)])
param_test_reparam[3] = atan(param_test[3]/2)
param_test_reparam[5] = atanh(param_test[5])
init_param = param_test_reparam


#optim_test = optim(param_test, loglik, X= data_set$x, theta = data_set$theta, control = list(trace = 6, REPORT = 1))
optim_test = optim(init_param, loglik_new, x = as.matrix(data_set), control = list(trace = 6, REPORT = 1), method = "BFGS")
param_est = c(exp(optim_test$par[1:2]), 2*tan(optim_test$par[3]), exp(optim_test$par[4]), tanh(optim_test$par[5]))
values = apply(vals, MARGIN= 1, FUN = dabeley, param=param_est)
values[which(values==Inf)]=0
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))


#theta = ifelse(u<0&v<0, theta-pi, theta)
#theta = ifelse(u<0&v>0, theta + pi, theta)
#r_sample = as.vector(r[711:730,651:670])
#theta_sample = as.vector(theta[711:730,651:670])
#df <- data.frame(x=rep((1:20),20)*0.05,y=rep((1:20),each=20)*0.05,dx=r_sample*cos(theta_sample),dy=r_sample*sin(theta_sample))
#ggplot(data=df, aes(x=x, y=y)) + geom_segment(aes(xend=x+dx, yend=y+dy), arrow = arrow(length = unit(0.2,"cm"))) + theme(legend.position="none") + coord_fixed(ratio=1)


library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
pdf(file="C:/Users/henri/Documents/GitHub/Master-Thesis/Images/test_image2.pdf")
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(20-16, 32+16), ylim = c(72-4, 76+4), expand = FALSE) + scale_x_continuous(breaks = seq(-32,40,8)) + scale_y_continuous(breaks = seq(64,80,2)) + geom_segment(data = df, aes(x = x, y = y, xend=x+dx, yend=y+dy), arrow = arrow(length = unit(0.1,"cm"))) + xlab("") + ylab("")
dev.off()
