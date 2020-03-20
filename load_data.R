#https://stackoverflow.com/questions/32183495/netcdf-files-in-r
library(ncdf4)

#nc = nc_open("C:/Users/henri/Documents/20170107090000-GLOBCURRENT-L4-CUReul_15m-ALT_SUM-v03.0-fv01.0.nc")
#nc = nc_open("C:/Users/henri/Documents/20180104000000-GLOBCURRENT-L4-CUReul_15m-ALT_SUM_NRT-v03.0-fv01.0.nc")
nc = nc_open("C:/Users/henri/Documents/19950106120000-GLOBCURRENT-L4-CUReul_hs-ALT_SUM-v03.0-fv01.0.nc")
nc = nc_open("C:/Users/henri/Documents/19950406120000-GLOBCURRENT-L4-CUReul_hs-ALT_SUM-v03.0-fv01.0.nc")
nc = nc_open("C:/Users/henri/Documents/20050705120000-GLOBCURRENT-L4-CUReul_hs-ALT_SUM-v03.0-fv01.0.nc")
nc = nc_open("C:/Users/henri/Documents/20051003120000-GLOBCURRENT-L4-CUReul_hs-ALT_SUM-v03.0-fv01.0.nc")

current = list()
current$x = ncvar_get(nc, "lon")
current$y = ncvar_get(nc, "lat")
current$u = ncvar_get(nc, "eastward_eulerian_current_velocity")
current$v = ncvar_get(nc, "northward_eulerian_current_velocity")
#current$u = ncvar_get(nc, "u")
#current$v = ncvar_get(nc, "v")
# y_points = which(current$y>=68 & current$y <= 76)
# x_points = which(current$x <= 8 & current$x >= 0)
y_points = which(current$y>=66 & current$y <= 72)
x_points = which(current$x <= 12 & current$x >= -3)[floor(1:24)*5/2]

# image(x=x_points, y = y_points, current$u[x_points,y_points])
current$theta = atan(current$v/current$u)
current$theta = ifelse(current$u<0&current$v<0, current$theta-pi, current$theta)
current$theta = ifelse(current$u<0&current$v>0, current$theta + pi, current$theta)
current$r = sqrt(current$u^2+current$v^2)
# image.plot(x=x_points, y = y_points, current$theta[x_points,y_points])
# image.plot(x=x_points, y = y_points, current$r[x_points,y_points])
r_sample = as.vector(current$r[x_points, y_points])
theta_sample = as.vector(current$theta[x_points, y_points])
df <- data.frame(x=rep(current$x[x_points],length(y_points)),y=rep(current$y[y_points],each=length(x_points)),dx=r_sample*cos(theta_sample)*3,dy=r_sample*sin(theta_sample)*3)
ggplot(data=df, aes(x=x, y=y)) + geom_segment(aes(xend=x+dx, yend=y+dy), arrow = arrow(length = unit(0.08,"cm"))) + theme_bw(base_size=20)# + coord_fixed(ratio=1) + xlab("Longitude") + ylab("Latitude")
#ggplot(data=df, aes(x=x, y=y, col = as.factor(apply(estimated_probabilities, 1, which.max)))) + geom_segment(aes(xend=x+dx, yend=y+dy), arrow = arrow(length = unit(0.08,"cm"))) + theme_bw(base_size=20) + theme(legend.position="none") + coord_fixed(ratio=1) + xlab("Longitude") + ylab("Latitude") + scale_color_manual(values = colorss)


data_set = data.frame(x = r_sample, theta = theta_sample)
ggplot(data_set) + geom_point(aes(x=x, y = theta)) + theme_classic(base_size=20) + xlab("X") + ylab(expression(paste(Phi))) + coord_cartesian(xlim = c(0,0.5))
#ggplot(data_set) + geom_point(aes(x=x, y = theta, col = as.factor(apply(estimated_probabilities, 1, which.max)))) + theme_classic(base_size=20) + theme(legend.position = "none") + xlab("X") + ylab(expression(paste(Phi))) + coord_cartesian(xlim = c(0,0.5)) + scale_color_manual(values = colorss)
simulated_sample = as.matrix(data_set)


par_est = read.table("C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_2015_fall_5.csv")[,1]
estimated_probabilities = find_back_probs(par_est, n_rows, simulated_sample, n_cols)


loglik_new = function(param1, x){
  return(-sum(log(apply(x, 1, dabeley_reparam, param=param1))))
}
#param_test = c(1.1,0.9,0,1.7,-0.8)
# X100:
#param_test = c(1.7,0.1,0,0.8,0)
#param_test = c(3,0.5,0,0,0)
# X1:
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


ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

colorss = ggplotColours(n=3)
colorss = ggplotColours(n=6)[c(1,3,5,2,6)]

# index_2 = 1:576
# for (i in 1:220){
#   swap = sample(1:576,2)
#   index_2[swap] = index_2[swap[2:1]]
#   if(length(unique(index_2))!=576){print(swap)}
# }
# 
# index_2_rev = rep(0,576)
# for (i in 1:576){
#   index_2_rev[i] = which(index_2==i)
# }
