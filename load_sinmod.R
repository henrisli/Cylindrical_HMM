#https://stackoverflow.com/questions/32183495/netcdf-files-in-r
library(ncdf4)

nc = nc_open("C:/Users/henri/Documents/currentdata/fall.nc")
nc_2 = nc_open("C:/Users/henri/Documents/currentdata/spring.nc")

t_start = 150

current = list()
current$x = ncvar_get(nc, "gridLons")
current$y = ncvar_get(nc, "gridLats")
current$u = ncvar_get(nc, "u_east", start = c(1,1,1,t_start), count = c(-1,-1,1,1))
current$v = ncvar_get(nc, "v_north", start = c(1,1,1,t_start), count = c(-1,-1,1,1))

current_2 = list()
current_2$x = ncvar_get(nc_2, "gridLons")
current_2$y = ncvar_get(nc_2, "gridLats")
current_2$u = ncvar_get(nc_2, "u_east", start = c(1,1,1,t_start), count = c(-1,-1,1,1))
current_2$v = ncvar_get(nc_2, "v_north", start = c(1,1,1,t_start), count = c(-1,-1,1,1))

# x_points = seq(106+36,175+36,3)
# y_points = seq(132,201,3)
x_points = seq(1,70,3)
y_points = seq(123,192,3)

# image(x=x_points, y = y_points, current$u[x_points,y_points])
current$theta = atan(current$v/current$u)
current$theta = ifelse(current$u<0&current$v<0, current$theta-pi, current$theta)
current$theta = ifelse(current$u<0&current$v>0, current$theta + pi, current$theta)
current$r = sqrt(current$u^2+current$v^2)

current_2$theta = atan(current_2$v/current_2$u)
current_2$theta = ifelse(current_2$u<0&current_2$v<0, current_2$theta-pi, current_2$theta)
current_2$theta = ifelse(current_2$u<0&current_2$v>0, current_2$theta + pi, current_2$theta)
current_2$r = sqrt(current_2$u^2+current_2$v^2)

r_sample = as.vector(current$r[x_points, y_points])
theta_sample = as.vector(current$theta[x_points, y_points])
r_sample_2 = as.vector(current_2$r[x_points, y_points])
theta_sample_2 = as.vector(current_2$theta[x_points, y_points])

df <- data.frame(x=as.vector(current$x[x_points,y_points]),y=as.vector(current$y[x_points,y_points]),dx=r_sample*cos(theta_sample)*0.2,dy=r_sample*sin(theta_sample)*0.2)
df_red_arrow = data.frame(x = 2.625, dx=3*0.5*2.5, dy = 0, y = 65.5)
ggplot(data=df) + geom_segment(aes(x=x, y=y, xend=x+dx, yend=y+dy), arrow = arrow(length = unit(0.08,"cm"))) + theme_bw(base_size=20) + xlab("Longitude (°E)") + ylab("Latitude (°N)") + coord_cartesian(xlim = c(3.48,5), ylim = c(61.6,62.26))# + coord_fixed(ratio=1) 
#ggplot(data=df, aes(x=x, y=y, col = as.factor(apply(estimated_probabilities, 1, which.max)))) + geom_segment(aes(xend=x+dx, yend=y+dy), arrow = arrow(length = unit(0.08,"cm"))) + theme_bw(base_size=20) + theme(legend.position="none") + xlab("Longitude (°E)") + ylab("Latitude (°N)") + scale_color_manual(values = colorss) + coord_cartesian(xlim = c(3.48,5), ylim = c(61.6,62.26))
#ggplot(data=df, aes(x=x, y=y, col = as.factor(apply(estimated_probabilities, 1, which.max)), alpha = apply(estimated_probabilities, 1, max))) + geom_segment(aes(xend=x+dx, yend=y+dy), arrow = arrow(length = unit(0.08,"cm"))) + theme_bw(base_size=20) + theme(legend.position="none") + xlab("Longitude (°E)") + ylab("Latitude (°N)") + scale_color_manual(values = colorss) + coord_cartesian(xlim = c(-3-1.875,12+1.875), ylim = c(65.25,72.75))

df_2 <- data.frame(x=as.vector(current_2$x[x_points,y_points]),y=as.vector(current_2$y[x_points,y_points]),dx=r_sample_2*cos(theta_sample_2)*0.2,dy=r_sample_2*sin(theta_sample_2)*0.2)
ggplot(data=df_2, aes(x=x, y=y)) + geom_segment(aes(xend=x+dx, yend=y+dy), arrow = arrow(length = unit(0.08,"cm"))) + theme_bw(base_size=20) + xlab("Longitude (°E)") + ylab("Latitude (°N)") + coord_cartesian(xlim = c(3.48,5), ylim = c(61.6,62.26))# + coord_fixed(ratio=1) 
#ggplot(data=df_2, aes(x=x, y=y, col = as.factor(apply(estimated_probabilities_2, 1, which.max)))) + geom_segment(aes(xend=x+dx, yend=y+dy), arrow = arrow(length = unit(0.08,"cm"))) + theme_bw(base_size=20) + theme(legend.position="none") + xlab("Longitude (°E)") + ylab("Latitude (°N)") + scale_color_manual(values = colorss) + coord_cartesian(xlim = c(3.48,5), ylim = c(61.6,62.26))
#ggplot(data=df_2, aes(x=x, y=y, col = as.factor(apply(estimated_probabilities_2, 1, which.max)), alpha = apply(estimated_probabilities_2, 1, max))) + geom_segment(aes(xend=x+dx, yend=y+dy), arrow = arrow(length = unit(0.08,"cm"))) + theme_bw(base_size=20) + theme(legend.position="none") + xlab("Longitude (°E)") + ylab("Latitude (°N)") + scale_color_manual(values = colorss) + coord_cartesian(xlim = c(-3-1.875,12+1.875), ylim = c(65.25,72.75))


data_set = data.frame(x = r_sample, theta = theta_sample)
ggplot(data_set) + geom_point(aes(x=x, y = theta)) + theme_classic(base_size=20) + xlab("X") + ylab(expression(paste(Phi))) + coord_cartesian(xlim = c(0,0.6))
#ggplot(data_set) + geom_point(aes(x=x, y = theta, col = as.factor(apply(estimated_probabilities, 1, which.max)))) + theme_classic(base_size=20) + theme(legend.position = "none") + xlab("X") + ylab(expression(paste(Phi))) + coord_cartesian(xlim = c(0,0.6)) + scale_color_manual(values = colorss)
#ggplot(data_set) + geom_point(aes(x=x, y = theta, col = as.factor(apply(estimated_probabilities, 1, which.max)), alpha = apply(estimated_probabilities, 1, max))) + theme_classic(base_size=20) + theme(legend.position = "none") + xlab("X") + ylab(expression(paste(Phi))) + coord_cartesian(xlim = c(0,0.5)) + scale_color_manual(values = colorss)

data_set_2 = data.frame(x = r_sample_2, theta = theta_sample_2)
ggplot(data_set_2) + geom_point(aes(x=x, y = theta)) + theme_classic(base_size=20) + xlab("X") + ylab(expression(paste(Phi))) + coord_cartesian(xlim = c(0,0.6))
#ggplot(data_set_2) + geom_point(aes(x=x, y = theta, col = as.factor(apply(estimated_probabilities_2, 1, which.max)))) + theme_classic(base_size=20) + theme(legend.position = "none") + xlab("X") + ylab(expression(paste(Phi))) + coord_cartesian(xlim = c(0,0.6)) + scale_color_manual(values = colorss)
#ggplot(data_set_2) + geom_point(aes(x=x, y = theta, col = as.factor(apply(estimated_probabilities_2, 1, which.max)), alpha = apply(estimated_probabilities_2, 1, max))) + theme_classic(base_size=20) + theme(legend.position = "none") + xlab("X") + ylab(expression(paste(Phi))) + coord_cartesian(xlim = c(0,0.5)) + scale_color_manual(values = colorss)


simulated_sample = as.matrix(data_set)
simulated_sample_2 = as.matrix(data_set_2)



par_est = read.table("C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_sinmod_3.csv")[,1]
estimated_probabilities = find_back_probs(par_est, n_rows, simulated_sample, n_cols)
estimated_probabilities_2 = find_back_probs(par_est, n_rows, simulated_sample_2, n_cols)
# estimated_probabilities = find_back_probs_htlp(par_est, n_rows, simulated_sample, n_cols)
# estimated_probabilities_2 = find_back_probs_htlp(par_est, n_rows, simulated_sample_2, n_cols)


loglik_new = function(param1, x){
  return(-sum(log(apply(x, 1, dabeley_reparam, param=param1))))
}

# WeiSSVM
#param_test = c(1.1,0.9,0,1.7,-0.8)
# X100:
#param_test = c(1.7,0.1,0,0.8,0)
#param_test = c(3,0.5,0,0,0)
# X1: 0.9-3  5-25  -pi - pi  0 - 3  -0.9-0.9
param_test = c(1.4,10,0,0.6,0)
param_test = c(1.9,15,0,0,0)
values = apply(vals, MARGIN= 1, FUN = dabeley, param=param_test)
values[which(values==Inf)]=0
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))

# HTLP
# 0.3-1.2  0.04-0.2  -pi-pi  0-3  0-0.99
param_test = c(0.6,0.1,0,3,0.85)
values = apply(vals, MARGIN= 1, FUN = dhtlp, param=param_test)
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
colorss = ggplotColours(n=6)[c(1,3,5,2)]
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
