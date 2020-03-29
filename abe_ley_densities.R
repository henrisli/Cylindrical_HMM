library(ggplot2)
library(geoR)
library(fields)
library(akima)
library(metR)

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
parameters = rbind(c(2,1,0,0,0), c(2,1,0,0,0.5), c(2,1,0,0,1), c(2,1,0,1,0), c(2,1,0,1,0.5), c(2,1,0,1,1), c(2,1,0,1.5,0), c(2,1,0,1.5,0.5), c(2,1,0,1.5,1))
values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[1,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[2,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[3,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[4,])
#df = data.frame(x=vals[,1], y = vals[,2], z = values)
#ggplot(df) + stat_contour(aes(x = x, y = y, z = z)) + theme_classic() + geom_text_contour(aes(x=x, y = y, z=z), stroke.color = "black", stroke = 0.2) + xlab("x") + ylab(expression(paste(phi)))
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[5,])
#df = data.frame(x=vals[,1], y = vals[,2], z = values)
#ggplot(df) + stat_contour(aes(x = x, y = y, z = z)) + theme_classic() + geom_text_contour(aes(x=x, y = y, z=z), stroke.color = "black", stroke = 0.2) + xlab("x") + ylab(expression(paste(phi)))
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[6,])
#df = data.frame(x=vals[,1], y = vals[,2], z = values)
#ggplot(df) + stat_contour(aes(x = x, y = y, z = z)) + theme_classic() + geom_text_contour(aes(x=x, y = y, z=z), stroke.color = "black", stroke = 0.2) + xlab("x") + ylab(expression(paste(phi)))
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[7,])
#df = data.frame(x=vals[,1], y = vals[,2], z = values)
#ggplot(df) + stat_contour(aes(x = x, y = y, z = z)) + theme_classic() + geom_text_contour(aes(x=x, y = y, z=z), stroke.color = "black", stroke = 0.2) + xlab("x") + ylab(expression(paste(phi)))
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[8,])
#df = data.frame(x=vals[,1], y = vals[,2], z = values)
#ggplot(df) + stat_contour(aes(x = x, y = y, z = z)) + theme_classic() + geom_text_contour(aes(x=x, y = y, z=z), stroke.color = "black", stroke = 0.2) + xlab("x") + ylab(expression(paste(phi)))
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[9,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

# Plot true densities from simulated samples
parameters = rbind(c(2,1,0,0,1), c(2,1,0,0,-1), c(2,0.6,0,1.5,0))
values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[1,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[2,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[3,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

parameters = rbind(c(3,1,0,0.21,0.8), c(5,5,0,0.21,0), c(1,0.8,0,1.7,-0.8))
values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[1,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[2,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[3,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)




# Contour of fitted models to real data and points
X_cor = seq(0,0.5,l=100)
y_cor = seq(-pi,pi,l=100)
vals = cbind(rep(X_cor,100), rep(y_cor,each=100))
parameters_test_reparam = read.table("C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_summer_2.csv")[,1]

parameters = parameters_test_reparam

parameters[c(2,3,4,5,8,9)] = exp(parameters[c(2,3,4,5,8,9)])
parameters[c(6,7)] = 2*atan(parameters[c(6,7)])
parameters[c(10,11)] = tanh(parameters[c(10,11)])

# parameters[c(2,3,4,5,6,7,11,12,13)] = exp(parameters[c(2,3,4,5,6,7,11,12,13)])
# parameters[c(8,9,10)] = 2*atan(parameters[c(8,9,10)])
# parameters[c(14,15,16)] = tanh(parameters[c(14,15,16)])

parameters = matrix(parameters[2:(5*ncolor_test+1)],nrow=ncolor_test)
parameters
values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[1,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 0.01, cex.lab = 1.7)
points(x = simulated_sample[,1], y = simulated_sample[,2], col = rgb(0,0,0,alpha=estimated_probabilities[,1]), pch=16, cex = 1)
points(x = simulated_sample_2[,1], y = simulated_sample_2[,2], col = rgb(0,0,0,alpha=estimated_probabilities_2[,1]), pch=15, cex = 1)


values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[2,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 0.01, cex.lab = 1.7)
points(x = simulated_sample[,1], y = simulated_sample[,2], col = rgb(0,0,0,alpha=estimated_probabilities[,2]), pch=16, cex = 1)
points(x = simulated_sample_2[,1], y = simulated_sample_2[,2], col = rgb(0,0,0,alpha=estimated_probabilities_2[,2]), pch=15, cex = 1)


values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[3,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 0.01, cex.lab = 1.7)
points(x = simulated_sample[,1], y = simulated_sample[,2], col = rgb(0,0,0,alpha=estimated_probabilities[,3]), pch=16, cex = 1)
points(x = simulated_sample_2[,1], y = simulated_sample_2[,2], col = rgb(0,0,0,alpha=estimated_probabilities_2[,3]), pch=15, cex = 1)




# Windrose and sample
parameters = c(2,1,0,1,1)
simulated_sample = data.frame(x = rep(NA,1000), theta = rep(NA,1000))
for(i in 1:(1000)){
  theta_1 = rwrappedcauchy(1, mu = circular(parameters[3]), rho = tanh(parameters[4]/2))
  theta_1 = ifelse(as.numeric(theta_1)>pi, as.numeric(theta_1)-2*pi, as.numeric(theta_1))
  
  u = runif(1)
  simulated_sample$theta[i] = ifelse(u<(1+parameters[5]*sin(theta_1-parameters[3]))/2, theta_1, -theta_1)
  simulated_sample$x[i] = rweibull(1, scale = 1/(parameters[2]*(1-tanh(parameters[4])*cos(simulated_sample$theta[i]-parameters[3]))^(1/parameters[1])), shape = parameters[1])
}
ggplot(simulated_sample) + geom_point(aes(x=x, y = theta)) + theme_classic(base_size=20) + theme(legend.position = "none") + xlab("X") + ylab(expression(paste(Phi))) + coord_cartesian(xlim = c(0,5), ylim = c(-pi,pi))
simulated_sample = as.matrix(simulated_sample)

df <- data.frame(x=rep(0,1000),y=rep(0,1000),dx=simulated_sample[,1]*cos(simulated_sample[,2]),dy=simulated_sample[,1]*sin(simulated_sample[,2]))
ggplot(data=df) + geom_segment(aes(x=x, y=y, xend=x+dx, yend=y+dy), arrow = arrow(length = unit(0.15,"cm"))) + theme_classic(base_size=20) + xlab("u") + ylab("v") + coord_cartesian(xlim = c(-5,5), ylim = c(-5,5)) 
