library(ggplot2)
library(geoR)
library(fields)
library(akima)
library(metR)

dhtlp <- function(param, x){
  alpha = param[1]
  beta=param[2]
  mu=param[3]
  tau=param[4]
  kappa=param[5]
  if (tau==0){
    return(sqrt(1-kappa^2)/(2*pi*beta*alpha)*(x[1]/beta)^(1/alpha-1)*exp(-(x[1]/beta)^(1/alpha)*(1-kappa*cos(x[2]-mu))))
  }else{
    return(sqrt(1-kappa^2)/(2*pi*beta*alpha)*(x[1]/beta)^(1/alpha-1)*(1+tau/alpha*(x[1]/beta)^(1/alpha)*(1-kappa*cos(x[2]-mu)))^(-alpha/tau-1))
  }
}
X_cor = seq(0,5,l=100)
y_cor = seq(-pi,pi,l=100)
vals = cbind(rep(X_cor,100), rep(y_cor,each=100))
parameters = rbind(c(0.25,1,0,0,0.75), c(0.25,1,0,0.3,0.75), c(0.25,1,0,0.6,0.75), c(0.5,1,0,0,0.75), c(0.5,1,0,0.3,0.75), c(0.5,1,0,0.6,0.75), c(0.75,1,0,0,0.75), c(0.75,1,0,0.3,0.75), c(0.75,1,0,0.6,0.75))
values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[1,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[2,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[3,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[4,])
#df = data.frame(x=vals[,1], y = vals[,2], z = values)
#ggplot(df) + stat_contour(aes(x = x, y = y, z = z)) + theme_classic() + geom_text_contour(aes(x=x, y = y, z=z), stroke.color = "black", stroke = 0.2) + xlab("x") + ylab(expression(paste(phi)))
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[5,])
#df = data.frame(x=vals[,1], y = vals[,2], z = values)
#ggplot(df) + stat_contour(aes(x = x, y = y, z = z)) + theme_classic() + geom_text_contour(aes(x=x, y = y, z=z), stroke.color = "black", stroke = 0.2) + xlab("x") + ylab(expression(paste(phi)))
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[6,])
#df = data.frame(x=vals[,1], y = vals[,2], z = values)
#ggplot(df) + stat_contour(aes(x = x, y = y, z = z)) + theme_classic() + geom_text_contour(aes(x=x, y = y, z=z), stroke.color = "black", stroke = 0.2) + xlab("x") + ylab(expression(paste(phi)))
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[7,])
#df = data.frame(x=vals[,1], y = vals[,2], z = values)
#ggplot(df) + stat_contour(aes(x = x, y = y, z = z)) + theme_classic() + geom_text_contour(aes(x=x, y = y, z=z), stroke.color = "black", stroke = 0.2) + xlab("x") + ylab(expression(paste(phi)))
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[8,])
#df = data.frame(x=vals[,1], y = vals[,2], z = values)
#ggplot(df) + stat_contour(aes(x = x, y = y, z = z)) + theme_classic() + geom_text_contour(aes(x=x, y = y, z=z), stroke.color = "black", stroke = 0.2) + xlab("x") + ylab(expression(paste(phi)))
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[9,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

# Plot true densities from simulated samples

parameters = rbind(c(0.25,1,-pi/2,0,0.7), c(0.25,1,pi/2,0,0.7), c(0.5,0.5,0,0.8,0.8))
values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[1,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[2,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[3,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

parameters = rbind(c(0.25,1,0,0,0.6), c(0.25,3,0,0,0.6), c(0.25,0.5,0,0,0.2))
values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[1,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[2,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)

values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[3,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 1, cex.lab = 1.7)




# Contour of fitted models to real data and points
X_cor = seq(0,0.5,l=100)
y_cor = seq(-pi,pi,l=100)
vals = cbind(rep(X_cor,100), rep(y_cor,each=100))
parameters_test_reparam = read.table("C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_summer_2_htlp.csv")[,1]

parameters = parameters_test_reparam

parameters[c(2,3,4,5,8,9)] = exp(parameters[c(2,3,4,5,8,9)])
parameters[c(6,7)] = 2*atan(parameters[c(6,7)])
parameters[c(10,11)] = exp(parameters[c(10,11)])/(1+exp(parameters[c(10,11)]))

# parameters[c(2,3,4,5,6,7,11,12,13)] = exp(parameters[c(2,3,4,5,6,7,11,12,13)])
# parameters[c(8,9,10)] = 2*atan(parameters[c(8,9,10)])
# parameters[c(14,15,16)] = exp(parameters[c(14,15,16)])/(1+exp(parameters[c(14,15,16)]))

parameters = matrix(parameters[2:(5*ncolor_test+1)],nrow=ncolor_test)
parameters
values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[1,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 0.01, cex.lab = 1.7)
points(x = simulated_sample[,1], y = simulated_sample[,2], col = rgb(0,0,0,alpha=estimated_probabilities[,1]), pch=16, cex = 1)
points(x = simulated_sample_2[,1], y = simulated_sample_2[,2], col = rgb(0,0,0,alpha=estimated_probabilities_2[,1]), pch=15, cex = 1)


values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[2,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 0.01, cex.lab = 1.7)
points(x = simulated_sample[,1], y = simulated_sample[,2], col = rgb(0,0,0,alpha=estimated_probabilities[,2]), pch=16, cex = 1)
points(x = simulated_sample_2[,1], y = simulated_sample_2[,2], col = rgb(0,0,0,alpha=estimated_probabilities_2[,2]), pch=15, cex = 1)


values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[3,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 0.01, cex.lab = 1.7)
points(x = simulated_sample[,1], y = simulated_sample[,2], col = rgb(0,0,0,alpha=estimated_probabilities[,3]), pch=16, cex = 1)
points(x = simulated_sample_2[,1], y = simulated_sample_2[,2], col = rgb(0,0,0,alpha=estimated_probabilities_2[,3]), pch=15, cex = 1)



# Windrose and sample
parameters = c(0.5,1,0,0,0.75)
#parameters = c(0.25,1,0,0.3,0.75)
simulated_sample = data.frame(x = rep(NA,1000), theta = rep(NA,1000))
for(i in 1:(1000)){
  theta_1 = rwrappedcauchy(1, mu = circular(parameters[3]), rho = parameters[5]/(1+sqrt(1-parameters[5]^2)))
  simulated_sample$theta[i] = ifelse(as.numeric(theta_1)>pi, as.numeric(theta_1)-2*pi, as.numeric(theta_1))
  
  u = runif(1)
  if (parameters[4]==0){
    simulated_sample$x[i] = parameters[2]*(-log(1-u)/(1-parameters[5]*cos(simulated_sample$theta[i] - parameters[3])))^parameters[1]
  }else{
    simulated_sample$x[i] = parameters[2]*(((1-u)^(-parameters[4]/parameters[1])-1)/(parameters[4]/parameters[1]*(1-parameters[5]*cos(simulated_sample$theta[i] - parameters[3]))))^parameters[1]
  }
}
ggplot(simulated_sample) + geom_point(aes(x=x, y = theta)) + theme_classic(base_size=20) + theme(legend.position = "none") + xlab("X") + ylab(expression(paste(Phi))) + coord_cartesian(xlim = c(0,10), ylim = c(-pi,pi))
simulated_sample = as.matrix(simulated_sample)

df <- data.frame(x=rep(0,1000),y=rep(0,1000),dx=simulated_sample[,1]*cos(simulated_sample[,2]),dy=simulated_sample[,1]*sin(simulated_sample[,2]))
ggplot(data=df) + geom_segment(aes(x=x, y=y, xend=x+dx, yend=y+dy), arrow = arrow(length = unit(0.15,"cm"))) + theme_classic(base_size=20) + xlab("u") + ylab("v") + coord_cartesian(xlim = c(-7,7), ylim = c(-7,7)) 
