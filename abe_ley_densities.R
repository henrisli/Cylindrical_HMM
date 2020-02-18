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
