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
