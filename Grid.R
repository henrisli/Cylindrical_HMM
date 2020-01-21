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



n=1000000
kappa=1
lambda=1
alpha=2
beta=1
mu=0
theta_1 = rwrappedcauchy(n, mu = mu, rho = tanh(kappa/2))
u = runif(n)
theta = ifelse(u<(1+lambda*sin(theta_1-mu))/2, theta_1, -theta_1)
#for(i in 1:10){
X = rweibull(n, shape = beta*(1-tanh(kappa)*cos(theta-mu))^(1/alpha))

df = data.frame(x=X, y=theta)
ggplot(df) + geom_bin2d(aes(x=x, y=y))





abeley = function(alpha, beta, mu, kappa, lambda, X){
  return(alpha*beta^alpha/(2*pi*cosh(kappa))*(1+lambda*sin(X[2]-mu))*X[1]^(alpha-1)*exp(-(beta*X[1])^alpha*(1-tanh(kappa)*cos(X[2]-mu))))}
X_cor = seq(0,5,l=100)
y_cor = seq(-3,3,l=100)
vals = cbind(rep(X_cor,100), rep(y_cor,each=100))
values = apply(X = vals, MARGIN= 1, FUN = abeley, alpha = alpha, beta = beta, mu=mu, kappa=kappa, lambda = lambda)
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
