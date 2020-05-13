par_est = read.table("C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_summer_2_htlp.csv")[,1]
par_est = read.table("C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_single_2.csv")[,1]
par_est = read.table("C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_sinmod_3.csv")[,1]

parameters = par_est

# parameters[c(2,3,4,5,8,9)] = exp(parameters[c(2,3,4,5,8,9)])
# parameters[c(6,7)] = 2*atan(parameters[c(6,7)])
# parameters[c(10,11)] = tanh(parameters[c(10,11)])

parameters[c(2,3,4,5,6,7,11,12,13)] = exp(parameters[c(2,3,4,5,6,7,11,12,13)])
parameters[c(8,9,10)] = 2*atan(parameters[c(8,9,10)])
parameters[c(14,15,16)] = tanh(parameters[c(14,15,16)])

# parameters[c(2,3,4,5,8,9)] = exp(parameters[c(2,3,4,5,8,9)])
# parameters[c(6,7)] = 2*atan(parameters[c(6,7)])
# parameters[c(10,11)] = exp(parameters[c(10,11)])/(1+exp(parameters[c(10,11)]))

# parameters[c(2,3,4,5,6,7,11,12,13)] = exp(parameters[c(2,3,4,5,6,7,11,12,13)])
# parameters[c(8,9,10)] = 2*atan(parameters[c(8,9,10)])
# parameters[c(14,15,16)] = exp(parameters[c(14,15,16)])/(1+exp(parameters[c(14,15,16)]))

parameters = matrix(parameters[2:(5*ncolor_test+1)],nrow=ncolor_test)
parameters

estimated_probabilities_hold_out = matrix(NA, nrow = 24*24, ncol = ncolor_test)
estimated_probabilities_hold_out_2 = matrix(NA, nrow = 24*24, ncol = ncolor_test)
for (i in 1:(24*24)){
  for (j in 1:ncolor_test){
    estimated_probabilities_hold_out[i,j] = estimated_probabilities[i,j]/dabeley(parameters[j,], as.vector(simulated_sample[i,]))
    estimated_probabilities_hold_out_2[i,j] = estimated_probabilities_2[i,j]/dabeley(parameters[j,], as.vector(simulated_sample_2[i,]))
  }
  estimated_probabilities_hold_out[i,] = estimated_probabilities_hold_out[i,]/sum(estimated_probabilities_hold_out[i,])
  estimated_probabilities_hold_out_2[i,] = estimated_probabilities_hold_out_2[i,]/sum(estimated_probabilities_hold_out_2[i,])
}

observation_density = function(param, x, hold_out_prob){
  sum = 0
  for (i in 1:ncolor_test){
    sum = sum + dabeley(param[i,], x)*hold_out_prob[i]
  }
  return(sum)
}



num_MCMC = 10000
pred_sample_1 = matrix(NA, ncol = 2, nrow = num_MCMC)
pred_sample_2 = matrix(NA, ncol = 2, nrow = num_MCMC)
energy_scores_1 = rep(0, 50)
energy_scores_2 = rep(0, 50)

sample_points = sample(1:(24*24), 50, replace = F)
for (j in sample_points){
  for (i in 1:num_MCMC){
    k_1 = which(rmultinom(1,1,estimated_probabilities_hold_out[j,])==1)
    k_2 = which(rmultinom(1,1,estimated_probabilities_hold_out_2[j,])==1)
  }
  
}


for (i in 1:10){
  j = sample(1:(24*24), 1)
  print(j)
  values = apply(vals, MARGIN= 1, FUN = observation_density, param=parameters, hold_out_prob = estimated_probabilities_hold_out[j,])
  par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
  contour(x=X_cor, y = y_cor, z = matrix(values,nrow=100), xlab = "x", ylab = expression(paste(phi)), labcex = 0.01, cex.lab = 1.7)
  points(x = simulated_sample[j,1], y = simulated_sample[j,2], pch=16, cex = 1.4)
  
}

calculate_expectation = function(k){
  
  legendreFUN = function(x){
    return(1/(cosh(parameters[k,4])-sqrt(cosh(parameters[k,4])^2-1)*cos(x))^(1+1/parameters[k,1]))
  }
  
  legendreVAL = 1/pi*integral(legendreFUN, -pi, 0)
  
  expectation = (cosh(parameters[k,4])/parameters[k,2]^parameters[k,1])^(1/parameters[k,1])*gamma((1+parameters[k,1])/parameters[k,1])*legendreVAL
  return(expectation)
}

expectations_classes = sapply(1:ncolor_test, calculate_expectation)

expectations = estimated_probabilities_hold_out%*%expectations_classes

test_df = data.frame(true = simulated_sample[,1], pred = expectations)

ggplot(test_df) + geom_point(aes(x=1:576, y = true, col = "true")) + geom_point(aes(x=1:576, y = pred, col = "pred"))


X_cor = seq(0,2,l=1000)
y_cor = seq(-pi,pi,l=1000)
vals = cbind(rep(X_cor,1000), rep(y_cor,each=1000))
# X_cor = seq(0,0.6,l=100)
# y_cor = seq(-pi,pi,l=100)
# vals = cbind(rep(X_cor,100), rep(y_cor,each=100))

values = apply(vals, MARGIN= 1, FUN = observation_density, param=parameters, hold_out_prob = estimated_probabilities_hold_out[j,])
par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=1000), xlab = "x", ylab = expression(paste(phi)), cex.lab = 1.7, nlevels = 6, drawlabels = F, method = "flattest")
points(x = simulated_sample[j,1], y = simulated_sample[j,2], pch=16, cex = 1.4)
