par_est = read.table("C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_summer_2.csv")[,1]
# par_est = read.table("C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_single_2.csv")[,1]
# par_est = read.table("C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_sinmod_3.csv")[,1]

parameters = par_est

parameters[c(2,3,4,5,8,9)] = exp(parameters[c(2,3,4,5,8,9)])
parameters[c(6,7)] = 2*atan(parameters[c(6,7)])
parameters[c(10,11)] = tanh(parameters[c(10,11)])

# parameters[c(2,3,4,5,6,7,11,12,13)] = exp(parameters[c(2,3,4,5,6,7,11,12,13)])
# parameters[c(8,9,10)] = 2*atan(parameters[c(8,9,10)])
# parameters[c(14,15,16)] = tanh(parameters[c(14,15,16)])

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



num_MCMC = 1000
pred_sample_1 = matrix(NA, ncol = 2, nrow = num_MCMC)
pred_sample_euclid_1 = pred_sample_1
pred_sample_2 = matrix(NA, ncol = 2, nrow = num_MCMC)
pred_sample_euclid_2 = pred_sample_2
energy_scores_1_u = rep(0, 50)
energy_scores_2_u = rep(0, 50)

energy_scores_1_v = rep(0,50)
energy_scores_2_v = rep(0,50)

sample_points = read.table("C://Users//henri//Documents//GitHub//Master-Thesis//Data//sample_points.csv")[,1]

for (num_iter in 1:50){
  j = sample_points[num_iter]
  obs_1 = c(simulated_sample[j,1]*cos(simulated_sample[j,2]), simulated_sample[j,1]*sin(simulated_sample[j,2]))
  obs_2 = c(simulated_sample_2[j,1]*cos(simulated_sample_2[j,2]), simulated_sample_2[j,1]*sin(simulated_sample_2[j,2]))
  for (i in 1:num_MCMC){
    k_1 = which(rmultinom(1,1,estimated_probabilities_hold_out[j,])==1)
    k_2 = which(rmultinom(1,1,estimated_probabilities_hold_out_2[j,])==1)
    
    theta_1 = rwrappedcauchy(1, mu = circular(parameters[k_1,3]), rho = tanh(parameters[k_1,4]/2))
    theta_1 = ifelse(as.numeric(theta_1)>pi, as.numeric(theta_1)-2*pi, as.numeric(theta_1))
    
    u = runif(1)
    pred_sample_1[i,2] = ifelse(u<(1+parameters[k_1,5]*sin(theta_1-parameters[k_1,3]))/2, theta_1, -theta_1)
    pred_sample_1[i,1] = rweibull(1, scale = 1/(parameters[k_1,2]*(1-tanh(parameters[k_1,4])*cos(pred_sample_1[i,2]-parameters[k_1,3]))^(1/parameters[k_1,1])), shape = parameters[k_1,1])
    pred_sample_euclid_1[i,1] = pred_sample_1[i,1]*cos(pred_sample_1[i,2])
    pred_sample_euclid_1[i,2] = pred_sample_1[i,1]*sin(pred_sample_1[i,2])
    
    
    theta_2 = rwrappedcauchy(1, mu = circular(parameters[k_2,3]), rho = tanh(parameters[k_2,4]/2))
    theta_2 = ifelse(as.numeric(theta_2)>pi, as.numeric(theta_2)-2*pi, as.numeric(theta_2))
    
    u = runif(1)
    pred_sample_2[i,2] = ifelse(u<(1+parameters[k_2,5]*sin(theta_2-parameters[k_2,3]))/2, theta_2, -theta_2)
    pred_sample_2[i,1] = rweibull(1, scale = 1/(parameters[k_2,2]*(1-tanh(parameters[k_2,4])*cos(pred_sample_2[i,2]-parameters[k_2,3]))^(1/parameters[k_2,1])), shape = parameters[k_2,1])
    pred_sample_euclid_2[i,1] = pred_sample_2[i,1]*cos(pred_sample_2[i,2])
    pred_sample_euclid_2[i,2] = pred_sample_2[i,1]*sin(pred_sample_2[i,2])
    
  }
  
  # u
  sum_a = 0
  for (i in 1:num_MCMC){
    sum_a = sum_a + 1/num_MCMC*abs(pred_sample_euclid_1[i,1] - obs_1[1])
  }
  
  
  sum_aa = 0
  for (num_1 in 1:num_MCMC){
    for (num_2 in 1:num_MCMC){
      sum_aa = sum_aa + 1/num_MCMC^2*abs(pred_sample_euclid_1[num_1,1] - pred_sample_euclid_1[num_2,1])
    }
  }
  
  energy_scores_1_u[num_iter] = sum_a - sum_aa/2
  
  sum_a = 0
  for (i in 1:num_MCMC){
    sum_a = sum_a + 1/num_MCMC*abs(pred_sample_euclid_2[i,1] - obs_2[1])
  }
  
  
  sum_aa = 0
  for (num_1 in 1:num_MCMC){
    for (num_2 in 1:num_MCMC){
      sum_aa = sum_aa + 1/num_MCMC^2*abs(pred_sample_euclid_2[num_1,1] - pred_sample_euclid_2[num_2,1])
    }
  }
  
  energy_scores_2_u[num_iter] = sum_a - sum_aa/2
  
  # v
  sum_a = 0
  for (i in 1:num_MCMC){
    sum_a = sum_a + 1/num_MCMC*abs(pred_sample_euclid_1[i,2] - obs_1[2])
  }
  
  
  sum_aa = 0
  for (num_1 in 1:num_MCMC){
    for (num_2 in 1:num_MCMC){
      sum_aa = sum_aa + 1/num_MCMC^2*abs(pred_sample_euclid_1[num_1,2] - pred_sample_euclid_1[num_2,2])
    }
  }
  
  energy_scores_1_v[num_iter] = sum_a - sum_aa/2
  
  sum_a = 0
  for (i in 1:num_MCMC){
    sum_a = sum_a + 1/num_MCMC*abs(pred_sample_euclid_2[i,2] - obs_2[2])
  }
  
  
  sum_aa = 0
  for (num_1 in 1:num_MCMC){
    for (num_2 in 1:num_MCMC){
      sum_aa = sum_aa + 1/num_MCMC^2*abs(pred_sample_euclid_2[num_1,2] - pred_sample_euclid_2[num_2,2])
    }
  }
  
  energy_scores_2_v[num_iter] = sum_a - sum_aa/2
  
  print(num_iter)
}

mean(energy_scores_1_u)
mean(energy_scores_1_v)
mean(energy_scores_2_u)
mean(energy_scores_2_v)



num_grid = 26
x_grid = seq(-1.5,1.5, l = num_grid)
y_grid = seq(-1.5,1.5, l = num_grid)
grid_size = x_grid[2]-x_grid[1]
cartesian_grid_vals = cbind(rep(x_grid,num_grid), rep(y_grid,each=num_grid))
polar_grid_vals = cartesian_grid_vals
polar_grid_vals[,2] = atan(cartesian_grid_vals[,2]/cartesian_grid_vals[,1])
polar_grid_vals[,2] = ifelse(cartesian_grid_vals[,1]<0&cartesian_grid_vals[,2]<0, polar_grid_vals[,2]-pi, polar_grid_vals[,2])
polar_grid_vals[,2] = ifelse(cartesian_grid_vals[,1]<0&cartesian_grid_vals[,2]>0, polar_grid_vals[,2] + pi, polar_grid_vals[,2])
polar_grid_vals[,1] = sqrt(cartesian_grid_vals[,1]^2+cartesian_grid_vals[,2]^2)

for (num_iter in 1:50){
  j = sample_points[num_iter]
  obs_1 = c(simulated_sample[j,1]*cos(simulated_sample[j,2]), simulated_sample[j,1]*sin(simulated_sample[j,2]))
  obs_2 = c(simulated_sample_2[j,1]*cos(simulated_sample_2[j,2]), simulated_sample_2[j,1]*sin(simulated_sample_2[j,2]))
  
  # u
  sum_a = 0
  for (i in 1:num_grid^2){
    sum_a = sum_a + grid_size^2*abs(caresian_grid_vals[i,1] - obs_1[1])*observation_density(parameters, polar_grid_vals[i,], estimated_probabilities_hold_out[j,])
  }
  
  
  sum_aa = 0
  for (num_1 in 1:num_grid^2){
    for (num_2 in 1:num_grid^2){
      sum_aa = sum_aa + grid_size^4*abs(cartesian_grid_vals[num_1,1] - cartesian_grid_vals[num_2,1])*observation_density(parameters, polar_grid_vals[num_1,], estimated_probabilities_hold_out[j,])*observation_density(parameters, polar_grid_vals[num_2,], estimated_probabilities_hold_out[j,])
    }
  }
  
  energy_scores_1_u[num_iter] = sum_a - sum_aa/2
  
  sum_a = 0
  for (i in 1:num_grid^2){
    sum_a = sum_a + grid_size^2*abs(cartesian_grid_vals[i,1] - obse_2[1])*observation_density(parameters, polar_grid_vals[i,], estimated_probabilities_hold_out_2[j,])
  }
  
  
  sum_aa = 0
  for (num_1 in 1:num_grid^2){
    for (num_2 in 1:num_grid^2){
      sum_aa = sum_aa + grid_size^4*abs(cartesian_grid_vals[num_1,1] - cartesian_grid_vals[num_2,1])*observation_density(parameters, polar_grid_vals[num_1,], estimated_probabilities_hold_out_2[j,])*observation_density(parameters, polar_grid_vals[num_2,], estimated_probabilities_hold_out_2[j,])
    }
  }
  
  energy_scores_2_u[num_iter] = sum_a - sum_aa/2
  
  # v
  sum_a = 0
  for (i in 1:num_grid^2){
    sum_a = sum_a + grid_size^2*abs(cartesian_grid_vals[i,2] - obs_1[2])*observation_density(parameters, polar_grid_vals[i,], estimated_probabilities_hold_out[j,])
  }
  
  
  sum_aa = 0
  for (num_1 in 1:num_grid^2){
    for (num_2 in 1:num_grid^2){
      sum_aa = sum_aa + grid_size^4*abs(cartesian_grid_vals[num_1,2] - cartesian_grid_vals[num_2,2])*observation_density(parameters, polar_grid_vals[num_1,], estimated_probabilities_hold_out[j,])*observation_density(parameters, polar_grid_vals[num_2,], estimated_probabilities_hold_out[j,])
    }
  }
  
  energy_scores_1_v[num_iter] = sum_a - sum_aa/2
  
  sum_a = 0
  for (i in 1:num_grid^2){
    sum_a = sum_a + grid_size^2*abs(cartesian_grid_vals[i,2] - obs_2[2])*observation_density(parameters, polar_grid_vals[i,], estimated_probabilities_hold_out_2[j,])
  }
  
  
  sum_aa = 0
  for (num_1 in 1:num_grid^2){
    for (num_2 in 1:num_grid^2){
      sum_aa = sum_aa + grid_size^4*abs(cartesian_grid_vals[num_1,2] - cartesian_grid_vals[num_2,2])*observation_density(parameters, polar_grid_vals[num_1,], estimated_probabilities_hold_out_2[j,])*observation_density(parameters, polar_grid_vals[num_2,], estimated_probabilities_hold_out_2[j,])
    }
  }
  
  energy_scores_2_v[num_iter] = sum_a - sum_aa/2
  
  
  
  print(num_iter)
}

mean(energy_scores_1_u)
mean(energy_scores_1_v)
mean(energy_scores_2_u)
mean(energy_scores_2_v)


sswc = function(phi, param){
  return((1+param[5]*sin(phi-param[3]))/(2*pi*(cosh(param[4])-sinh(param[4])*cos(phi-param[3]))))
}


for (i in 1:10){
  j = sample(1:(24*24), 1)
  print(j)
  values = apply(vals, MARGIN= 1, FUN = observation_density, param=parameters, hold_out_prob = estimated_probabilities_hold_out[j,])
  par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
  contour(x=X_cor, y = y_cor, z = matrix(values,nrow=1000), xlab = "x", ylab = expression(paste(phi)), labcex = 0.01, cex.lab = 1.7)
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


X_cor = seq(0,2,l=200)
y_cor = seq(-pi,pi,l=200)
vals = cbind(rep(X_cor,200), rep(y_cor,each=200))
# X_cor = seq(0,0.6,l=100)
# y_cor = seq(-pi,pi,l=100)
# vals = cbind(rep(X_cor,100), rep(y_cor,each=100))

bins_1 = rep(0,576)
bins_2 = rep(0,576)

for (j in 1:576){
  values = apply(vals, MARGIN= 1, FUN = observation_density, param=parameters, hold_out_prob = estimated_probabilities_hold_out[j,])
  values = matrix(values,nrow=200)
  x_val = which.min((X_cor-simulated_sample[j,1])^2)
  y_val = which.min((y_cor-simulated_sample[j,2])^2)
  bins_1[j] = sum(values[which(values>values[x_val,y_val])])/sum(values)
  
  values = apply(vals, MARGIN= 1, FUN = observation_density, param=parameters, hold_out_prob = estimated_probabilities_hold_out_2[j,])
  values = matrix(values,nrow=200)
  x_val = which.min((X_cor-simulated_sample_2[j,1])^2)
  y_val = which.min((y_cor-simulated_sample_2[j,2])^2)
  bins_2[j] = sum(values[which(values>values[x_val,y_val])])/sum(values)
  
  print(j)
}
hist(bins_1)
hist(bins_2)

par(mgp=c(1.8,0.7,0),mar=c(3.7,3.7,2,2)+0.1)
contour(x=X_cor, y = y_cor, z = matrix(values,nrow=length(X_cor)), xlab = "x", ylab = expression(paste(phi)), cex.lab = 1.7, nlevels = 6, drawlabels = F, method = "flattest")
points(x = simulated_sample[j,1], y = simulated_sample[j,2], pch=16, cex = 1.4)
