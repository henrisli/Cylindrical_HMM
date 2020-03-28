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


init_param = parameters_test_reparam


n_replicates = 200

elapsed_time_1 = rep(NA, n_replicates)

parameter_estimates_1 = matrix(NA, nrow = n_replicates, ncol = ncolor_test*5+1)

true_field_1 = matrix(NA, nrow = n_replicates, ncol = n_grid^2)
true_field_2 = matrix(NA, nrow = n_replicates, ncol = n_grid^2)
estimated_field_1 = matrix(NA, nrow = n_replicates, ncol = n_grid^2)
estimated_field_2 = matrix(NA, nrow = n_replicates, ncol = n_grid^2)
reverse = as.vector(t(matrix(1:(n_cols^2), nrow=n_cols)))

for (rep_num in 101:n_replicates){
  ## Draw random field ##
  true_field_1[rep_num,1:n_cols] = simulate_backward_probs_1(parameters_test_reparam, n_rows, simulated_sample[1:n_cols,], n_cols)
  true_field_2[rep_num,1:n_cols] = simulate_backward_probs_1(parameters_test_reparam, n_rows, simulated_sample_2[1:n_cols,], n_cols)
  for (i in 1:n_grid){
    true_field_1[rep_num,reverse[((i-1)*n_grid+1):(n_grid*i)]] = simulate_backward_probs_1_seed(parameters_test_reparam, n_rows, simulated_sample[reverse[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i)],], n_cols, x_seed = true_field_1[rep_num,i])
    true_field_2[rep_num,reverse[((i-1)*n_grid+1):(n_grid*i)]] = simulate_backward_probs_1_seed(parameters_test_reparam, n_rows, simulated_sample_2[reverse[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i)],], n_cols, x_seed = true_field_2[rep_num,i])
  }
  # for (i in 1:n_grid){
  #   true_field_1[rep_num,((i-1)*n_grid+1):(n_grid*i)] = simulate_backward_probs_1(parameters_test_reparam, n_rows, simulated_sample[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),], n_cols)
  #   true_field_2[rep_num,((i-1)*n_grid+1):(n_grid*i)] = simulate_backward_probs_1(parameters_test_reparam, n_rows, simulated_sample_2[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),], n_cols)
  # }
  
  ## Draw observations ##
  simulated_sample_drawn = data.frame(x = rep(NA,n_grid*n_grid), theta = rep(NA,n_grid*n_grid))
  for(i in 1:(n_grid*n_grid)){
    k_i = true_field_1[rep_num,i]
    theta_1 = rwrappedcauchy(1, mu = circular(parameters[k_i,3]), rho = tanh(parameters[k_i,4]/2))
    theta_1 = ifelse(as.numeric(theta_1)>pi, as.numeric(theta_1)-2*pi, as.numeric(theta_1))
    
    u = runif(1)
    simulated_sample_drawn$theta[i] = ifelse(u<(1+parameters[k_i,5]*sin(theta_1-parameters[k_i,3]))/2, theta_1, -theta_1)
    simulated_sample_drawn$x[i] = rweibull(1, scale = 1/(parameters[k_i,2]*(1-tanh(parameters[k_i,4])*cos(simulated_sample_drawn$theta[i]-parameters[k_i,3]))^(1/parameters[k_i,1])), shape = parameters[k_i,1])
  }
  #ggplot(simulated_sample_drawn) + geom_point(aes(x=x, y = theta, col = as.factor(true_field_1[rep_num,]))) + theme_classic(base_size=20) + theme(legend.position = "none") + xlab("X") + ylab(expression(paste(Phi)))
  simulated_sample_drawn = as.matrix(simulated_sample_drawn)
  
  simulated_sample_drawn_2 = data.frame(x = rep(NA,n_grid*n_grid), theta = rep(NA,n_grid*n_grid))
  for(i in 1:(n_grid*n_grid)){
    k_i = true_field_2[rep_num,i]
    theta_1 = rwrappedcauchy(1, mu = circular(parameters[k_i,3]), rho = tanh(parameters[k_i,4]/2))
    theta_1 = ifelse(as.numeric(theta_1)>pi, as.numeric(theta_1)-2*pi, as.numeric(theta_1))
    
    u = runif(1)
    simulated_sample_drawn_2$theta[i] = ifelse(u<(1+parameters[k_i,5]*sin(theta_1-parameters[k_i,3]))/2, theta_1, -theta_1)
    simulated_sample_drawn_2$x[i] = rweibull(1, scale = 1/(parameters[k_i,2]*(1-tanh(parameters[k_i,4])*cos(simulated_sample_drawn_2$theta[i]-parameters[k_i,3]))^(1/parameters[k_i,1])), shape = parameters[k_i,1])
  }
  #ggplot(simulated_sample_drawn_2) + geom_point(aes(x=x, y = theta, col = as.factor(true_field_2[rep_num,]))) + theme_classic(base_size=20) + theme(legend.position = "none") + xlab("X") + ylab(expression(paste(Phi)))
  simulated_sample_drawn_2 = as.matrix(simulated_sample_drawn_2)
  
  ttime = proc.time()
  optimal = tryCatch(
    optim(init_param, neg_likelihood_exact_real, method = "BFGS", control = list(trace=6, REPORT = 1, reltol = 1e-5), n_rows = n_rows, data_sample = simulated_sample_drawn, data_sample_2 = simulated_sample_drawn_2, n_cols = n_cols),
    error = function(e){ 
      T
    }, finally = {}
  )
  if(class(optimal)!="list"){next}
  elapsed_time_1[rep_num] = (proc.time() - ttime)[[3]]
  
  
  
  estimated_param = rep(optimal$par[1],1+5*ncolor_test)
  
  estimated_param[c(2,3,4,5,8,9)] = exp(optimal$par[c(2,3,4,5,8,9)])
  estimated_param[c(6,7)] = 2*atan(optimal$par[c(6,7)])
  estimated_param[c(10,11)] = tanh(optimal$par[c(10,11)])
  
  # estimated_param[c(2,3,4,5,6,7,11,12,13)] = exp(optimal$par[c(2,3,4,5,6,7,11,12,13)])
  # estimated_param[c(8,9,10)] = 2*atan(optimal$par[c(8,9,10)])
  # estimated_param[c(14,15,16)] = tanh(optimal$par[c(14,15,16)])
  
  parameter_estimates_1[rep_num, ] = estimated_param
  
  
  estimated_probabilities_1 = find_back_probs(optimal$par, n_rows, simulated_sample_drawn, n_cols)
  estimated_probabilities_2 = find_back_probs(optimal$par, n_rows, simulated_sample_drawn_2, n_cols)

  estimated_field_1[rep_num, ] = apply(estimated_probabilities_1, 1, which.max)
  estimated_field_2[rep_num, ] = apply(estimated_probabilities_2, 1, which.max)
  write.table(parameter_estimates_1[1:rep_num,], "C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_summer_quantile.csv")
  print(rep_num)
}


mean(elapsed_time_1/60, na.rm = T)

matrix(apply(parameter_estimates_1[,2:(ncolor_test*5+1)], 2, mean, na.rm=T),nrow=ncolor_test)
mean(parameter_estimates_1[,1], na.rm = T)

1-(length(which(estimated_field_1 != true_field_1))/n_replicates)/(24*24)
1-(length(which(estimated_field_2 != true_field_2))/n_replicates)/(24*24)

true_parameters = c(parameters_test_reparam[1], as.vector(parameters))
mean(apply(parameter_estimates_1, 1, function(i) sqrt(mean((i-true_parameters)^2))), na.rm = T)


df = data.frame(y = as.vector(parameter_estimates_1))
df$x=c(rep(NA,n_replicates),rep(rep(1:ncolor_test,each=n_replicates),5))
df$x = as.factor(df$x)
df$supp = c(rep("rho", n_replicates), rep(c("alpha", "beta", "mu", "kappa", "lambda"), each =n_replicates*ncolor_test))

#pdf(file="C:/Users/henri/Documents/GitHub/Master-Thesis/Images/Fall_quantile_weissvm.pdf")
ggplot(df, aes(x=supp, y=y, fill=x)) +  geom_boxplot() + labs(fill = "Class") + scale_x_discrete(labels=c(expression(paste(alpha)),expression(paste(beta)),expression(paste(mu)),expression(paste(kappa)), expression(paste(lambda)), expression(paste(rho))), limits = c("alpha", "beta", "mu", "kappa", "lambda", "rho")) + xlab("Parameter") + ylab("Value") + theme_classic(base_size=20) + coord_cartesian(ylim=c(-2,6))
#dev.off()


parameter_estimates_sorted = apply(parameter_estimates_1,2,sort)

matrix(parameter_estimates_sorted[5,2:(ncolor_test*5+1)],ncol=5)
parameters
matrix(parameter_estimates_sorted[195,2:(ncolor_test*5+1)],ncol=5)
