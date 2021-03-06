parameters_test_reparam = read.table("C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_single_3_htlp.csv")[,1]

parameters = parameters_test_reparam

parameters[c(2,3,4,5,6,7,11,12,13)] = exp(parameters[c(2,3,4,5,6,7,11,12,13)])
parameters[c(8,9,10)] = 2*atan(parameters[c(8,9,10)])
parameters[c(14,15,16)] = exp(parameters[c(14,15,16)])/(1+exp(parameters[c(14,15,16)]))

parameters = matrix(parameters[2:(5*ncolor_test+1)],nrow=ncolor_test)
parameters

init_param = parameters_test_reparam


n_replicates = 200

elapsed_time_1 = rep(NA, n_replicates)

parameter_estimates_1 = matrix(NA, nrow = n_replicates, ncol = ncolor_test*5+1)

true_field_1 = matrix(NA, nrow = n_replicates, ncol = n_grid^2)
true_field_2 = matrix(NA, nrow = n_replicates, ncol = n_grid^2)
true_field_3 = matrix(NA, nrow = n_replicates, ncol = n_grid^2)
true_field_4 = matrix(NA, nrow = n_replicates, ncol = n_grid^2)
estimated_field_1 = matrix(NA, nrow = n_replicates, ncol = n_grid^2)
estimated_field_2 = matrix(NA, nrow = n_replicates, ncol = n_grid^2)
estimated_field_3 = matrix(NA, nrow = n_replicates, ncol = n_grid^2)
estimated_field_4 = matrix(NA, nrow = n_replicates, ncol = n_grid^2)
reverse = as.vector(t(matrix(1:(n_cols^2), nrow=n_cols)))

for (rep_num in 85:n_replicates){
  ## Draw random field ##
  true_field_1[rep_num,1:n_cols] = simulate_backward_probs_1_htlp(parameters_test_reparam, n_rows, simulated_sample[1:n_cols,], n_cols)
  true_field_2[rep_num,1:n_cols] = simulate_backward_probs_1_htlp(parameters_test_reparam, n_rows, simulated_sample_2[1:n_cols,], n_cols)
  true_field_3[rep_num,1:n_cols] = simulate_backward_probs_1_htlp(parameters_test_reparam, n_rows, simulated_sample_3[1:n_cols,], n_cols)
  true_field_4[rep_num,1:n_cols] = simulate_backward_probs_1_htlp(parameters_test_reparam, n_rows, simulated_sample_4[1:n_cols,], n_cols)
  for (i in 1:n_grid){
    true_field_1[rep_num,reverse[((i-1)*n_grid+1):(n_grid*i)]] = simulate_backward_probs_1_seed_htlp(parameters_test_reparam, n_rows, simulated_sample[reverse[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i)],], n_cols, x_seed = true_field_1[rep_num,i])
    true_field_2[rep_num,reverse[((i-1)*n_grid+1):(n_grid*i)]] = simulate_backward_probs_1_seed_htlp(parameters_test_reparam, n_rows, simulated_sample_2[reverse[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i)],], n_cols, x_seed = true_field_2[rep_num,i])
    true_field_3[rep_num,reverse[((i-1)*n_grid+1):(n_grid*i)]] = simulate_backward_probs_1_seed_htlp(parameters_test_reparam, n_rows, simulated_sample_3[reverse[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i)],], n_cols, x_seed = true_field_3[rep_num,i])
    true_field_4[rep_num,reverse[((i-1)*n_grid+1):(n_grid*i)]] = simulate_backward_probs_1_seed_htlp(parameters_test_reparam, n_rows, simulated_sample_4[reverse[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i)],], n_cols, x_seed = true_field_4[rep_num,i])
  }

  ## Draw observations ##
  simulated_sample_drawn = data.frame(x = rep(NA,n_grid*n_grid), theta = rep(NA,n_grid*n_grid))
  for(i in 1:(n_grid*n_grid)){
    k_i = true_field_1[rep_num,i]
    theta_1 = rwrappedcauchy(1, mu = circular(parameters[k_i,3]), rho = parameters[k_i,5]/(1+sqrt(1-parameters[k_i,5]^2)))
    simulated_sample_drawn$theta[i] = ifelse(as.numeric(theta_1)>pi, as.numeric(theta_1)-2*pi, as.numeric(theta_1))
    
    u = runif(1)
    if (parameters[k_i,4]==0){
      simulated_sample_drawn$x[i] = parameters[k_i,2]*(-log(1-u)/(1-parameters[k_i,5]*cos(simulated_sample_drawn$theta[i] - parameters[k_i,3])))^parameters[k_i,1]
    }else{
      simulated_sample_drawn$x[i] = parameters[k_i,2]*(((1-u)^(-parameters[k_i,4]/parameters[k_i,1])-1)/(parameters[k_i,4]/parameters[k_i,1]*(1-parameters[k_i,5]*cos(simulated_sample_drawn$theta[i] - parameters[k_i,3]))))^parameters[k_i,1]
    }
  }
  #ggplot(simulated_sample_drawn) + geom_point(aes(x=x, y = theta, col = as.factor(true_field_1[rep_num,]))) + theme_classic(base_size=20) + theme(legend.position = "none") + xlab("X") + ylab(expression(paste(Phi)))
  simulated_sample_drawn = as.matrix(simulated_sample_drawn)
  
  simulated_sample_drawn_2 = data.frame(x = rep(NA,n_grid*n_grid), theta = rep(NA,n_grid*n_grid))
  for(i in 1:(n_grid*n_grid)){
    k_i = true_field_2[rep_num,i]
    theta_1 = rwrappedcauchy(1, mu = circular(parameters[k_i,3]), rho = parameters[k_i,5]/(1+sqrt(1-parameters[k_i,5]^2)))
    simulated_sample_drawn_2$theta[i] = ifelse(as.numeric(theta_1)>pi, as.numeric(theta_1)-2*pi, as.numeric(theta_1))
    
    u = runif(1)
    if (parameters[k_i,4]==0){
      simulated_sample_drawn_2$x[i] = parameters[k_i,2]*(-log(1-u)/(1-parameters[k_i,5]*cos(simulated_sample_drawn_2$theta[i] - parameters[k_i,3])))^parameters[k_i,1]
    }else{
      simulated_sample_drawn_2$x[i] = parameters[k_i,2]*(((1-u)^(-parameters[k_i,4]/parameters[k_i,1])-1)/(parameters[k_i,4]/parameters[k_i,1]*(1-parameters[k_i,5]*cos(simulated_sample_drawn_2$theta[i] - parameters[k_i,3]))))^parameters[k_i,1]
    }
  }
  #ggplot(simulated_sample_drawn_2) + geom_point(aes(x=x, y = theta, col = as.factor(true_field_2[rep_num,]))) + theme_classic(base_size=20) + theme(legend.position = "none") + xlab("X") + ylab(expression(paste(Phi)))
  simulated_sample_drawn_2 = as.matrix(simulated_sample_drawn_2)
  
  simulated_sample_drawn_3 = data.frame(x = rep(NA,n_grid*n_grid), theta = rep(NA,n_grid*n_grid))
  for(i in 1:(n_grid*n_grid)){
    k_i = true_field_3[rep_num,i]
    theta_1 = rwrappedcauchy(1, mu = circular(parameters[k_i,3]), rho = parameters[k_i,5]/(1+sqrt(1-parameters[k_i,5]^2)))
    simulated_sample_drawn_3$theta[i] = ifelse(as.numeric(theta_1)>pi, as.numeric(theta_1)-2*pi, as.numeric(theta_1))
    
    u = runif(1)
    if (parameters[k_i,4]==0){
      simulated_sample_drawn_3$x[i] = parameters[k_i,2]*(-log(1-u)/(1-parameters[k_i,5]*cos(simulated_sample_drawn_3$theta[i] - parameters[k_i,3])))^parameters[k_i,1]
    }else{
      simulated_sample_drawn_3$x[i] = parameters[k_i,2]*(((1-u)^(-parameters[k_i,4]/parameters[k_i,1])-1)/(parameters[k_i,4]/parameters[k_i,1]*(1-parameters[k_i,5]*cos(simulated_sample_drawn_3$theta[i] - parameters[k_i,3]))))^parameters[k_i,1]
    }
  }
  #ggplot(simulated_sample_drawn_3) + geom_point(aes(x=x, y = theta, col = as.factor(true_field_3[rep_num,]))) + theme_classic(base_size=20) + theme(legend.position = "none") + xlab("X") + ylab(expression(paste(Phi)))
  simulated_sample_drawn_3 = as.matrix(simulated_sample_drawn_3)
  
  simulated_sample_drawn_4 = data.frame(x = rep(NA,n_grid*n_grid), theta = rep(NA,n_grid*n_grid))
  for(i in 1:(n_grid*n_grid)){
    k_i = true_field_4[rep_num,i]
    theta_1 = rwrappedcauchy(1, mu = circular(parameters[k_i,3]), rho = parameters[k_i,5]/(1+sqrt(1-parameters[k_i,5]^2)))
    simulated_sample_drawn_4$theta[i] = ifelse(as.numeric(theta_1)>pi, as.numeric(theta_1)-2*pi, as.numeric(theta_1))
    
    u = runif(1)
    if (parameters[k_i,4]==0){
      simulated_sample_drawn_4$x[i] = parameters[k_i,2]*(-log(1-u)/(1-parameters[k_i,5]*cos(simulated_sample_drawn_4$theta[i] - parameters[k_i,3])))^parameters[k_i,1]
    }else{
      simulated_sample_drawn_4$x[i] = parameters[k_i,2]*(((1-u)^(-parameters[k_i,4]/parameters[k_i,1])-1)/(parameters[k_i,4]/parameters[k_i,1]*(1-parameters[k_i,5]*cos(simulated_sample_drawn_4$theta[i] - parameters[k_i,3]))))^parameters[k_i,1]
    }
  }
  #ggplot(simulated_sample_drawn_4) + geom_point(aes(x=x, y = theta, col = as.factor(true_field_4[rep_num,]))) + theme_classic(base_size=20) + theme(legend.position = "none") + xlab("X") + ylab(expression(paste(Phi)))
  simulated_sample_drawn_4 = as.matrix(simulated_sample_drawn_4)
  
  ttime = proc.time()
  optimal = tryCatch(
    optim(init_param, neg_likelihood_exact_real_htlp_single, method = "BFGS", control = list(trace=6, REPORT = 1, reltol = 1e-5), n_rows = n_rows, data_sample = simulated_sample_drawn, data_sample_2 = simulated_sample_drawn_2, data_sample_3 = simulated_sample_drawn_3, data_sample_4 = simulated_sample_drawn_4, n_cols = n_cols),
    error = function(e){ 
      T
    }, finally = {}
  )
  if(class(optimal)!="list"){next}
  elapsed_time_1[rep_num] = (proc.time() - ttime)[[3]]
  
  
  
  estimated_param = rep(optimal$par[1],1+5*ncolor_test)

  estimated_param[c(2,3,4,5,6,7,11,12,13)] = exp(optimal$par[c(2,3,4,5,6,7,11,12,13)])
  estimated_param[c(8,9,10)] = 2*atan(optimal$par[c(8,9,10)])
  estimated_param[c(14,15,16)] = exp(optimal$par[c(14,15,16)])/(1+exp(optimal$par[c(14,15,16)]))
  
  parameter_estimates_1[rep_num, ] = estimated_param
  
  
  estimated_probabilities_1 = find_back_probs_htlp(optimal$par, n_rows, simulated_sample_drawn, n_cols)
  estimated_probabilities_2 = find_back_probs_htlp(optimal$par, n_rows, simulated_sample_drawn_2, n_cols)
  estimated_probabilities_3 = find_back_probs_htlp(optimal$par, n_rows, simulated_sample_drawn_3, n_cols)
  estimated_probabilities_4 = find_back_probs_htlp(optimal$par, n_rows, simulated_sample_drawn_4, n_cols)
  
  estimated_field_1[rep_num, ] = apply(estimated_probabilities_1, 1, which.max)
  estimated_field_2[rep_num, ] = apply(estimated_probabilities_2, 1, which.max)
  estimated_field_3[rep_num, ] = apply(estimated_probabilities_3, 1, which.max)
  estimated_field_4[rep_num, ] = apply(estimated_probabilities_4, 1, which.max)
  write.table(parameter_estimates_1[1:rep_num,], "C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_single_quantile_htlp.csv")
  print(rep_num)
}


mean(elapsed_time_1/60, na.rm = T)

matrix(apply(parameter_estimates_1[,2:(ncolor_test*5+1)], 2, mean, na.rm=T),nrow=ncolor_test)
mean(parameter_estimates_1[,1], na.rm = T)

1-(length(which(estimated_field_1 != true_field_1))/n_replicates)/(24*24)
1-(length(which(estimated_field_2 != true_field_2))/n_replicates)/(24*24)
1-(length(which(estimated_field_3 != true_field_3))/n_replicates)/(24*24)
1-(length(which(estimated_field_4 != true_field_4))/n_replicates)/(24*24)

true_parameters = c(parameters_test_reparam[1], as.vector(parameters))
mean(apply(parameter_estimates_1, 1, function(i) sqrt(mean((i-true_parameters)^2))), na.rm = T)


df = data.frame(y = as.vector(parameter_estimates_1))
df$x=c(rep(NA,n_replicates),rep(rep(1:ncolor_test,each=n_replicates),5))
df$x = as.factor(df$x)
df$supp = c(rep("rho", n_replicates), rep(c("alpha", "beta", "mu", "kappa", "lambda"), each =n_replicates*ncolor_test))

#pdf(file="C:/Users/henri/Documents/GitHub/Master-Thesis/Images/Fall_quantile_weissvm.pdf")
ggplot(df, aes(x=supp, y=y, fill=x)) +  geom_boxplot() + labs(fill = "Class") + scale_x_discrete(labels=c(expression(paste(alpha)),expression(paste(beta)),expression(paste(mu)),expression(paste(kappa)), expression(paste(lambda)), expression(paste(rho))), limits = c("alpha", "beta", "mu", "kappa", "lambda", "rho")) + xlab("Parameter") + ylab("Value") + theme_classic(base_size=20)# + coord_cartesian(ylim=c(-2,6))
#dev.off()


parameter_estimates_1 = as.matrix(read.table("C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_single_quantile_htlp.csv"))
parameter_estimates_sorted = apply(parameter_estimates_1,2,sort)

parameter_estimates_sorted[5,1] + parameters_test_reparam[1] - mean(parameter_estimates_sorted[,1])
parameter_estimates_sorted[195,1] + parameters_test_reparam[1] - mean(parameter_estimates_sorted[,1])
matrix(parameter_estimates_sorted[5,2:(ncolor_test*5+1)],ncol=5) + parameters - matrix(apply(parameter_estimates_sorted[,2:16], 2, mean),ncol=5)
parameters
matrix(parameter_estimates_sorted[195,2:(ncolor_test*5+1)],ncol=5) + parameters - matrix(apply(parameter_estimates_sorted[,2:16], 2, mean),ncol=5)
