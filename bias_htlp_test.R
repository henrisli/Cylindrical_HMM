rho = 0.2
potts_param <- c(rep(0, ncolor), rho)


# Case 1: "Hard"
#parameters = rbind(c(0.25,1,-pi/2,0,0.7), c(0.25,1,pi/2,0,0.7), c(0.5,0.5,0,0.8,0.8))

# Case 2: "Easy"
parameters = rbind(c(0.25,1,0,0,0.6), c(0.25,3,0,0,0.6), c(0.25,0.5,0,0,0.2))

values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[1,])
values[which(values==Inf)]=0
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[2,])
values[which(values==Inf)]=0
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
values = apply(vals, MARGIN= 1, FUN = dhtlp, param=parameters[3,])
values[which(values==Inf)]=0
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))


parameters_test_reparam = c(rho, as.vector(parameters)*(0.99)+0.001)
parameters_test_reparam[c(2,3,4,5,6,7,11,12,13)] = log(parameters_test_reparam[c(2,3,4,5,6,7,11,12,13)])
parameters_test_reparam[c(8,9,10)] = atan(parameters_test_reparam[c(8,9,10)]/2)
parameters_test_reparam[c(14,15,16)] = log(parameters_test_reparam[c(14,15,16)]/(1-parameters_test_reparam[c(14,15,16)]))
init_param = parameters_test_reparam


n_replicates = 200

elapsed_time_1 = rep(NA, n_replicates)

parameter_estimates_1 = matrix(NA, nrow = n_replicates, ncol = 16)

true_field = matrix(NA, nrow = n_replicates, ncol = n_grid^2)
estimated_field_1 = matrix(NA, nrow = n_replicates, ncol = n_grid^2)

for (rep_num in 60:n_replicates){
  # Draw random field
  x_potts <- matrix(1, nrow = n_grid, ncol = n_grid)
  foo <- packPotts(x_potts, ncolor)
  out <- potts(foo, potts_param, nbatch = 10)
  spat_pros = as.vector(unpackPotts(out$final))
  true_field[rep_num,] = spat_pros
  
  # Draw observations
  simulated_sample = data.frame(x = rep(NA,n_grid*n_grid), theta = rep(NA,n_grid*n_grid))
  for(i in 1:(n_grid*n_grid)){
    k_i = spat_pros[i]
    theta_1 = rwrappedcauchy(1, mu = circular(parameters[k_i,3]), rho = parameters[k_i,5]/(1+sqrt(1-parameters[k_i,5]^2)))
    simulated_sample$theta[i] = ifelse(as.numeric(theta_1)>pi, as.numeric(theta_1)-2*pi, as.numeric(theta_1))
    
    u = runif(1)
    if (parameters[k_i,4]==0){
      simulated_sample$x[i] = parameters[k_i,2]*(-log(1-u)/(1-parameters[k_i,5]*cos(simulated_sample$theta[i] - parameters[k_i,3])))^parameters[k_i,1]
    }else{
      simulated_sample$x[i] = parameters[k_i,2]*(((1-u)^(-parameters[k_i,4]/parameters[k_i,1])-1)/(parameters[k_i,4]/parameters[k_i,1]*(1-parameters[k_i,5]*cos(simulated_sample$theta[i] - parameters[k_i,3]))))^parameters[k_i,1]
    }
  }
  #ggplot(simulated_sample) + geom_point(aes(x=x, y = theta, col = as.factor(spat_pros))) + theme_classic(base_size=20) + theme(legend.position = "none") + xlab("X") + ylab(expression(paste(Phi)))
  simulated_sample = as.matrix(simulated_sample)
  
  
  
  # n_rows = 1
  n_rows = 1
  
  xi_A_n = matrix(1:3, ncol = 1)
  for (i in 1:(n_rows)){
    xi_A_n = cbind(rbind(xi_A_n, xi_A_n, xi_A_n), rep(1:3, each = 3^i))
  }
  moving_window_summation = matrix(0, nrow = ncolor_test^(n_rows), ncol = ncolor_test)
  for (j in 1:(ncolor_test^(n_rows))){
    moving_window_summation[j,] = seq((j-1)%/%(ncolor_test^2)*ncolor_test^2*(ncolor_test-1)+(j), (j-1)%/%(ncolor_test^2)*ncolor_test^2*2+j+ncolor_test^2*(ncolor_test-1), l=ncolor_test)}
  moving_window_summation = unique(moving_window_summation)
  seq_n_2 = 1:(ncolor_test^2)
  seq_n_3 = 1:(ncolor_test^3)
  conf_n_2 = c()
  for (i in 1:ncolor_test){
    conf_n_2 = c(conf_n_2, rep(seq_n_2[((i-1)*ncolor_test + 1):(ncolor_test*i)], ncolor_test^(n_rows-1)))
  }
  conf_n_3 = c()
  for (i in 1:(ncolor_test^2)){
    conf_n_3 = c(conf_n_3, rep(seq_n_3[((i-1)*ncolor_test + 1):(ncolor_test*i)], ncolor_test^(n_rows-2)))
  }
  
  ttime = Sys.time()
  optimal = tryCatch(
    optim(init_param, neg_likelihood_exact_htlp, method = "BFGS", control = list(trace=6, REPORT = 1, reltol = 1e-5), n_rows = n_rows, data_sample = simulated_sample, n_cols = n_cols),
    error = function(e){ 
      T
    }, finally = {}
  )
  if(class(optimal)!="list"){next}
  elapsed_time_1[rep_num] = Sys.time() - ttime
  
  
  estimated_param = rep(optimal$par[1],1+5*ncolor_test)
  estimated_param[c(2,3,4,5,6,7,11,12,13)] = exp(optimal$par[c(2,3,4,5,6,7,11,12,13)])
  estimated_param[c(8,9,10)] = 2*atan(optimal$par[c(8,9,10)])
  estimated_param[c(14,15,16)] = exp(optimal$par[c(14,15,16)])/(1+exp(optimal$par[c(14,15,16)]))
  parameter_estimates_1[rep_num, ] = estimated_param
  
  estimated_probabilities = find_back_probs_htlp(optimal$par, n_rows, simulated_sample, n_cols)
  estimated_field_1[rep_num, ] = apply(estimated_probabilities, 1, which.max)
  write.table(parameter_estimates_1[1:rep_num,], "C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_case2_02_htlp.csv")
  print(rep_num)
}

elapsed_time_1[which(elapsed_time_1>10)] = elapsed_time_1[which(elapsed_time_1>10)]/60

mean(elapsed_time_1, na.rm = T)

matrix(apply(parameter_estimates_1[,2:16], 2, mean),nrow=3)
mean(parameter_estimates_1[,1])

1-(length(which(estimated_field_1 != true_field))/n_replicates)/(24*24)

true_parameters = c(rho, as.vector(parameters))
mean(apply(parameter_estimates_1, 1, function(i) sqrt(mean((i-true_parameters)^2))), na.rm = T)

df = data.frame(y = as.vector(parameter_estimates_1))
df$x=c(rep(NA,n_replicates),rep(rep(1:3,each=n_replicates),5))
df$x = as.factor(df$x)
df$supp = c(rep("rho", n_replicates), rep(c("alpha", "beta", "mu", "tau", "kappa"), each =n_replicates*3))

#pdf(file="C:/Users/henri/Documents/GitHub/Master-Thesis/Images/Case1_02_bias_htlp.pdf")
ggplot(df, aes(x=supp, y=y, fill=x)) +  geom_boxplot() + labs(fill = "Class") + scale_x_discrete(labels=c(expression(paste(alpha)),expression(paste(beta)),expression(paste(mu)),expression(paste(tau)), expression(paste(kappa)), expression(paste(rho))), limits = c("alpha", "beta", "mu", "tau", "kappa", "rho")) + xlab("Parameter") + ylab("Value") + theme_classic(base_size=20) + coord_cartesian(ylim=c(-2,3.5))
#dev.off()

parameter_estimates_1_1_05 = as.matrix(read.table("C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_case1_05_htlp.csv"))
parameter_estimates_1_1_08 = as.matrix(read.table("C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_case1_08_htlp.csv"))
parameter_estimates_1_2_05 = as.matrix(read.table("C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_case2_05_htlp.csv"))
parameter_estimates_1_2_08 = as.matrix(read.table("C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_case2_08_htlp.csv"))
parameter_estimates_1_2_02 = as.matrix(read.table("C://Users//henri//Documents//GitHub//Master-Thesis//Data//parameter_estimates_case2_02_htlp.csv"))
df_test = data.frame(y105 = apply(parameter_estimates_1_1_05,2,var))
df_test$y205 = apply(parameter_estimates_1_2_05,2,var)
df_test$y108 = apply(parameter_estimates_1_1_08,2,var)
df_test$y208 = apply(parameter_estimates_1_2_08,2,var)
df_test$y202 = apply(parameter_estimates_1_2_02,2,var)
ggplot(df_test, aes(x=1:16)) + geom_point(aes(y=y202, col = "2, 0.2"), size = 3) + geom_point(aes(y=y105, col = "1, 0.5"), size = 3) + geom_point(aes(y=y108, col = "1, 0.8"), size = 3) + geom_point(aes(y=y205, col = "2, 0.5"), size = 3) + geom_point(aes(y=y208, col = "2, 0.8"), size = 3) + scale_x_discrete(labels=c(expression(paste(rho)), expression(paste(alpha)[1]), expression(paste(alpha)[2]), expression(paste(alpha)[3]), expression(paste(beta)[1]), expression(paste(beta)[2]), expression(paste(beta)[3]), expression(paste(mu)[1]), expression(paste(mu)[2]), expression(paste(mu)[3]), expression(paste(tau)[1]), expression(paste(tau)[2]), expression(paste(tau)[3]), expression(paste(kappa)[1]), expression(paste(kappa)[2]), expression(paste(kappa)[3])), limits = 1:16) + xlab("Parameter") + ylab("Variance") + theme_classic(base_size=20) + labs(col = expression(paste("Case, ",rho)))
