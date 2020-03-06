full_likelihood_RMSE = function(parameter, data_sample){
  rho_val = parameter[1]
  theta_val = matrix(parameter[2:(5*ncolor+1)], nrow=ncolor)
  potts_prob = apply(xi_A, 1, function(i) exp(rho_val*ifelse(i[1]==i[2],1,0)))
  potts_prob = potts_prob/sum(potts_prob)
  value = 0
  for (i in 1:A){
    temp = 0
    #value = value + log(potts_prob%*%apply(xi_A, 1, function(j) dabeley_reparam(theta_val[j[1],], as.vector(simulated_sample[A_list[[i]][1],]))*dabeley_reparam(theta_val[j[2],], as.vector(simulated_sample[A_list[[i]][2],]))))
    for (j in 1:ncolor^2){
      temp = temp + potts_prob[j]*dabeley_reparam(theta_val[xi_A[j,1],], as.vector(data_sample[A_list[[i]][1],]))*dabeley_reparam(theta_val[xi_A[j,2],], as.vector(data_sample[A_list[[i]][2],]))
    }
    value = value + log(temp)
  }
  return(-value)
}

n_replicates = 50

elapsed_time_1 = rep(NA, n_replicates)
elapsed_time_cl = rep(NA, n_replicates)

parameter_estimates_1 = matrix(NA, nrow = n_replicates, ncol = 16)
parameter_estimates_cl = matrix(NA, nrow = n_replicates, ncol = 16)

true_field = matrix(NA, nrow = n_replicates, ncol = n_grid^2)
estimated_field_1 = matrix(NA, nrow = n_replicates, ncol = n_grid^2)
estimated_field_cl = matrix(NA, nrow = n_replicates, ncol = n_grid^2)
start_RMSE = rep(NA, n_replicates)

for (rep_num in 1:n_replicates){
  
  rho = 0.8
  potts_param <- c(rep(0, ncolor), rho)
  
  discrepancy = 0.5
  # theta_start = matrix(NA, nrow=ncolor, ncol = 5)
  # theta_start[,1] = runif(ncolor,1,5) # 1-3
  # theta_start[,2] = runif(ncolor,0.5,3) # 0.5-1
  # theta_start[,3] = runif(ncolor,-0.5,0.5) # -0.5-0.5
  # theta_start[,4] = runif(ncolor,0.05,3) # 0.05-3
  # theta_start[,5] = runif(ncolor,-0.9,0.9) # -0.9-0.9
  # 
  
  parameters = rbind(c(2,1,0,0,1), c(2,1,0,0,-1), c(2,0.6,0,1.5,0))
  #parameters = rbind(c(3,1,0,0.21,0.8), c(5,5,0,0.21,0), c(1,0.8,0,1.7,-0.8))
  parameters_test_reparam = c(rho, as.vector(parameters)) + runif(16,c(-0.8,-1.5,-1.5,-1.5,-0.6,-0.6,-0.2, -0.7,-0.7,-0.7, 0.01, 0.01, -1.2, -0.99,0.01,-0.5), c(0.2,1.5,1.5,1.5,1,1,1.4, 0.7, 0.7, 0.7, 1.5, 1.5, 1.5, 0.01,0.99,0.5))
  
  test1 = sqrt(mean((parameters_test_reparam-c(rho, as.vector(parameters)))^2))
  test2 = sqrt(mean((parameters_test_reparam-c(rho, as.vector(parameters[c(1,3,2),])))^2))
  test3 = sqrt(mean((parameters_test_reparam-c(rho, as.vector(parameters[c(2,3,1),])))^2))
  test4 = sqrt(mean((parameters_test_reparam-c(rho, as.vector(parameters[c(2,1,3),])))^2))
  test5 = sqrt(mean((parameters_test_reparam-c(rho, as.vector(parameters[c(3,2,1),])))^2))
  test6 = sqrt(mean((parameters_test_reparam-c(rho, as.vector(parameters[c(3,1,2),])))^2))
  start_RMSE[rep_num] = min(test1,test2,test3,test4,test5,test6)
  
  
  
  
  #parameters_test_reparam = c(rho, as.vector(theta_start))
  parameters_test_reparam[c(2,3,4,5,6,7,11,12,13)] = log(parameters_test_reparam[c(2,3,4,5,6,7,11,12,13)])
  parameters_test_reparam[c(8,9,10)] = atan(parameters_test_reparam[c(8,9,10)]/2)
  parameters_test_reparam[c(14,15,16)] = atanh(parameters_test_reparam[c(14,15,16)])
  init_param = parameters_test_reparam
  #init_param = c(runif(1,0,log(1+sqrt(ncolor_test))), runif(ncolor_test*5, parameters_test_reparam[2:(ncolor_test*5+1)]-discrepancy, parameters_test_reparam[2:(ncolor_test*5+1)]+discrepancy))
  
  
  
  
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
    theta_1 = rwrappedcauchy(1, mu = circular(parameters[k_i,3]), rho = tanh(parameters[k_i,4]/2))
    theta_1 = ifelse(as.numeric(theta_1)>pi, as.numeric(theta_1)-2*pi, as.numeric(theta_1))
    
    u = runif(1)
    simulated_sample$theta[i] = ifelse(u<(1+parameters[k_i,5]*sin(theta_1-parameters[k_i,3]))/2, theta_1, -theta_1)
    simulated_sample$x[i] = rweibull(1, scale = 1/(parameters[k_i,2]*(1-tanh(parameters[k_i,4])*cos(simulated_sample$theta[i]-parameters[k_i,3]))^(1/parameters[k_i,1])), shape = parameters[k_i,1])
  }
  simulated_sample = as.matrix(simulated_sample)
  
  print(full_likelihood_RMSE(parameters_test_reparam, simulated_sample))
  # Composite likelihood
  ttime = Sys.time()
  opt_test = tryCatch(
    optim(init_param, method = "BFGS", fn = full_likelihood_RMSE, control = list(trace = 6, REPORT = 1, reltol =1e-5), data_sample = simulated_sample),
    error = function(e){ 
      T
    }, finally = {}
  )
  if(class(opt_test)!="list"){next}
  elapsed_time_cl[rep_num] = Sys.time() - ttime
  
  estimated_param = rep(opt_test$par[1],1+5*ncolor_test)
  estimated_param[c(2,3,4,5,6,7,11,12,13)] = exp(opt_test$par[c(2,3,4,5,6,7,11,12,13)])
  estimated_param[c(8,9,10)] = 2*atan(opt_test$par[c(8,9,10)])
  estimated_param[c(14,15,16)] = tanh(opt_test$par[c(14,15,16)])
  parameter_estimates_cl[rep_num, ] = estimated_param
  
  theta_est = matrix(estimated_param[2:16], ncol = 5)
  rho_est = estimated_param[1]
  potts_prob = apply(xi_A, 1, function(i) exp(rho_est*ifelse(i[1]==i[2],1,0)))
  potts_prob = potts_prob/sum(potts_prob)
  xi_probs_i = matrix(0, nrow = n_grid^2, ncol = ncolor)
  for (i in 1:A){
    for (j in 1:ncolor^2){
      xi_probs[i,j] = potts_prob[j]*dabeley(theta_est[xi_A[j,1],], as.vector(simulated_sample[A_list[[i]][1],]))*dabeley(theta_est[xi_A[j,2],], as.vector(simulated_sample[A_list[[i]][2],]))
    }
    xi_probs[i,] = xi_probs[i,]/sum(xi_probs[i,])
  }
  for (i in 1:A){
    for (j in 1:ncolor^2){
      xi_probs_i[A_list[[i]][1],xi_A[j,1]] = xi_probs_i[A_list[[i]][1],xi_A[j,1]] + xi_probs[i,j]
      xi_probs_i[A_list[[i]][2],xi_A[j,2]] = xi_probs_i[A_list[[i]][2],xi_A[j,2]] + xi_probs[i,j]
    }
  }
  estimated_field_cl[rep_num, ] = apply(xi_probs_i,1,which.max)
  
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
  
  print(neg_likelihood_exact(parameters_test_reparam, n_rows = n_rows, data_sample = simulated_sample, n_cols = n_cols))
  ttime = Sys.time()
  optimal = tryCatch(
    optim(init_param, neg_likelihood_exact, method = "BFGS", control = list(trace=6, REPORT = 1, reltol = 1e-5), n_rows = n_rows, data_sample = simulated_sample, n_cols = n_cols),
    error = function(e){ 
      T
    }, finally = {}
  )
  if(class(optimal)!="list"){next}
  elapsed_time_1[rep_num] = Sys.time() - ttime
  
  
  
  estimated_param = rep(optimal$par[1],1+5*ncolor_test)
  estimated_param[c(2,3,4,5,6,7,11,12,13)] = exp(optimal$par[c(2,3,4,5,6,7,11,12,13)])
  estimated_param[c(8,9,10)] = 2*atan(optimal$par[c(8,9,10)])
  estimated_param[c(14,15,16)] = tanh(optimal$par[c(14,15,16)])
  parameter_estimates_1[rep_num, ] = estimated_param
  
  estimated_probabilities = find_back_probs(optimal$par, n_rows, simulated_sample, n_cols)
  estimated_field_1[rep_num, ] = apply(estimated_probabilities, 1, which.max)
  
  print(rep_num)
}

elapsed_time_1[which(elapsed_time_1>10)] = elapsed_time_1[which(elapsed_time_1>10)]/60
elapsed_time_cl[which(elapsed_time_cl>20)] = elapsed_time_cl[which(elapsed_time_cl>20)]/60
mean(elapsed_time_1, na.rm = T)
mean(elapsed_time_cl, na.rm = T)

1-(length(which(estimated_field_1 != true_field))/(n_replicates))/(24*24)
1-(length(which(estimated_field_cl != true_field))/(n_replicates))/(24*24)

true_parameters = c(rho, as.vector(parameters))
mean(apply(parameter_estimates_1, 1, function(i) sqrt(mean((i-true_parameters)^2))), na.rm = T)
mean(apply(parameter_estimates_cl, 1, function(i) sqrt(mean((i-true_parameters)^2))), na.rm = T)

diff_in_est_1 = rep(0,n_replicates)
diff_in_est_cl = rep(0,n_replicates)

for (j in 1:n_replicates){
  #print("1:")
  #print(matrix(parameter_estimates_1[j,2:16], nrow=3))
  #print("cl:")
  #print(matrix(parameter_estimates_cl[j,2:16], nrow=3))
  #print("true:")
  #print(matrix(true_parameters[2:16], nrow=3))
  #print(j)
  #print(sqrt(mean((parameter_estimates_1[j,]-c(rho, as.vector(parameters)))^2)))
  test1 = sqrt(mean((parameter_estimates_1[j,]-c(rho, as.vector(parameters)))^2))
  test2 = sqrt(mean((parameter_estimates_1[j,]-c(rho, as.vector(parameters[c(1,3,2),])))^2))
  test3 = sqrt(mean((parameter_estimates_1[j,]-c(rho, as.vector(parameters[c(2,3,1),])))^2))
  test4 = sqrt(mean((parameter_estimates_1[j,]-c(rho, as.vector(parameters[c(2,1,3),])))^2))
  test5 = sqrt(mean((parameter_estimates_1[j,]-c(rho, as.vector(parameters[c(3,2,1),])))^2))
  test6 = sqrt(mean((parameter_estimates_1[j,]-c(rho, as.vector(parameters[c(3,1,2),])))^2))
  diff_in_est_1[j] = min(test1,test2,test3,test4,test5,test6)
  
  test1 = sqrt(mean((parameter_estimates_cl[j,]-c(rho, as.vector(parameters)))^2))
  test2 = sqrt(mean((parameter_estimates_cl[j,]-c(rho, as.vector(parameters[c(1,3,2),])))^2))
  test3 = sqrt(mean((parameter_estimates_cl[j,]-c(rho, as.vector(parameters[c(2,3,1),])))^2))
  test4 = sqrt(mean((parameter_estimates_cl[j,]-c(rho, as.vector(parameters[c(2,1,3),])))^2))
  test5 = sqrt(mean((parameter_estimates_cl[j,]-c(rho, as.vector(parameters[c(3,2,1),])))^2))
  test6 = sqrt(mean((parameter_estimates_cl[j,]-c(rho, as.vector(parameters[c(3,1,2),])))^2))
  diff_in_est_cl[j] = min(test1,test2,test3,test4,test5,test6)
  # 
  # test <- matrix(parameter_estimates_cl[j,2:16], nrow=3)
  #   
  # values = apply(vals, MARGIN= 1, FUN = dabeley, param=test[1,])
  # image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
  # values = apply(vals, MARGIN= 1, FUN = dabeley, param=test[2,])
  # image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
  # values = apply(vals, MARGIN= 1, FUN = dabeley, param=test[3,])
  # image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
}
