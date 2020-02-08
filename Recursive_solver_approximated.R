n_grid_test = 30
ncolor_test = 3
rho_test = 0.5
potts_param_test <- c(rep(0, ncolor_test), rho_test)
x_potts_test <- matrix(1, nrow = n_grid_test, ncol = n_grid_test)
foo_test <- packPotts(x_potts_test, ncolor_test)
out_test <- potts(foo_test, potts_param_test, nbatch = 10)
image.plot(unpackPotts(out_test$final))

spat_pros_test = as.vector(unpackPotts(out_test$final))
parameters_test = rbind(c(3,1,0,0.21,0.8), c(5,5,0,0.21,0), c(1,0.8,0,1.7,-0.8))
values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters_test[1,])
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters_test[2,])
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters_test[3,])
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))

simulated_sample_test = data.frame(x = rep(NA,n_grid_test*n_grid_test), theta = rep(NA,n_grid_test*n_grid_test))
for(i in 1:(n_grid_test*n_grid_test)){
  k_i = spat_pros_test[i]
  theta_1_test = rwrappedcauchy(1, mu = circular(parameters_test[k_i,3]), rho = tanh(parameters_test[k_i,4]/2))
  theta_1_test = ifelse(as.numeric(theta_1_test)>pi, as.numeric(theta_1_test)-2*pi, as.numeric(theta_1_test))
  
  u = runif(1)
  simulated_sample_test$theta[i] = ifelse(u<(1+parameters_test[k_i,5]*sin(theta_1_test-parameters_test[k_i,3]))/2, theta_1_test, -theta_1_test)
  simulated_sample_test$x[i] = rweibull(1, scale = 1/(parameters_test[k_i,2]*(1-tanh(parameters_test[k_i,4])*cos(simulated_sample_test$theta[i]-parameters_test[k_i,3]))^(1/parameters_test[k_i,1])), shape = parameters_test[k_i,1])
}
ggplot(simulated_sample_test) + geom_point(aes(x=x, y = theta, col = as.factor(spat_pros_test))) + theme_bw() + theme(legend.position = "none")
simulated_sample_test = as.matrix(simulated_sample_test)

xi_A_2 = matrix(NA,nrow = ncolor_test^2,ncol = 2)
it = 1
for (row in 1:ncolor_test){
  for (column in 1:ncolor_test){
    xi_A_2[it, ] = c(row,column)
    it = it + 1
  }
}

xi_A_3 = matrix(NA,nrow = ncolor_test^3,ncol = 3)
it = 1
for (row in 1:ncolor_test){
  for (column in 1:ncolor_test){
    for (depth in 1:ncolor_test){
      xi_A_3[it, ] = c(row,column, depth)
      it = it + 1}
  }
}


likelihood <- function(parameters_input){
  rho_val = parameters_input[1]
  theta_val = matrix(parameters_input[2:(5*ncolor_test+1)], nrow=ncolor_test)
  
  # First initiate normality constants and forward probabilities
  norm_const = 0
  for (j in 1:ncolor_test){
    norm_const = norm_const + (1/ncolor_test)*dabeley_reparam(theta_val[j,], as.vector(simulated_sample_test[1,]))
  }
  norm_const = 1/norm_const
  norm_const = c(norm_const, rep(NA, n_grid_test^2-1))
  #norm_const = c(1/(dnorm(y_test[1],mu[1],1)*0.5+dnorm(y_test[1],mu[2],1)*0.5), rep(NA,8))
  forward_prob = matrix(0,nrow=n_grid_test^2,ncol = ncolor_test)
  for (j in 1:ncolor_test){
    forward_prob[1,j] = norm_const[1]*(1/ncolor_test)*dabeley_reparam(theta_val[j,], as.vector(simulated_sample_test[1,]))
  }
  #forward_prob[1,1] = norm_const[1]*0.5*dnorm(y_test[1],mu[1],1)
  #forward_prob[1,2] = norm_const[1]*0.5*dnorm(y_test[1],mu[2],1)
  transition = matrix(0, nrow = n_grid_test^2, ncol = ncolor_test^3)
  potts_prob_2 = apply(xi_A_2, 1, function(i) exp(rho_val*ifelse(i[1]==i[2],1,0)))
  potts_prob_2 = potts_prob_2/sum(potts_prob_2)
  potts_prob_3 = apply(xi_A_3, 1, function(i) exp(rho_val*(length(which(i[1:2]==i[3])))))
  potts_prob_3 = potts_prob_3/sum(potts_prob_3)
  
  # Iterate to find constants and probabilities
  for (t in 2:n_grid_test){
    #norm_const[t] = 1/(dnorm(y_test[t],0,tau)*p*forward_prob[t-1,1] + (1-p)*dnorm(y_test[t],0,tau)*forward_prob[t-1,2] + dnorm(y_test[t],1,tau)*(1-p)*forward_prob[t-1,1] + dnorm(y_test[t],1,tau)*p*forward_prob[t-1,2])
    for (j in 1:(ncolor_test^2)){
      transition[t,j] = potts_prob_2[j]*forward_prob[t-1,xi_A_2[j,1]]*dabeley_reparam(theta_val[xi_A_2[j,2],], as.vector(simulated_sample_test[t,]))}
    #transition[t,1] = potts_prob_2[1]*dnorm(y_test[t],mu[1],1)*forward_prob[t-1,1]
    #transition[t,2] = potts_prob_2[2]*dnorm(y_test[t],mu[2],1)*forward_prob[t-1,1]
    #transition[t,3] = potts_prob_2[3]*dnorm(y_test[t],mu[1],1)*forward_prob[t-1,2]
    #transition[t,4] = potts_prob_2[4]*dnorm(y_test[t],mu[2],1)*forward_prob[t-1,2]
    norm_const[t] = 1/sum(transition[t,])
    transition[t,] = transition[t,]*norm_const[t]
    for (j in 1:ncolor_test){
      forward_prob[t,j] = sum(transition[t,which(xi_A_2[,2]==j)])
    }
    #forward_prob[t,1] = transition[t,1] + transition[t,3]
    #forward_prob[t,2] = transition[t,2] + transition[t,4]
  }
  for (t in (n_grid_test+1):(n_grid_test^2)){
    if ((t-1)%%n_grid_test==0){
      for (j in 1:(ncolor_test^2)){
        transition[t,j] = potts_prob_2[j]*dabeley_reparam(theta_val[xi_A_2[j,2],], as.vector(simulated_sample_test[t,]))*forward_prob[t-n_grid_test,xi_A_2[j,1]]}
      #transition[t,1] = potts_prob_2[1]*dnorm(y_test[t],mu[1],1)*forward_prob[t-1,1]
      #transition[t,2] = potts_prob_2[2]*dnorm(y_test[t],mu[2],1)*forward_prob[t-1,1]
      #transition[t,3] = potts_prob_2[3]*dnorm(y_test[t],mu[1],1)*forward_prob[t-1,2]
      #transition[t,4] = potts_prob_2[4]*dnorm(y_test[t],mu[2],1)*forward_prob[t-1,2]
      norm_const[t] = 1/sum(transition[t,])
      transition[t,] = transition[t,]*norm_const[t]
      for (j in 1:ncolor_test){
        forward_prob[t,j] = sum(transition[t,which(xi_A_2[,2]==j)])
      }
      #forward_prob[t,1] = transition[t,1] + transition[t,3]
      #forward_prob[t,2] = transition[t,2] + transition[t,4]
    }else{
      for (j in 1:(ncolor_test^3)){
        transition[t,j] = potts_prob_3[j]*dabeley_reparam(theta_val[xi_A_3[j,3],], as.vector(simulated_sample_test[t,]))*forward_prob[t-1,xi_A_3[j,2]]*forward_prob[t-n_grid_test, xi_A_3[j,1]]
      }
      #transition[t,1] = potts_prob_3[1]*dnorm(y_test[t],mu[1],1)*forward_prob[t-1,1]*forward_prob[t-n_grid_test, 1]
      #transition[t,2] = potts_prob_3[2]*dnorm(y_test[t],mu[2],1)*forward_prob[t-1,1]*forward_prob[t-n_grid_test, 1]
      #transition[t,3] = potts_prob_3[3]*dnorm(y_test[t],mu[1],1)*forward_prob[t-1,2]*forward_prob[t-n_grid_test, 1]
      #transition[t,4] = potts_prob_3[4]*dnorm(y_test[t],mu[2],1)*forward_prob[t-1,2]*forward_prob[t-n_grid_test, 1]
      #transition[t,5] = potts_prob_3[5]*dnorm(y_test[t],mu[1],1)*forward_prob[t-1,1]*forward_prob[t-n_grid_test, 2]
      #transition[t,6] = potts_prob_3[6]*dnorm(y_test[t],mu[2],1)*forward_prob[t-1,1]*forward_prob[t-n_grid_test, 2]
      #transition[t,7] = potts_prob_3[7]*dnorm(y_test[t],mu[1],1)*forward_prob[t-1,2]*forward_prob[t-n_grid_test, 2]
      #transition[t,8] = potts_prob_3[8]*dnorm(y_test[t],mu[2],1)*forward_prob[t-1,2]*forward_prob[t-n_grid_test, 2]
      norm_const[t] = 1/sum(transition[t,])
      transition[t,] = transition[t,]*norm_const[t]
      for (j in 1:ncolor_test){
        forward_prob[t,j] = sum(transition[t,which(xi_A_3[,3]==j)])
      }
      #forward_prob[t,1] = transition[t,1] + transition[t,3] + transition[t,5] + transition[t,7]
      #forward_prob[t,2] = transition[t,2] + transition[t,4] + transition[t,6] + transition[t,8]
    }
  }
  #print(norm_const)
  #print(forward_prob)
  #print(transition)
  return(-sum(log(norm_const)))
  #return(norm_const)
}

neg_likelihood <- function(parameters){
  return(-likelihood(parameters))
}


discrepancy = 0.7
parameters_test_reparam = c(rho_test, as.vector(parameters_test[1:ncolor_test,]))
parameters_test_reparam[c(2,3,4,5,6,7,11,12,13)] = log(parameters_test_reparam[c(2,3,4,5,6,7,11,12,13)])
parameters_test_reparam[c(8,9,10)] = atan(parameters_test_reparam[c(8,9,10)]/2)
parameters_test_reparam[c(14,15,16)] = atanh(parameters_test_reparam[c(14,15,16)])
init_param = parameters_test_reparam
init_param = c(runif(1,0,log(1+sqrt(ncolor_test))), runif(ncolor_test*5, parameters_test_reparam[2:(ncolor_test*5+1)]-discrepancy, parameters_test_reparam[2:(ncolor_test*5+1)]+discrepancy))
#init_param = c(0.8, as.vector(parameters))
#init_param = c(0.5, as.vector(parameters_test[1:2,]))
#init_param[c(2,3,4,5,6,7,11,12,13)] = log(init_param[c(2,3,4,5,6,7,11,12,13)])
#init_param[c(8,9,10)] = tan(init_param[c(8,9,10)]/2)
#init_param[c(14,15,16)] = atanh(init_param[c(14,15,16)])
#init_param[c(2,3,4,5,8,9)] = log(init_param[c(2,3,4,5,8,9)])
#init_param[c(6,7)] = tan(init_param[c(6,7)]/2)
#init_param[c(10,11)] = atanh(init_param[c(10,11)])
likelihood(init_param)

time_test = Sys.time()
optimal = optim(init_param, neg_likelihood, method = "BFGS", control = list(trace=6, REPORT = 1))
#likelihood(c(0.5,-5,5,15,25, 35))
Sys.time() - time_test
optimal$par
estimated_param = rep(optimal$par[1],1+5*ncolor_test)
estimated_param[c(2,3,4,5,6,7,11,12,13)] = exp(optimal$par[c(2,3,4,5,6,7,11,12,13)])
estimated_param[c(8,9,10)] = 2*atan(optimal$par[c(8,9,10)])
estimated_param[c(14,15,16)] = tanh(optimal$par[c(14,15,16)])
#estimated_param[c(2,3,4,5,8,9)] = exp(optimal$par[c(2,3,4,5,8,9)])
#estimated_param[c(6,7)] = 2*atan(optimal$par[c(6,7)])
#estimated_param[c(10,11)] = tanh(optimal$par[c(10,11)])

estimated_param[1]
estimated_param = matrix(estimated_param[2:16],nrow=3)
estimated_param
parameters_test[1:ncolor_test,]



values = apply(vals, MARGIN= 1, FUN = dabeley, param=estimated_param[1,])
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
values = apply(vals, MARGIN= 1, FUN = dabeley, param=estimated_param[2,])
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
values = apply(vals, MARGIN= 1, FUN = dabeley, param=estimated_param[3,])
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
