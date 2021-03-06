library(prodlim)
n_rows = 1
n_cols = 24
ncolor_test = 3
rho_test = 0.8
potts_param_test <- c(rep(0, ncolor_test), rho_test)
x_potts_test <- matrix(1, nrow = n_rows, ncol = n_cols)
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

simulated_sample_test = data.frame(x = rep(NA,n_cols*n_cols), theta = rep(NA,n_cols*n_cols))
for(i in 1:(n_cols*n_cols)){
  k_i = spat_pros_test[i]
  theta_1_test = rwrappedcauchy(1, mu = circular(parameters_test[k_i,3]), rho = tanh(parameters_test[k_i,4]/2))
  theta_1_test = ifelse(as.numeric(theta_1_test)>pi, as.numeric(theta_1_test)-2*pi, as.numeric(theta_1_test))
  
  u = runif(1)
  simulated_sample_test$theta[i] = ifelse(u<(1+parameters_test[k_i,5]*sin(theta_1_test-parameters_test[k_i,3]))/2, theta_1_test, -theta_1_test)
  simulated_sample_test$x[i] = rweibull(1, scale = 1/(parameters_test[k_i,2]*(1-tanh(parameters_test[k_i,4])*cos(simulated_sample_test$theta[i]-parameters_test[k_i,3]))^(1/parameters_test[k_i,1])), shape = parameters_test[k_i,1])
}
ggplot(simulated_sample_test) + geom_point(aes(x=x, y = theta, col = as.factor(spat_pros_test))) + theme_classic() + theme(legend.position = "none") + xlab("X") + ylab(expression(paste(Phi)))
simulated_sample_test = as.matrix(simulated_sample_test)

xi_A_2 = matrix(NA,nrow = ncolor_test^2,ncol = 2)
it = 1
for (column in 1:ncolor_test){
  for (row in 1:ncolor_test){
    xi_A_2[it, ] = c(row,column)
    it = it + 1
  }
}

xi_A_3 = matrix(NA,nrow = ncolor_test^3,ncol = 3)
it = 1
for (depth in 1:ncolor_test){
  for (column in 1:ncolor_test){
    for (row in 1:ncolor_test){
      xi_A_3[it, ] = c(row,column, depth)
      it = it + 1}
  }
}

xi_A_n = matrix(1:3, ncol = 1)
for (i in 1:(n_rows)){
  xi_A_n = cbind(rbind(xi_A_n, xi_A_n, xi_A_n), rep(1:3, each = 3^i))
}
moving_window_summation = matrix(0, nrow = ncolor_test^(n_rows), ncol = ncolor_test)
for (j in 1:(ncolor_test^(n_rows))){
  moving_window_summation[j,] = seq((j-1)%/%(ncolor_test^2)*ncolor_test^2*(ncolor_test-1)+(j), (j-1)%/%(ncolor_test^2)*ncolor_test^2*2+j+ncolor_test^2*(ncolor_test-1), l=ncolor_test)}
  #print(j)}
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

# Exact likelihood for m>1

likelihood_exact <- function(parameters_input, n_rows, data_sample, n_cols){
  rho_val = parameters_input[1]
  theta_val = matrix(parameters_input[2:(5*ncolor_test+1)], nrow=ncolor_test)
  
  # First initiate normality constants and forward probabilities
  norm_const = 0
  for (j in 1:ncolor_test){
    norm_const = norm_const + (1/ncolor_test)*dabeley_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  norm_const = 1/norm_const
  norm_const = c(norm_const, rep(NA, n_rows*n_cols-1))
  #norm_const = c(1/(dnorm(y_test[1],mu[1],1)*0.5+dnorm(y_test[1],mu[2],1)*0.5), rep(NA,8))
  forward_prob = matrix(0,nrow=n_rows*n_cols,ncol = ncolor_test)
  for (j in 1:ncolor_test){
    forward_prob[1,j] = norm_const[1]*(1/ncolor_test)*dabeley_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  #forward_prob[1,1] = norm_const[1]*0.5*dnorm(y_test[1],mu[1],1)
  #forward_prob[1,2] = norm_const[1]*0.5*dnorm(y_test[1],mu[2],1)
  transition = matrix(0, nrow = n_cols*n_rows, ncol = ncolor_test^(n_rows+1))
  moving_window = matrix(0, nrow = n_cols*n_rows, ncol = ncolor_test^(n_rows))
  moving_window[1, 1:ncolor_test] = forward_prob[1,]
  potts_prob_2 = apply(xi_A_2, 1, function(i) exp(rho_val*ifelse(i[1]==i[2],1,0)))
  potts_prob_2 = potts_prob_2/sum(potts_prob_2)
  potts_prob_3 = apply(xi_A_3, 1, function(i) exp(rho_val*(length(which(i[2:3]==i[1])))))
  potts_prob_3 = potts_prob_3/sum(potts_prob_3)
  
  # Iterate to find constants and probabilities
  for (t in 2:n_rows){
    #norm_const[t] = 1/(dnorm(y_test[t],0,tau)*p*forward_prob[t-1,1] + (1-p)*dnorm(y_test[t],0,tau)*forward_prob[t-1,2] + dnorm(y_test[t],1,tau)*(1-p)*forward_prob[t-1,1] + dnorm(y_test[t],1,tau)*p*forward_prob[t-1,2])
    for (j in 1:(ncolor_test^t)){
      transition[t,j] = potts_prob_2[(j-1)%/%(ncolor_test^(t-2))+1]*moving_window[t-1,(j-1)%%(ncolor_test^(t-1))+1]*dabeley_reparam(theta_val[xi_A_n[j,t],], as.vector(data_sample[t,]))}
    #transition[t,1] = potts_prob_2[1]*dnorm(y_test[t],mu[1],1)*forward_prob[t-1,1]
    #transition[t,2] = potts_prob_2[2]*dnorm(y_test[t],mu[2],1)*forward_prob[t-1,1]
    #transition[t,3] = potts_prob_2[3]*dnorm(y_test[t],mu[1],1)*forward_prob[t-1,2]
    #transition[t,4] = potts_prob_2[4]*dnorm(y_test[t],mu[2],1)*forward_prob[t-1,2]
    norm_const[t] = 1/sum(transition[t,])
    transition[t,] = transition[t,]*norm_const[t]
    #for (j in 1:(ncolor_test^(t-1))){
    #  moving_window[t,j] = sum(transition[t, which(apply(xi_A_n[1:(ncolor_test^t),1:(t-1)], 1, function(x) identical(x, xi_A_n[j,1:(t-1)])))])
    #}
    moving_window[t,1:ncolor_test^t] = transition[t,1:ncolor_test^t]
    for (j in 1:ncolor_test){
      forward_prob[t,j] = sum(transition[t,which(xi_A_n[,t]==j)])
    }
    #forward_prob[t,1] = transition[t,1] + transition[t,3]
    #forward_prob[t,2] = transition[t,2] + transition[t,4]
  }
  for (t in (n_rows+1):(n_rows*n_cols)){
    if ((t-1)%%n_rows==0){
      for (j in 1:(ncolor_test^(n_rows+1))){
        transition[t,j] = potts_prob_2[conf_n_2[j]]*dabeley_reparam(theta_val[xi_A_n[j,n_rows+1],], as.vector(data_sample[t,]))*moving_window[t-1,(j-1)%%(ncolor_test^(n_rows))+1]
      }
      # rekkefølge: t-n, t-(n-1), ... , t-1, t for telle oppover i potts_prob og utsummering for moving window
      norm_const[t] = 1/sum(transition[t,])
      transition[t,] = transition[t,]*norm_const[t]
      for (j in 1:(ncolor_test^n_rows)){
        moving_window[t,j] = sum(transition[t, (j-1)*ncolor_test+1:ncolor_test])
      }
      for (j in 1:ncolor_test){
        forward_prob[t,j] = sum(transition[t,which(xi_A_n[,n_rows+1]==j)])
      }
    }else{
      for (j in 1:(ncolor_test^(n_rows+1))){
        transition[t,j] = potts_prob_3[conf_n_3[j]]*dabeley_reparam(theta_val[xi_A_n[j,n_rows+1],], as.vector(data_sample[t,]))*moving_window[t-1,(j-1)%%(ncolor_test^(n_rows))+1]
      }
     
      norm_const[t] = 1/sum(transition[t,])
      transition[t,] = transition[t,]*norm_const[t]
      for (j in 1:(ncolor_test^n_rows)){
        moving_window[t,j] = sum(transition[t, (j-1)*ncolor_test+1:ncolor_test])
      }
      for (j in 1:ncolor_test){
        forward_prob[t,j] = sum(transition[t,which(xi_A_n[,n_rows+1]==j)])
      }
      #forward_prob[t,1] = transition[t,1] + transition[t,3] + transition[t,5] + transition[t,7]
      #forward_prob[t,2] = transition[t,2] + transition[t,4] + transition[t,6] + transition[t,8]
    }
    #print(t)
  }
  #print(norm_const)
  #print(forward_prob)
  #print(transition)
  return(-sum(log(norm_const)))
  #return(norm_const)
}

neg_likelihood_exact <- function(parameters_input, n_rows, data_sample, n_cols){
  ret_val = 0
  for (i in 1:(n_cols-n_rows+1)){
    if (n_rows==1){
    ret_val = ret_val -likelihood_exact_1(parameters_input, n_rows, data_sample[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),], n_cols)
    }else{ret_val = ret_val -likelihood_exact(parameters_input, n_rows, data_sample[1:(n_rows*n_cols) + (i-1)*n_cols,], n_cols)}
  }
  reverse = as.vector(t(matrix(1:(n_cols^2), nrow=n_cols)))
  data_sample = data_sample[reverse,]
  for (i in 1:(n_cols-n_rows+1)){
    if (n_rows==1){
      ret_val = ret_val -likelihood_exact_1(parameters_input, n_rows, data_sample[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),], n_cols)
    }else{ret_val = ret_val -likelihood_exact(parameters_input, n_rows, data_sample[1:(n_rows*n_cols) + (i-1)*n_cols,], n_cols)}
  }
  return(ret_val)
}

neg_likelihood_exact_hess <- function(parameters_input, n_rows, data_sample, n_cols){
  ret_val = 0
  for (i in 1:(n_cols-n_rows+1)){
    if (n_rows==1){
      ret_val = ret_val -likelihood_exact_1(parameters_input, n_rows, data_sample[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),], n_cols)
    }else{ret_val = ret_val -likelihood_exact(parameters_input, n_rows, data_sample[1:(n_rows*n_cols) + (i-1)*n_cols,], n_cols)}
  }
  return(ret_val)
}



# Exact likelihood with m=1

likelihood_exact_1 <- function(parameters_input, n_rows, data_sample, n_cols){
  rho_val = parameters_input[1]
  theta_val = matrix(parameters_input[2:(5*ncolor_test+1)], nrow=ncolor_test)
  
  # First initiate normality constants and forward probabilities
  norm_const = 0
  for (j in 1:ncolor_test){
    norm_const = norm_const + (1/ncolor_test)*dabeley_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  norm_const = 1/norm_const
  norm_const = c(norm_const, rep(NA, n_rows*n_cols-1))
  #norm_const = c(1/(dnorm(y_test[1],mu[1],1)*0.5+dnorm(y_test[1],mu[2],1)*0.5), rep(NA,8))
  forward_prob = matrix(0,nrow=n_rows*n_cols,ncol = ncolor_test)
  for (j in 1:ncolor_test){
    forward_prob[1,j] = norm_const[1]*(1/ncolor_test)*dabeley_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  #forward_prob[1,1] = norm_const[1]*0.5*dnorm(y_test[1],mu[1],1)
  #forward_prob[1,2] = norm_const[1]*0.5*dnorm(y_test[1],mu[2],1)
  transition = matrix(0, nrow = n_cols*n_rows, ncol = ncolor_test^(n_rows+1))
  moving_window = matrix(0, nrow = n_cols*n_rows, ncol = ncolor_test^(n_rows))
  moving_window[1, 1:ncolor_test] = forward_prob[1,]
  potts_prob_2 = apply(xi_A_2, 1, function(i) exp(rho_val*ifelse(i[1]==i[2],1,0)))
  potts_prob_2 = potts_prob_2/sum(potts_prob_2)
  
  # Iterate to find constants and probabilities
  for (t in 2:n_cols){
    #norm_const[t] = 1/(dnorm(y_test[t],0,tau)*p*forward_prob[t-1,1] + (1-p)*dnorm(y_test[t],0,tau)*forward_prob[t-1,2] + dnorm(y_test[t],1,tau)*(1-p)*forward_prob[t-1,1] + dnorm(y_test[t],1,tau)*p*forward_prob[t-1,2])
    for (j in 1:(ncolor_test^2)){
      transition[t,j] = potts_prob_2[j]*moving_window[t-1,(j-1)%%(ncolor_test^(2-1))+1]*dabeley_reparam(theta_val[xi_A_n[j,2],], as.vector(data_sample[t,]))}
    #transition[t,1] = potts_prob_2[1]*dnorm(y_test[t],mu[1],1)*forward_prob[t-1,1]
    #transition[t,2] = potts_prob_2[2]*dnorm(y_test[t],mu[2],1)*forward_prob[t-1,1]
    #transition[t,3] = potts_prob_2[3]*dnorm(y_test[t],mu[1],1)*forward_prob[t-1,2]
    #transition[t,4] = potts_prob_2[4]*dnorm(y_test[t],mu[2],1)*forward_prob[t-1,2]
    norm_const[t] = 1/sum(transition[t,])
    transition[t,] = transition[t,]*norm_const[t]
    #for (j in 1:(ncolor_test^(t-1))){
    #  moving_window[t,j] = sum(transition[t, which(apply(xi_A_n[1:(ncolor_test^t),1:(t-1)], 1, function(x) identical(x, xi_A_n[j,1:(t-1)])))])
    #}
    for (j in 1:(ncolor_test^n_rows)){
      moving_window[t,j] = sum(transition[t, (j-1)*ncolor_test+1:ncolor_test])
    }
    #moving_window[t,1:ncolor_test^2] = transition[t,1:ncolor_test^2]
    for (j in 1:ncolor_test){
      forward_prob[t,j] = sum(transition[t,which(xi_A_n[,2]==j)])
    }
    #forward_prob[t,1] = transition[t,1] + transition[t,3]
    #forward_prob[t,2] = transition[t,2] + transition[t,4]
  }
  #print(norm_const)
  #print(forward_prob)
  #print(transition)
  return(-sum(log(norm_const)))
  #return(norm_const)
}
