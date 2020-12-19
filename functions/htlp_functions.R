dhtlp_reparam <- function(param, x){
  alpha = exp(param[1])
  beta=exp(param[2])
  mu=2*atan(param[3])
  tau=exp(param[4])
  kappa=exp(param[5])/(1+exp(param[5]))
  if (tau==0){
    return(sqrt(1-kappa^2)/(2*pi*beta*alpha)*(x[1]/beta)^(1/alpha-1)*exp(-(x[1]/beta)^(1/alpha)*(1-kappa*cos(x[2]-mu))))
  }else{
    return(sqrt(1-kappa^2)/(2*pi*beta*alpha)*(x[1]/beta)^(1/alpha-1)*(1+tau/alpha*(x[1]/beta)^(1/alpha)*(1-kappa*cos(x[2]-mu)))^(-alpha/tau-1))
  }
}

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


likelihood_exact_1_htlp <- function(parameters_input, n_rows, data_sample, n_cols){
  rho_val = parameters_input[1]
  theta_val = matrix(parameters_input[2:(5*ncolor_test+1)], nrow=ncolor_test)
  
  # First initiate normality constants and forward probabilities
  norm_const = 0
  for (j in 1:ncolor_test){
    norm_const = norm_const + (1/ncolor_test)*dhtlp_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  norm_const = 1/norm_const
  norm_const = c(norm_const, rep(NA, n_rows*n_cols-1))
  #norm_const = c(1/(dnorm(y_test[1],mu[1],1)*0.5+dnorm(y_test[1],mu[2],1)*0.5), rep(NA,8))
  forward_prob = matrix(0,nrow=n_rows*n_cols,ncol = ncolor_test)
  for (j in 1:ncolor_test){
    forward_prob[1,j] = norm_const[1]*(1/ncolor_test)*dhtlp_reparam(theta_val[j,], as.vector(data_sample[1,]))
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
      transition[t,j] = potts_prob_2[j]*moving_window[t-1,(j-1)%%(ncolor_test^(2-1))+1]*dhtlp_reparam(theta_val[xi_A_n[j,2],], as.vector(data_sample[t,]))}
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

likelihood_exact_htlp <- function(parameters_input, n_rows, data_sample, n_cols){
  rho_val = parameters_input[1]
  theta_val = matrix(parameters_input[2:(5*ncolor_test+1)], nrow=ncolor_test)
  
  # First initiate normality constants and forward probabilities
  norm_const = 0
  for (j in 1:ncolor_test){
    norm_const = norm_const + (1/ncolor_test)*dhtlp_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  norm_const = 1/norm_const
  norm_const = c(norm_const, rep(NA, n_rows*n_cols-1))
  #norm_const = c(1/(dnorm(y_test[1],mu[1],1)*0.5+dnorm(y_test[1],mu[2],1)*0.5), rep(NA,8))
  forward_prob = matrix(0,nrow=n_rows*n_cols,ncol = ncolor_test)
  for (j in 1:ncolor_test){
    forward_prob[1,j] = norm_const[1]*(1/ncolor_test)*dhtlp_reparam(theta_val[j,], as.vector(data_sample[1,]))
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
      transition[t,j] = potts_prob_2[(j-1)%/%(ncolor_test^(t-2))+1]*moving_window[t-1,(j-1)%%(ncolor_test^(t-1))+1]*dhtlp_reparam(theta_val[xi_A_n[j,t],], as.vector(data_sample[t,]))}
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
        transition[t,j] = potts_prob_2[conf_n_2[j]]*dhtlp_reparam(theta_val[xi_A_n[j,n_rows+1],], as.vector(data_sample[t,]))*moving_window[t-1,(j-1)%%(ncolor_test^(n_rows))+1]
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
        transition[t,j] = potts_prob_3[conf_n_3[j]]*dhtlp_reparam(theta_val[xi_A_n[j,n_rows+1],], as.vector(data_sample[t,]))*moving_window[t-1,(j-1)%%(ncolor_test^(n_rows))+1]
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

neg_likelihood_exact_htlp <- function(parameters_input, n_rows, data_sample, n_cols){
  ret_val = 0
  for (i in 1:(n_cols-n_rows+1)){
    if (n_rows==1){
      ret_val = ret_val -likelihood_exact_1_htlp(parameters_input, n_rows, data_sample[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),], n_cols)
    }else{ret_val = ret_val -likelihood_exact_htlp(parameters_input, n_rows, data_sample[1:(n_rows*n_cols) + (i-1)*n_cols,], n_cols)}
  }
  reverse = as.vector(t(matrix(1:(n_cols^2), nrow=n_cols)))
  data_sample = data_sample[reverse,]
  for (i in 1:(n_cols-n_rows+1)){
    if (n_rows==1){
      ret_val = ret_val -likelihood_exact_1_htlp(parameters_input, n_rows, data_sample[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),], n_cols)
    }else{ret_val = ret_val -likelihood_exact_htlp(parameters_input, n_rows, data_sample[1:(n_rows*n_cols) + (i-1)*n_cols,], n_cols)}
  }
  return(ret_val)
}

neg_likelihood_exact_hess_htlp <- function(parameters_input, n_rows, data_sample, n_cols){
  ret_val = 0
  for (i in 1:(n_cols-n_rows+1)){
    if (n_rows==1){
      ret_val = ret_val -likelihood_exact_1_htlp(parameters_input, n_rows, data_sample[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),], n_cols)
    }else{ret_val = ret_val -likelihood_exact_htlp(parameters_input, n_rows, data_sample[1:(n_rows*n_cols) + (i-1)*n_cols,], n_cols)}
  }
  return(ret_val)
}

backward_probs_htlp <- function(parameters_input, n_rows, data_sample, n_cols){
  rho_val = parameters_input[1]
  theta_val = matrix(parameters_input[2:(5*ncolor_test+1)], nrow=ncolor_test)
  
  # First initiate normality constants and forward probabilities
  norm_const = 0
  for (j in 1:ncolor_test){
    norm_const = norm_const + (1/ncolor_test)*dhtlp_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  norm_const = 1/norm_const
  norm_const = c(norm_const, rep(NA, n_rows*n_cols-1))
  forward_prob = matrix(0,nrow=n_rows*n_cols,ncol = ncolor_test)
  for (j in 1:ncolor_test){
    forward_prob[1,j] = norm_const[1]*(1/ncolor_test)*dhtlp_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  transition = matrix(0, nrow = n_cols*n_rows, ncol = ncolor_test^(n_rows+1))
  moving_window = matrix(0, nrow = n_cols*n_rows, ncol = ncolor_test^(n_rows))
  moving_window[1, 1:ncolor_test] = forward_prob[1,]
  potts_prob_2 = apply(xi_A_2, 1, function(i) exp(rho_val*ifelse(i[1]==i[2],1,0)))
  potts_prob_2 = potts_prob_2/sum(potts_prob_2)
  potts_prob_3 = apply(xi_A_3, 1, function(i) exp(rho_val*(length(which(i[2:3]==i[1])))))
  potts_prob_3 = potts_prob_3/sum(potts_prob_3)
  
  # Iterate to find constants and probabilities
  for (t in 2:n_rows){
    for (j in 1:(ncolor_test^t)){
      transition[t,j] = potts_prob_2[(j-1)%/%(ncolor_test^(t-2))+1]*moving_window[t-1,(j-1)%%(ncolor_test^(t-1))+1]*dhtlp_reparam(theta_val[xi_A_n[j,t],], as.vector(data_sample[t,]))}
    norm_const[t] = 1/sum(transition[t,])
    transition[t,] = transition[t,]*norm_const[t]
    moving_window[t,1:ncolor_test^t] = transition[t,1:ncolor_test^t]
    for (j in 1:ncolor_test){
      forward_prob[t,j] = sum(transition[t,which(xi_A_n[,t]==j)])
    }
  }
  for (t in (n_rows+1):(n_rows*n_cols)){
    if ((t-1)%%n_rows==0){
      for (j in 1:(ncolor_test^(n_rows+1))){
        transition[t,j] = potts_prob_2[conf_n_2[j]]*dhtlp_reparam(theta_val[xi_A_n[j,n_rows+1],], as.vector(data_sample[t,]))*moving_window[t-1,(j-1)%%(ncolor_test^(n_rows))+1]
      }
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
        transition[t,j] = potts_prob_3[conf_n_3[j]]*dhtlp_reparam(theta_val[xi_A_n[j,n_rows+1],], as.vector(data_sample[t,]))*moving_window[t-1,(j-1)%%(ncolor_test^(n_rows))+1]
      }
      
      norm_const[t] = 1/sum(transition[t,])
      transition[t,] = transition[t,]*norm_const[t]
      for (j in 1:(ncolor_test^n_rows)){
        moving_window[t,j] = sum(transition[t, (j-1)*ncolor_test+1:ncolor_test])
      }
      for (j in 1:ncolor_test){
        forward_prob[t,j] = sum(transition[t,which(xi_A_n[,n_rows+1]==j)])
      }
    }
  }
  back_transition = matrix(NA, ncol = ncolor_test^2, nrow = n_rows*n_cols)
  forward_transition = matrix(NA, ncol = ncolor_test^2, nrow = n_rows*n_cols)
  backward_prob = matrix(NA, ncol = ncolor_test, nrow = n_rows*n_cols)
  backward_prob[n_rows*n_cols,] = forward_prob[n_rows*n_cols, ]
  for (t in (n_rows*n_cols):(n_rows+1)){
    for (j in 1:(ncolor_test^2)){
      forward_transition[t,j] = sum(transition[t, (j-1)*ncolor_test^(n_rows-1) + 1:(ncolor_test^(n_rows-1))])
    }
    for (j in 1:(ncolor_test^2)){
      if (forward_prob[t,(j-1)%/%ncolor_test+1]==0 & forward_transition[t,j]==0){back_transition[t,j] = backward_prob[t, (j-1)%/%ncolor_test+1]
      }else{back_transition[t,j] = forward_transition[t,j]*backward_prob[t, (j-1)%/%ncolor_test+1]/forward_prob[t, (j-1)%/%ncolor_test+1]}
    }
    for (j in 1:ncolor_test){
      backward_prob[t-1, j] = sum(back_transition[t, j+0:(ncolor_test-1)*ncolor_test])
    }
  }
  for (t in n_rows:2){
    for (j in 1:(ncolor_test^2)){
      forward_transition[t,j] = sum(transition[t, (j-1)*ncolor_test^(t-2) + 1:(ncolor_test^(t-2))])
    }
    for (j in 1:(ncolor_test^2)){
      if (forward_prob[t, (j-1)%/%ncolor_test+1]==0 & forward_transition[t,j]==0){back_transition[t,j] = backward_prob[t, (j-1)%/%ncolor_test+1]
      }else{back_transition[t,j] = forward_transition[t,j]*backward_prob[t, (j-1)%/%ncolor_test+1]/forward_prob[t, (j-1)%/%ncolor_test+1]}
    }
    for (j in 1:ncolor_test){
      backward_prob[t-1,j] = sum(back_transition[t, j+0:(ncolor_test-1)*ncolor_test])
    }
  }
  return(backward_prob)
}

backward_probs_1_htlp <- function(parameters_input, n_rows, data_sample, n_cols){
  rho_val = parameters_input[1]
  theta_val = matrix(parameters_input[2:(5*ncolor_test+1)], nrow=ncolor_test)
  
  # First initiate normality constants and forward probabilities
  norm_const = 0
  for (j in 1:ncolor_test){
    norm_const = norm_const + (1/ncolor_test)*dhtlp_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  norm_const = 1/norm_const
  norm_const = c(norm_const, rep(NA, n_rows*n_cols-1))
  forward_prob = matrix(0,nrow=n_rows*n_cols,ncol = ncolor_test)
  for (j in 1:ncolor_test){
    forward_prob[1,j] = norm_const[1]*(1/ncolor_test)*dhtlp_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  transition = matrix(0, nrow = n_cols*n_rows, ncol = ncolor_test^(n_rows+1))
  moving_window = matrix(0, nrow = n_cols*n_rows, ncol = ncolor_test^(n_rows))
  moving_window[1, 1:ncolor_test] = forward_prob[1,]
  potts_prob_2 = apply(xi_A_2, 1, function(i) exp(rho_val*ifelse(i[1]==i[2],1,0)))
  potts_prob_2 = potts_prob_2/sum(potts_prob_2)
  
  # Iterate to find constants and probabilities
  for (t in 2:n_cols){
    for (j in 1:(ncolor_test^2)){
      transition[t,j] = potts_prob_2[j]*moving_window[t-1,(j-1)%%(ncolor_test^(2-1))+1]*dhtlp_reparam(theta_val[xi_A_n[j,2],], as.vector(data_sample[t,]))}
    norm_const[t] = 1/sum(transition[t,])
    transition[t,] = transition[t,]*norm_const[t]
    for (j in 1:(ncolor_test^n_rows)){
      moving_window[t,j] = sum(transition[t, (j-1)*ncolor_test+1:ncolor_test])
    }
    
    for (j in 1:ncolor_test){
      forward_prob[t,j] = sum(transition[t,which(xi_A_n[,2]==j)])
    }
  }
  back_transition = matrix(NA, ncol = ncolor_test^2, nrow = n_rows*n_cols)
  backward_prob = matrix(NA, ncol = ncolor_test, nrow = n_rows*n_cols)
  backward_prob[n_rows*n_cols,] = forward_prob[n_rows*n_cols, ]
  for (t in n_cols:2){
    for (j in 1:(ncolor_test^2)){
      if (forward_prob[t, (j-1)%/%ncolor_test+1]==0 & transition[t,j]==0){back_transition[t,j] = backward_prob[t, (j-1)%/%ncolor_test+1]
      }else{back_transition[t,j] = transition[t,j]*backward_prob[t, (j-1)%/%ncolor_test+1]/forward_prob[t, (j-1)%/%ncolor_test+1]}
    }
    for (j in 1:ncolor_test){
      backward_prob[t-1,j] = sum(back_transition[t, j+0:(ncolor_test-1)*ncolor_test])
    }
  }
  return(backward_prob)
}

# Simulate from the forward-backward algorithm
simulate_backward_probs_1_htlp <- function(parameters_input, n_rows, data_sample, n_cols){
  rho_val = parameters_input[1]
  theta_val = matrix(parameters_input[2:(5*ncolor_test+1)], nrow=ncolor_test)
  
  # First initiate normality constants and forward probabilities
  norm_const = 0
  for (j in 1:ncolor_test){
    norm_const = norm_const + (1/ncolor_test)*dhtlp_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  norm_const = 1/norm_const
  norm_const = c(norm_const, rep(NA, n_rows*n_cols-1))
  forward_prob = matrix(0,nrow=n_rows*n_cols,ncol = ncolor_test)
  for (j in 1:ncolor_test){
    forward_prob[1,j] = norm_const[1]*(1/ncolor_test)*dhtlp_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  transition = matrix(0, nrow = n_cols*n_rows, ncol = ncolor_test^(n_rows+1))
  moving_window = matrix(0, nrow = n_cols*n_rows, ncol = ncolor_test^(n_rows))
  moving_window[1, 1:ncolor_test] = forward_prob[1,]
  potts_prob_2 = apply(xi_A_2, 1, function(i) exp(rho_val*ifelse(i[1]==i[2],1,0)))
  potts_prob_2 = potts_prob_2/sum(potts_prob_2)
  
  # Iterate to find constants and probabilities
  for (t in 2:n_cols){
    for (j in 1:(ncolor_test^2)){
      transition[t,j] = potts_prob_2[j]*moving_window[t-1,(j-1)%%(ncolor_test^(2-1))+1]*dhtlp_reparam(theta_val[xi_A_n[j,2],], as.vector(data_sample[t,]))}
    norm_const[t] = 1/sum(transition[t,])
    transition[t,] = transition[t,]*norm_const[t]
    for (j in 1:(ncolor_test^n_rows)){
      moving_window[t,j] = sum(transition[t, (j-1)*ncolor_test+1:ncolor_test])
    }
    
    for (j in 1:ncolor_test){
      forward_prob[t,j] = sum(transition[t,which(xi_A_n[,2]==j)])
    }
  }
  back_transition = matrix(NA, ncol = ncolor_test^2, nrow = n_rows*n_cols)
  backward_prob = matrix(NA, ncol = ncolor_test, nrow = n_rows*n_cols)
  backward_prob[n_rows*n_cols,] = forward_prob[n_rows*n_cols, ]
  propagation = matrix(NA, ncol = ncolor_test^2, nrow = n_rows*n_cols)
  for (t in n_cols:2){
    for (j in 1:(ncolor_test^2)){
      if (forward_prob[t, (j-1)%/%ncolor_test+1]==0 & transition[t,j]==0){back_transition[t,j] = backward_prob[t, (j-1)%/%ncolor_test+1]
      }else{back_transition[t,j] = transition[t,j]*backward_prob[t, (j-1)%/%ncolor_test+1]/forward_prob[t, (j-1)%/%ncolor_test+1]}
    }
    for (j in 1:ncolor_test){
      backward_prob[t-1,j] = sum(back_transition[t, j+0:(ncolor_test-1)*ncolor_test])
    }
    for (j in 1:(ncolor_test^2)){
      propagation[t,j] = back_transition[t,j]/backward_prob[t-1, (j-1)%%ncolor_test+1]
    }
  }
  x_sample = rep(NA, n_cols*n_rows)
  x_sample[1] = which(rmultinom(1,1,backward_prob[1,])==1)
  for (i in 2:n_cols){
    x_sample[i] = which(rmultinom(1,1,propagation[i,x_sample[i-1]+0:(ncolor_test-1)*ncolor_test])==1)
  }
  return(x_sample)
}

# Simulate from the forward-backward algorithm with seed
simulate_backward_probs_1_seed_htlp <- function(parameters_input, n_rows, data_sample, n_cols, x_seed){
  rho_val = parameters_input[1]
  theta_val = matrix(parameters_input[2:(5*ncolor_test+1)], nrow=ncolor_test)
  
  # First initiate normality constants and forward probabilities
  norm_const = 0
  for (j in 1:ncolor_test){
    norm_const = norm_const + (1/ncolor_test)*dhtlp_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  norm_const = 1/norm_const
  norm_const = c(norm_const, rep(NA, n_rows*n_cols-1))
  forward_prob = matrix(0,nrow=n_rows*n_cols,ncol = ncolor_test)
  for (j in 1:ncolor_test){
    forward_prob[1,j] = norm_const[1]*(1/ncolor_test)*dhtlp_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  transition = matrix(0, nrow = n_cols*n_rows, ncol = ncolor_test^(n_rows+1))
  moving_window = matrix(0, nrow = n_cols*n_rows, ncol = ncolor_test^(n_rows))
  moving_window[1, 1:ncolor_test] = forward_prob[1,]
  potts_prob_2 = apply(xi_A_2, 1, function(i) exp(rho_val*ifelse(i[1]==i[2],1,0)))
  potts_prob_2 = potts_prob_2/sum(potts_prob_2)
  
  # Iterate to find constants and probabilities
  for (t in 2:n_cols){
    for (j in 1:(ncolor_test^2)){
      transition[t,j] = potts_prob_2[j]*moving_window[t-1,(j-1)%%(ncolor_test^(2-1))+1]*dhtlp_reparam(theta_val[xi_A_n[j,2],], as.vector(data_sample[t,]))}
    norm_const[t] = 1/sum(transition[t,])
    transition[t,] = transition[t,]*norm_const[t]
    for (j in 1:(ncolor_test^n_rows)){
      moving_window[t,j] = sum(transition[t, (j-1)*ncolor_test+1:ncolor_test])
    }
    
    for (j in 1:ncolor_test){
      forward_prob[t,j] = sum(transition[t,which(xi_A_n[,2]==j)])
    }
  }
  back_transition = matrix(NA, ncol = ncolor_test^2, nrow = n_rows*n_cols)
  backward_prob = matrix(NA, ncol = ncolor_test, nrow = n_rows*n_cols)
  backward_prob[n_rows*n_cols,] = forward_prob[n_rows*n_cols, ]
  propagation = matrix(NA, ncol = ncolor_test^2, nrow = n_rows*n_cols)
  for (t in n_cols:2){
    for (j in 1:(ncolor_test^2)){
      if (forward_prob[t, (j-1)%/%ncolor_test+1]==0 & transition[t,j]==0){back_transition[t,j] = backward_prob[t, (j-1)%/%ncolor_test+1]
      }else{back_transition[t,j] = transition[t,j]*backward_prob[t, (j-1)%/%ncolor_test+1]/forward_prob[t, (j-1)%/%ncolor_test+1]}
    }
    for (j in 1:ncolor_test){
      backward_prob[t-1,j] = sum(back_transition[t, j+0:(ncolor_test-1)*ncolor_test])
    }
    for (j in 1:(ncolor_test^2)){
      propagation[t,j] = back_transition[t,j]/backward_prob[t-1, (j-1)%%ncolor_test+1]
    }
  }
  x_sample = rep(NA, n_cols*n_rows)
  x_sample[1] = x_seed
  for (i in 2:n_cols){
    x_sample[i] = which(rmultinom(1,1,propagation[i,x_sample[i-1]+0:(ncolor_test-1)*ncolor_test])==1)
  }
  return(x_sample)
}



# backward_probabilities = matrix(0, nrow = n_cols^2, ncol = ncolor_test)
# for (i in 1:(n_cols/n_rows)){
#   backward_probabilities[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),] = backward_probs_1_htlp(init_param, n_rows, simulated_sample_temp[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),], n_cols)
# }

find_back_probs_htlp <- function(parameters_input, n_rows, data_sample, n_cols){
  backward_probabilities = matrix(0, nrow = n_cols^2, ncol = ncolor_test)
  for (i in 1:(n_cols/n_rows)){
    if (n_rows==1){
      backward_probabilities[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),] = backward_probs_1_htlp(parameters_input, n_rows, data_sample[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),], n_cols)
    }else{backward_probabilities[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),] = backward_probs_htlp(parameters_input, n_rows, data_sample[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),], n_cols)}
  }
  reverse = as.vector(t(matrix(1:(n_cols^2), nrow=n_cols)))
  data_sample = data_sample[reverse,]
  for (i in 1:(n_cols/n_rows)){
    if (n_rows==1){
      backward_probabilities[reverse[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i)],] = backward_probabilities[reverse[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i)],] + backward_probs_1_htlp(parameters_input, n_rows, data_sample[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),], n_cols)
    }else{backward_probabilities[reverse[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i)],] = backward_probabilities[reverse[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i)],] + backward_probs_htlp(parameters_input, n_rows, data_sample[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),], n_cols)}
  }
  return(backward_probabilities/2)
}

full_likelihood_htlp = function(parameter){
  rho_val = parameter[1]
  theta_val = matrix(parameter[2:(5*ncolor+1)], nrow=ncolor)
  potts_prob = apply(xi_A, 1, function(i) exp(rho_val*ifelse(i[1]==i[2],1,0)))
  potts_prob = potts_prob/sum(potts_prob)
  value = 0
  for (i in 1:A){
    temp = 0
    #value = value + log(potts_prob%*%apply(xi_A, 1, function(j) dhtlp_reparam(theta_val[j[1],], as.vector(simulated_sample[A_list[[i]][1],]))*dhtlp_reparam(theta_val[j[2],], as.vector(simulated_sample[A_list[[i]][2],]))))
    for (j in 1:ncolor^2){
      temp = temp + potts_prob[j]*dhtlp_reparam(theta_val[xi_A[j,1],], as.vector(simulated_sample[A_list[[i]][1],]))*dhtlp_reparam(theta_val[xi_A[j,2],], as.vector(simulated_sample[A_list[[i]][2],]))
    }
    value = value + log(temp)
  }
  return(-value)
}

theta_function_htlp = function(xi_probs_i_est_and_sample, theta_val){
  theta_val = matrix(theta_val, nrow = ncolor)
  #value = 0
  #for (i in 1:A){
  #  j = A_list[[i]][1]
  #  k = A_list[[i]][2]
  #  value = value + sum(xi_probs_i_est_and_sample[j,1:ncolor]*log(apply(theta_val, 1, dhtlp_reparam, x=simulated_sample[j,])))
  #  value = value + sum(xi_probs_i_est_and_sample[k,1:ncolor]*log(apply(theta_val, 1, dhtlp_reparam, x=simulated_sample[k,])))
  #}
  value = -sum(apply(xi_probs_i_est_and_sample, 1, function(i) i[1:ncolor]%*%log(apply(theta_val, 1, dhtlp_reparam, x = i[(ncolor+1):(ncolor+2)]))), na.rm=T)
  if(is.na(value)){return(10000000)}
  if (value<=0.000){return(10000000)}else{return(value)}
  #return(-value)
}
