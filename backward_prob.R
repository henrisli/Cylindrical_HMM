backward_probs <- function(parameters_input, n_rows, data_sample, n_cols){
  rho_val = parameters_input[1]
  theta_val = matrix(parameters_input[2:(5*ncolor_test+1)], nrow=ncolor_test)
  
  # First initiate normality constants and forward probabilities
  norm_const = 0
  for (j in 1:ncolor_test){
    norm_const = norm_const + (1/ncolor_test)*dabeley_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  norm_const = 1/norm_const
  norm_const = c(norm_const, rep(NA, n_rows*n_cols-1))
  forward_prob = matrix(0,nrow=n_rows*n_cols,ncol = ncolor_test)
  for (j in 1:ncolor_test){
    forward_prob[1,j] = norm_const[1]*(1/ncolor_test)*dabeley_reparam(theta_val[j,], as.vector(data_sample[1,]))
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
      transition[t,j] = potts_prob_2[(j-1)%/%(ncolor_test^(t-2))+1]*moving_window[t-1,(j-1)%%(ncolor_test^(t-1))+1]*dabeley_reparam(theta_val[xi_A_n[j,t],], as.vector(data_sample[t,]))}
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
        transition[t,j] = potts_prob_2[conf_n_2[j]]*dabeley_reparam(theta_val[xi_A_n[j,n_rows+1],], as.vector(data_sample[t,]))*moving_window[t-1,(j-1)%%(ncolor_test^(n_rows))+1]
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

backward_probs_1 <- function(parameters_input, n_rows, data_sample, n_cols){
  rho_val = parameters_input[1]
  theta_val = matrix(parameters_input[2:(5*ncolor_test+1)], nrow=ncolor_test)
  
  # First initiate normality constants and forward probabilities
  norm_const = 0
  for (j in 1:ncolor_test){
    norm_const = norm_const + (1/ncolor_test)*dabeley_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  norm_const = 1/norm_const
  norm_const = c(norm_const, rep(NA, n_rows*n_cols-1))
  forward_prob = matrix(0,nrow=n_rows*n_cols,ncol = ncolor_test)
  for (j in 1:ncolor_test){
    forward_prob[1,j] = norm_const[1]*(1/ncolor_test)*dabeley_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  transition = matrix(0, nrow = n_cols*n_rows, ncol = ncolor_test^(n_rows+1))
  moving_window = matrix(0, nrow = n_cols*n_rows, ncol = ncolor_test^(n_rows))
  moving_window[1, 1:ncolor_test] = forward_prob[1,]
  potts_prob_2 = apply(xi_A_2, 1, function(i) exp(rho_val*ifelse(i[1]==i[2],1,0)))
  potts_prob_2 = potts_prob_2/sum(potts_prob_2)
  
  # Iterate to find constants and probabilities
  for (t in 2:n_cols){
    for (j in 1:(ncolor_test^2)){
      transition[t,j] = potts_prob_2[j]*moving_window[t-1,(j-1)%%(ncolor_test^(2-1))+1]*dabeley_reparam(theta_val[xi_A_n[j,2],], as.vector(data_sample[t,]))}
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


# backward_probabilities = matrix(0, nrow = n_cols^2, ncol = ncolor_test)
# for (i in 1:(n_cols/n_rows)){
#   backward_probabilities[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),] = backward_probs_1(init_param, n_rows, simulated_sample_temp[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),], n_cols)
# }

find_back_probs <- function(parameters_input, n_rows, data_sample, n_cols){
  backward_probabilities = matrix(0, nrow = n_cols^2, ncol = ncolor_test)
  for (i in 1:(n_cols/n_rows)){
    if (n_rows==1){
      backward_probabilities[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),] = backward_probs_1(parameters_input, n_rows, data_sample[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),], n_cols)
    }else{backward_probabilities[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),] = backward_probs(parameters_input, n_rows, data_sample[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),], n_cols)}
  }
  reverse = as.vector(t(matrix(1:(n_cols^2), nrow=n_cols)))
  data_sample = data_sample[reverse,]
  for (i in 1:(n_cols/n_rows)){
    if (n_rows==1){
      backward_probabilities[reverse[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i)],] = backward_probabilities[reverse[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i)],] + backward_probs_1(parameters_input, n_rows, data_sample[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),], n_cols)
    }else{backward_probabilities[reverse[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i)],] = backward_probabilities[reverse[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i)],] + backward_probs(parameters_input, n_rows, data_sample[((i-1)*n_cols*n_rows+1):(n_rows*n_cols*i),], n_cols)}
  }
  return(backward_probabilities/2)
}



# Simulate from the forward-backward algorithm
simulate_backward_probs_1 <- function(parameters_input, n_rows, data_sample, n_cols){
  rho_val = parameters_input[1]
  theta_val = matrix(parameters_input[2:(5*ncolor_test+1)], nrow=ncolor_test)
  
  # First initiate normality constants and forward probabilities
  norm_const = 0
  for (j in 1:ncolor_test){
    norm_const = norm_const + (1/ncolor_test)*dabeley_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  norm_const = 1/norm_const
  norm_const = c(norm_const, rep(NA, n_rows*n_cols-1))
  forward_prob = matrix(0,nrow=n_rows*n_cols,ncol = ncolor_test)
  for (j in 1:ncolor_test){
    forward_prob[1,j] = norm_const[1]*(1/ncolor_test)*dabeley_reparam(theta_val[j,], as.vector(data_sample[1,]))
  }
  transition = matrix(0, nrow = n_cols*n_rows, ncol = ncolor_test^(n_rows+1))
  moving_window = matrix(0, nrow = n_cols*n_rows, ncol = ncolor_test^(n_rows))
  moving_window[1, 1:ncolor_test] = forward_prob[1,]
  potts_prob_2 = apply(xi_A_2, 1, function(i) exp(rho_val*ifelse(i[1]==i[2],1,0)))
  potts_prob_2 = potts_prob_2/sum(potts_prob_2)
  
  # Iterate to find constants and probabilities
  for (t in 2:n_cols){
    for (j in 1:(ncolor_test^2)){
      transition[t,j] = potts_prob_2[j]*moving_window[t-1,(j-1)%%(ncolor_test^(2-1))+1]*dabeley_reparam(theta_val[xi_A_n[j,2],], as.vector(data_sample[t,]))}
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
