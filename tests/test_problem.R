n_grid_test = 30
ncolor_test = 2
rho_test = 0.5
potts_param_test <- c(rep(0, ncolor_test), rho_test)
x_potts_test <- matrix(1, nrow = n_grid_test, ncol = n_grid_test)
foo_test <- packPotts(x_potts_test, ncolor_test)
out_test <- potts(foo_test, potts_param_test, nbatch = 10)
image.plot(unpackPotts(out_test$final))

spat_pros_test = as.vector(unpackPotts(out_test$final))
y_test = rnorm(n_grid_test^2, spat_pros_test*2-3, rep(1,9))

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


likelihood <- function(parameters){
  rho = parameters[1]
  mu = parameters[2:(ncolor_test+1)]
  
  # First initiate normality constants and forward probabilities
  norm_const = 0
  for (j in 1:ncolor_test){
    norm_const = norm_const + dnorm(y_test[1], mu[j],1)*(1/ncolor_test)
  }
  norm_const = 1/norm_const
  norm_const = c(norm_const, rep(NA, n_grid_test^2-1))
  #norm_const = c(1/(dnorm(y_test[1],mu[1],1)*0.5+dnorm(y_test[1],mu[2],1)*0.5), rep(NA,8))
  forward_prob = matrix(0,nrow=n_grid_test^2,ncol = ncolor_test)
  for (j in 1:ncolor_test){
    forward_prob[1,j] = norm_const[1]*(1/ncolor_test)*dnorm(y_test[1], mu[j], 1)
  }
  #forward_prob[1,1] = norm_const[1]*0.5*dnorm(y_test[1],mu[1],1)
  #forward_prob[1,2] = norm_const[1]*0.5*dnorm(y_test[1],mu[2],1)
  transition = matrix(0, nrow = n_grid_test^2, ncol = ncolor_test^3)
  potts_prob_2 = apply(xi_A_2, 1, function(i) exp(rho*ifelse(i[1]==i[2],1,0)))
  potts_prob_2 = potts_prob_2/sum(potts_prob_2)
  potts_prob_3 = apply(xi_A_3, 1, function(i) exp(rho*(length(which(i[1:2]==i[3])))))
  potts_prob_3 = potts_prob_3/sum(potts_prob_3)
  
  # Iterate to find constants and probabilities
  for (t in 2:n_grid_test){
    #norm_const[t] = 1/(dnorm(y_test[t],0,tau)*p*forward_prob[t-1,1] + (1-p)*dnorm(y_test[t],0,tau)*forward_prob[t-1,2] + dnorm(y_test[t],1,tau)*(1-p)*forward_prob[t-1,1] + dnorm(y_test[t],1,tau)*p*forward_prob[t-1,2])
    for (j in 1:(ncolor_test^2)){
      transition[t,j] = potts_prob_2[j]*dnorm(y_test[t],mu[xi_A_2[j,2]],1)*forward_prob[t-1,xi_A_2[j,1]]}
    #transition[t,1] = potts_prob_2[1]*dnorm(y_test[t],mu[1],1)*forward_prob[t-1,1]
    #transition[t,2] = potts_prob_2[2]*dnorm(y_test[t],mu[2],1)*forward_prob[t-1,1]
    #transition[t,3] = potts_prob_2[3]*dnorm(y_test[t],mu[1],1)*forward_prob[t-1,2]
    #transition[t,4] = potts_prob_2[4]*dnorm(y_test[t],mu[2],1)*forward_prob[t-1,2]
    norm_const[t] = 1/sum(transition[t,])
    transition[t,] = transition[t,]*norm_const[t]
    for (j in 1:ncolor_test){
      forward_prob[t,j] = sum(transition[t,(0:(ncolor_test-1))*ncolor_test+j])
    }
    #forward_prob[t,1] = transition[t,1] + transition[t,3]
    #forward_prob[t,2] = transition[t,2] + transition[t,4]
  }
  for (t in (n_grid_test+1):(n_grid_test^2)){
    if ((t-1)%%n_grid_test==0){
      for (j in 1:(ncolor_test^2)){
        transition[t,j] = potts_prob_2[j]*dnorm(y_test[t],mu[xi_A_2[j,2]],1)*forward_prob[t-n_grid_test,xi_A_2[j,1]]}
      #transition[t,1] = potts_prob_2[1]*dnorm(y_test[t],mu[1],1)*forward_prob[t-1,1]
      #transition[t,2] = potts_prob_2[2]*dnorm(y_test[t],mu[2],1)*forward_prob[t-1,1]
      #transition[t,3] = potts_prob_2[3]*dnorm(y_test[t],mu[1],1)*forward_prob[t-1,2]
      #transition[t,4] = potts_prob_2[4]*dnorm(y_test[t],mu[2],1)*forward_prob[t-1,2]
      norm_const[t] = 1/sum(transition[t,])
      transition[t,] = transition[t,]*norm_const[t]
      for (j in 1:ncolor_test){
        forward_prob[t,j] = sum(transition[t,(0:(ncolor_test-1))*ncolor_test+j])
      }
      #forward_prob[t,1] = transition[t,1] + transition[t,3]
      #forward_prob[t,2] = transition[t,2] + transition[t,4]
    }else{
      for (j in 1:(ncolor_test^3)){
        transition[t,j] = potts_prob_3[j]*dnorm(y_test[t],mu[xi_A_3[j,3]],1)*forward_prob[t-1,xi_A_3[j,2]]*forward_prob[t-n_grid_test, xi_A_3[j,1]]
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
  return(-sum(log(norm_const)))
  #return(norm_const)
}

neg_likelihood <- function(parameters){
  return(-likelihood(parameters))
}



time_test = Sys.time()
optimal = optim(c(0.3,-1,1), neg_likelihood, method = "BFGS", control = list(trace=6, REPORT = 1))
#likelihood(c(0.5,-5,5,15,25, 35))
Sys.time() - time_test
optimal$par
