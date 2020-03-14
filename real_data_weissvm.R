library(ggpubr)
library(grid)
library(reshape2)
library(circular)
library(ggplot2)
library(fields)
library(potts)
library(rootSolve)

n_grid = 32
mat_grid <- matrix(seq(n_grid^2), n_grid)

addresses <- expand.grid(x = 1:n_grid, y = 1:n_grid)

# Relative addresses
z <- rbind(c(0,-1,1,0),c(-1,0,0,1))

get.neighbors <- function(rw) {
  # Convert to absolute addresses 
  z2 <- t(z + unlist(rw))
  # Choose those with indices within mat_grid 
  b.good <- rowSums(z2 > 0)==2  &  z2[,1] <= nrow(mat_grid)  &  z2[,2] <=ncol(mat_grid)
  mat_grid[z2[b.good,]]
}

# Not exactly sure if the simulation is working...

neighbor_list = apply(addresses,1, get.neighbors) # Returns a list with neighbors
k=4
ncolor = 4


dabeley <- function(param, x){
  alpha = param[1]
  beta=param[2]
  mu=param[3]
  kappa=param[4]
  lambda=param[5]
  return(alpha*beta^alpha/(2*pi*cosh(kappa))*(1+lambda*sin(x[2]-mu))*x[1]^(alpha-1)*exp(-(beta*x[1])^alpha*(1-tanh(kappa)*cos(x[2]-mu))))}
X_cor = seq(0,0.35,l=100)
y_cor = seq(-pi,pi,l=100)
vals = cbind(rep(X_cor,100), rep(y_cor,each=100))
parameters = rbind(c(1.4,10,0,0.6,0), c(1.9,15,0,0,0))
#parameters = rbind(c(0.5,0.1,0,0.9,0.5), c(0.8,0.9,1.5,0.9,0.5), c(1.3,1.7,3,0.9,0.5))
#parameters = rbind(c(0.5,0.3,0,0.75,0.5), c(0.75,0.6,1,0.75,0.5), c(1,0.9,2,0.75,0.5))
values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[1,])
values[which(values==Inf)]=0
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[2,])
values[which(values==Inf)]=0
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))

simulated_sample = as.matrix(data_set)


# Estimate rho by maximum pseudo likelihood if we consider the hidden markov RF to be known

neighbor_matrix = list(matrix(0,n_grid,n_grid))
for (mat_count in 2:k){
  neighbor_matrix[[mat_count]] = matrix(0,n_grid,n_grid)}
for (grid_count in 1:n_grid^2){
  neighbor_matrix[[spat_pros[grid_count]]][neighbor_list[[grid_count]]] =  neighbor_matrix[[spat_pros[grid_count]]][neighbor_list[[grid_count]]]+1
}
equal.neighbors = rep(NA,n_grid^2)
for(count in 1:n_grid^2){
  equal.neighbors[count] = neighbor_matrix[[spat_pros[count]]][count]
}
loglikelihood <- function(beta){
  result = sum(beta*equal.neighbors - log(exp(beta*as.vector(neighbor_matrix[[1]])) + exp(beta*as.vector(neighbor_matrix[[2]])) + exp(beta*as.vector(neighbor_matrix[[3]]))))
  return(-result)
}
optimal = optim(1,loglikelihood, lower = 0, method = "Brent", upper = 2*log(1+sqrt(ncolor)))
optimal$par


ncolor = 4

# Initialize parameters

iteration = 1
A = n_grid*(n_grid-1)*2
A_list = list()
i = 1
for (j in 1:n_grid){
  for (k in 1:(n_grid-1)){
    A_list[[i]] = c((j-1)*n_grid + k, (j-1)*n_grid + k + 1)
    i = i+1
  }
}
for (j in 1:n_grid){
  for (k in 1:(n_grid-1)){
    A_list[[i]] = c((k-1)*n_grid + j, k*n_grid + j)
    i = i+1
  }
}
xi_probs = matrix(NA, ncol = ncolor^2, nrow = A)
xi_A = matrix(NA,nrow = ncolor^2,ncol = 2)
it = 1
for (row in 1:ncolor){
  for (column in 1:ncolor){
    xi_A[it, ] = c(row,column)
    it = it + 1
  }
}

rho_function = function(xi_probs_est, rho_val){
  value = 0
  potts_prob = apply(xi_A, 1, function(i) exp(rho_val*ifelse(i[1]==i[2],1,0)))
  potts_prob = potts_prob/sum(potts_prob)
  return(-sum(apply(xi_probs_est, 1, function(i) i%*%log(potts_prob))))
  #for (i in 1:A){
  #  value = value + xi_probs_est[i,]%*%log(potts_prob)
  #for (j in 1:ncolor^2){
  #  value = value + xi_probs_est[i,j]*log(potts_prob[j])
  #}
  #}
  #return(-value)
}


dabeley_reparam <- function(param, x){
  alpha = exp(param[1])
  beta=exp(param[2])
  mu=2*atan(param[3])
  kappa=exp(param[4])
  lambda=tanh(param[5])
  return(alpha*beta^alpha/(2*pi*cosh(kappa))*(1+lambda*sin(x[2]-mu))*x[1]^(alpha-1)*exp(-(beta*x[1])^alpha*(1-tanh(kappa)*cos(x[2]-mu))))
}



theta_function = function(xi_probs_i_est_and_sample, theta_val){
  theta_val = matrix(theta_val, nrow = ncolor)
  #value = 0
  #for (i in 1:A){
  #  j = A_list[[i]][1]
  #  k = A_list[[i]][2]
  #  value = value + sum(xi_probs_i_est_and_sample[j,1:ncolor]*log(apply(theta_val, 1, dabeley_reparam, x=simulated_sample[j,])))
  #  value = value + sum(xi_probs_i_est_and_sample[k,1:ncolor]*log(apply(theta_val, 1, dabeley_reparam, x=simulated_sample[k,])))
  #}
  value = -sum(apply(xi_probs_i_est_and_sample, 1, function(i) i[1:ncolor]%*%log(apply(theta_val, 1, dabeley_reparam, x = i[(ncolor+1):(ncolor+2)]))), na.rm=T)
  if(is.na(value)){return(10000000)}#else{return(value)}
  if (value<=0.000){return(10000000)}else{return(value)}
  #return(-value)
}

theta_function_converge = function(xi_probs_i_est_and_sample, theta_val){
  theta_val = matrix(theta_val, nrow = ncolor)
  #value = 0
  #for (i in 1:A){
  #  j = A_list[[i]][1]
  #  k = A_list[[i]][2]
  #  value = value + sum(xi_probs_i_est_and_sample[j,1:ncolor]*log(apply(theta_val, 1, dabeley_reparam, x=simulated_sample[j,])))
  #  value = value + sum(xi_probs_i_est_and_sample[k,1:ncolor]*log(apply(theta_val, 1, dabeley_reparam, x=simulated_sample[k,])))
  #}
  value = -sum(apply(xi_probs_i_est_and_sample, 1, function(i) i[1:ncolor]%*%log(apply(theta_val, 1, dabeley_reparam, x = i[(ncolor+1):(ncolor+2)]))), na.rm=T)
  if(is.na(value)){return(10000000)}else{return(value)}
  #if (value<=0.000){return(10000000)}else{return(value)}
  #return(-value)
}

full_likelihood = function(parameter){
  rho_val = parameter[1]
  theta_val = matrix(parameter[2:(5*ncolor+1)], nrow=ncolor)
  potts_prob = apply(xi_A, 1, function(i) exp(rho_val*ifelse(i[1]==i[2],1,0)))
  potts_prob = potts_prob/sum(potts_prob)
  value = 0
  for (i in 1:A){
    temp = 0
    #value = value + log(potts_prob%*%apply(xi_A, 1, function(j) dabeley_reparam(theta_val[j[1],], as.vector(simulated_sample[A_list[[i]][1],]))*dabeley_reparam(theta_val[j[2],], as.vector(simulated_sample[A_list[[i]][2],]))))
    for (j in 1:ncolor^2){
      temp = temp + potts_prob[j]*dabeley_reparam(theta_val[xi_A[j,1],], as.vector(simulated_sample[A_list[[i]][1],]))*dabeley_reparam(theta_val[xi_A[j,2],], as.vector(simulated_sample[A_list[[i]][2],]))
    }
    value = value + log(temp)
  }
  return(-value)
}


#composite_likelihood = list(Inf)
ttime = Sys.time()
n_start = 10
composite_likelihood<- rep(1000000, n_start)
rho_vec = runif(n_start, 0, log(1+sqrt(ncolor)))
theta_list = list()
for (iter in 1:n_start){
  theta_list[[iter]] = matrix(NA, nrow=ncolor, ncol = 5)
  theta_list[[iter]][,1] = runif(ncolor,0.9,3) # 0.9-5 x100, 0.9-3 x1
  theta_list[[iter]][,2] = runif(ncolor,5,25) # 0-0.5 x100, 5-25 x1
  theta_list[[iter]][,3] = runif(ncolor,-pi,pi) # -pi - pi
  theta_list[[iter]][,4] = runif(ncolor,0,3) # 0 - 3
  theta_list[[iter]][,5] = runif(ncolor,-0.9,0.9) # -0.9-0.9
}

converged = 1
theta_list_converged = list()
likelihood_converged = list()
for (start_point in 1:n_start){
  conv_exit = F
  while(T){
    exit = F
    
    theta_est = theta_list[[start_point]]
    rho_est = rho_vec[start_point]
    
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
    xi_probs_i_normal = matrix(0, nrow = n_grid^2, ncol = ncolor)
    for (i in 1:n_grid^2){xi_probs_i_normal[i,] = xi_probs_i[i,]/sum(xi_probs_i[i,])}
    iteration = iteration + 1
    if(any(is.na(xi_probs))){break}
    opt_rho = optim(rho_est,rho_function, xi_probs_est = xi_probs, method = "L-BFGS-B", lower = 0, upper = log(1+sqrt(ncolor)))
    
    rho_vec[start_point] = opt_rho$par
    theta_est_reparam=theta_est
    theta_est_reparam[,c(1,2,4)] = log(theta_est_reparam[,c(1,2,4)])
    theta_est_reparam[,3] = tan(theta_est_reparam[,3]/2)
    theta_est_reparam[,5] = atanh(theta_est_reparam[,5])
    
    opt_theta = tryCatch(
      optim(as.vector(theta_est_reparam), fn = theta_function, method = "BFGS", xi_probs_i_est_and_sample = cbind(xi_probs_i,simulated_sample), control = list(reltol = 0.01)),
      error = function(e){ 
        exit = T
      }, finally = {}
    )
    if(class(opt_theta)!="list"){break}
    theta_est_new=matrix(opt_theta$par,nrow=ncolor)
    theta_est_new[,c(1,2,4)] = exp(theta_est_new[,c(1,2,4)])
    theta_est_new[,3] = 2*atan(theta_est_new[,3])
    theta_est_new[,5] = tanh(theta_est_new[,5])
    
    theta_list[[start_point]] = theta_est_new
    
    if(abs(full_likelihood(c(opt_rho$par, opt_theta$par)) - composite_likelihood[start_point])/(composite_likelihood[start_point])<0.01){
      composite_likelihood[start_point] = full_likelihood(c(opt_rho$par, opt_theta$par))
      conv_exit = T
      break}else{composite_likelihood[start_point] = full_likelihood(c(opt_rho$par, opt_theta$par))}
    # if(abs(composite_likelihood[[iteration]] - composite_likelihood[[iteration-1]])<10){break}
  }
  if (conv_exit){
    theta_list_converged[[converged]] = c(opt_rho$par, opt_theta$par)
    likelihood_converged[[converged]] = composite_likelihood[start_point]
    converged = converged + 1
  }
  print(start_point)
}

for (start_point in 1:n_start){
  conv_exit = F
  while(T){
    exit = F
    
    theta_est = theta_list[[start_point]]
    rho_est = rho_vec[start_point]
    
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
    xi_probs_i_normal = matrix(0, nrow = n_grid^2, ncol = ncolor)
    for (i in 1:n_grid^2){xi_probs_i_normal[i,] = xi_probs_i[i,]/sum(xi_probs_i[i,])}
    iteration = iteration + 1
    if(any(is.na(xi_probs))){break}
    opt_rho = optim(rho_est,rho_function, xi_probs_est = xi_probs, method = "L-BFGS-B", lower = 0, upper = log(1+sqrt(ncolor)))
    
    rho_vec[start_point] = opt_rho$par
    theta_est_reparam=theta_est
    theta_est_reparam[,c(1,2,4)] = log(theta_est_reparam[,c(1,2,4)])
    theta_est_reparam[,3] = tan(theta_est_reparam[,3]/2)
    theta_est_reparam[,5] = atanh(theta_est_reparam[,5])
    if (theta_function(as.vector(theta_est_reparam), xi_probs_i_est_and_sample = cbind(xi_probs_i,simulated_sample))!=1e7){break}
    opt_theta = tryCatch(
      optim(as.vector(theta_est_reparam), fn = theta_function_converge, method = "BFGS", xi_probs_i_est_and_sample = cbind(xi_probs_i,simulated_sample), control = list(reltol = 0.01)),
      error = function(e){ 
        exit = T
      }, finally = {}
    )
    if(class(opt_theta)!="list"){break}
    theta_est_new=matrix(opt_theta$par,nrow=ncolor)
    theta_est_new[,c(1,2,4)] = exp(theta_est_new[,c(1,2,4)])
    theta_est_new[,3] = 2*atan(theta_est_new[,3])
    theta_est_new[,5] = tanh(theta_est_new[,5])
    
    theta_list[[start_point]] = theta_est_new
    
    if(abs(full_likelihood(c(opt_rho$par, opt_theta$par)) - composite_likelihood[start_point])/(composite_likelihood[start_point])<0.01){
      composite_likelihood[start_point] = full_likelihood(c(opt_rho$par, opt_theta$par))
      conv_exit = T
      break}else{composite_likelihood[start_point] = full_likelihood(c(opt_rho$par, opt_theta$par))}
    # if(abs(composite_likelihood[[iteration]] - composite_likelihood[[iteration-1]])<10){break}
  }
  if (conv_exit){
    theta_list_converged[[converged]] = c(opt_rho$par, opt_theta$par)
    likelihood_converged[[converged]] = composite_likelihood[start_point]
    converged = converged + 1
  }
  print(start_point)
}

#plot(sapply(theta_list, function(i) sum(abs(i-parameters)^2)))
#plot(sapply(theta_list, function(i) sum(abs(i-parameters))))
plot(composite_likelihood)
min.val = which.min(composite_likelihood)
theta_iter = list(theta_list[[min.val]])
rho_iter = list(rho_vec[min.val])
composite_likelihood_iter = list(composite_likelihood[min.val])
iteration = 1
rho_est = rho_iter[[iteration]]
while(T){
  theta_est = theta_iter[[iteration]]
  rho_est = rho_iter[[iteration]]
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
  xi_probs_i_normal = matrix(0, nrow = n_grid^2, ncol = ncolor)
  for (i in 1:n_grid^2){xi_probs_i_normal[i,] = xi_probs_i[i,]/sum(xi_probs_i[i,])}
  #image.plot(matrix(apply(xi_probs_i_normal,1,function(i) max(i)*(1-max(i))), nrow = n_grid))
  iteration = iteration + 1
  opt_rho = optim(rho_est,rho_function, xi_probs_est = xi_probs, method = "Brent", lower = 0, upper = log(1+sqrt(ncolor)))
  rho_iter[[iteration]] = opt_rho$par
  print(opt_rho$par)
  theta_est_reparam=theta_est
  theta_est_reparam[,c(1,2,4)] = log(theta_est_reparam[,c(1,2,4)])
  theta_est_reparam[,3] = tan(theta_est_reparam[,3]/2)
  theta_est_reparam[,5] = atanh(theta_est_reparam[,5])
  opt_theta = optim(as.vector(theta_est_reparam), fn = theta_function, method = "BFGS", xi_probs_i_est_and_sample = cbind(xi_probs_i,simulated_sample), control = list(trace=6, REPORT=1, reltol = 1e-5))
  #opt_theta = optim(as.vector(theta_est_reparam), fn = theta_function, xi_probs_i_est_and_sample = cbind(xi_probs_i,simulated_sample), control = list(reltol = 1e-4))
  theta_est_new=matrix(opt_theta$par,nrow=ncolor)
  theta_est_new[,c(1,2,4)] = exp(theta_est_new[,c(1,2,4)])
  theta_est_new[,3] = 2*atan(theta_est_new[,3])
  theta_est_new[,5] = tanh(theta_est_new[,5])
  theta_iter[[iteration]] = theta_est_new
  composite_likelihood_iter[[iteration]] = full_likelihood(c(opt_rho$par, opt_theta$par))
  #composite_likelihood_iter[[iteration]] = rho_function(xi_probs_est = xi_probs, rho_val = rho_est) + opt_theta$value
  print(composite_likelihood_iter[[iteration]])
  if(abs(composite_likelihood_iter[[iteration]] - composite_likelihood_iter[[iteration-1]])/composite_likelihood_iter[[iteration-1]]<0.000001){break}
}
#opt_test = optim(c(opt_rho$par, opt_theta$par), method = "BFGS", fn = full_likelihood, control = list(trace = 6, REPORT = 1, reltol =1e-5))
opt_test = optim(theta_list_converged[[which.min(likelihood_converged)]], method = "BFGS", fn = full_likelihood, control = list(trace = 6, REPORT = 1, reltol =1e-5))
#opt_test_rho = optim(opt_rho$par, method = "BFGS", fn = full_likelihood_rho_only, control = list(trace = 6, REPORT = 1, reltol =1e-5), parameter = opt_theta$par)
rho_estimate = opt_test$par[1+5*ncolor]
theta_estimate = matrix(opt_test$par[1:(5*ncolor)], ncol = 5)
theta_estimate[,c(1,2,4)] = exp(theta_estimate[,c(1,2,4)])
theta_estimate[,3]=2*atan(theta_estimate[,3])
theta_estimate[,5]=tanh(theta_estimate[,5])
Sys.time() - ttime

grad_1 = numDeriv::grad(full_likelihood, c(opt_rho$par, opt_theta$par))
hess_1 = numDeriv::hessian(full_likelihood, c(opt_rho$par, opt_theta$par))
sum(diag(grad_1%*%t(grad_1)%*%solve(hess_1)))
grad_2 = numDeriv::grad(full_likelihood, c(opt_test_rho$par, opt_theta$par))
hess_2 = numDeriv::hessian(full_likelihood, c(opt_test_rho$par, opt_theta$par))
sum(diag((grad_2)%*%t(grad_2)%*%solve(hess_2)))
grad_3 = numDeriv::grad(full_likelihood, opt_test$par)
hess_3 = numDeriv::hessian(full_likelihood, opt_test$par)
sum(diag((grad_3)%*%t(grad_3)%*%solve(hess_3)))
grad_4 = rootSolve::gradient(theta_function, opt_theta$par, xi_probs_i_est_and_sample = cbind(xi_probs_i,simulated_sample))
#hess_4 = jacobian(gradient_func, opt_theta$par, xi_probs_i_est_and_sample = cbind(xi_probs_i,simulated_sample))
sum(diag((grad_4)%*%t(grad_4)%*%solve(hess_4)))



pdf(file="C:/Users/henri/Documents/GitHub/Master-Thesis/Images/Case1_latent_composite.pdf")
image(matrix(apply(xi_probs_i_normal,1,which.max),nrow=n_grid), x = 1:24, y = 1:24, xlab = "", ylab = "", col = tim.colors(64))
dev.off()
plot_image(apply(xi_probs_i_normal,1,which.max))
plot_image(spat_pros_test)
spat_pros_2 = spat_pros_test
spat_pros_2[which(spat_pros_test==1)] = 2
spat_pros_2[which(spat_pros_test==2)] = 3
spat_pros_2[which(spat_pros_test==3)] = 1
plot_image(spat_pros_2)

# 3: blue, 2: green, 1: red
estimated_xi = apply(xi_probs_i_normal,1,which.max)
ggplot(data.frame(x = simulated_sample[,1], theta = simulated_sample[,2])) + geom_point(aes(x=x, y = theta, col = as.factor(estimated_xi), shape = as.factor(spat_pros))) + theme_bw() + theme(legend.position = "none")

plot(sapply(theta_iter, function(i) mean(abs(i-parameters)^2)))
plot(sapply(theta_iter, function(i) mean(abs(i-parameters))))
plot(unlist(composite_likelihood_iter))

theta_max = theta_iter[[iteration]]
theta_max = theta_estimate

X_cor = seq(0,5,l=100)
y_cor = seq(-pi,pi,l=100)
vals = cbind(rep(X_cor,100), rep(y_cor,each=100))
for (i in 1:ncolor){
  values = apply(vals, MARGIN= 1, FUN = dabeley, param=theta_max[i,])
  image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
  if(i<=3){
    values = apply(vals, MARGIN= 1, FUN = dabeley, param=parameters[i,])
    image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))}}

# With only one class

loglik <- function(param, X, theta){
  alpha = param[1]
  beta = param[2]
  mu = param[3]
  kappa = param[4]
  lambda = param[5]
  if(alpha<0|beta<0|mu<(-pi)|mu>pi|kappa<0|abs(lambda)>1){return(Inf)}
  return(-((alpha-1)*sum(log(X)) - beta^alpha*sum(X^alpha*(1-tanh(kappa)*cos(theta-mu)))+sum(log(1+lambda*sin(theta-mu)))+n_grid*(alpha*log(beta)+log(alpha)-log(2*pi*cosh(kappa)))))
}

init_par = c(1.4,10,0,0.6,0)
ans = optim(init_par, loglik, X = simulated_sample[,1], theta = simulated_sample[,2], control = list(trace=6, maxit=10000))
values = apply(vals, MARGIN= 1, FUN = dabeley, param=ans$par)
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))




# Exact likelihood

library(prodlim)
n_rows = 1
n_cols = 32
ncolor_test = 3
# xi_A_n = matrix(1:2, ncol = 1)
# for (i in 1:(n_rows)){
#   xi_A_n = cbind(rbind(xi_A_n, xi_A_n), rep(1:2, each = 2^i))
# }
# xi_A_n = matrix(1:3, ncol = 1)
# for (i in 1:(n_rows)){
#   xi_A_n = cbind(rbind(xi_A_n, xi_A_n, xi_A_n), rep(1:3, each = 3^i))
# }
xi_A_n = matrix(1:4, ncol = 1)
for (i in 1:(n_rows)){
  xi_A_n = cbind(rbind(xi_A_n, xi_A_n,xi_A_n, xi_A_n), rep(1:4, each = 4^i))
}


parameters_test_reparam = c(rho_est, as.vector(theta_iter[[1]]))

# parameters_test_reparam[c(2,3,4,5,8,9)] = log(parameters_test_reparam[c(2,3,4,5,8,9)])
# parameters_test_reparam[c(6,7)] = atan(parameters_test_reparam[c(6,7)]/2)
# parameters_test_reparam[c(10,11)] = atanh(parameters_test_reparam[c(10,11)])

parameters_test_reparam[c(2,3,4,5,6,7,11,12,13)] = log(parameters_test_reparam[c(2,3,4,5,6,7,11,12,13)])
parameters_test_reparam[c(8,9,10)] = atan(parameters_test_reparam[c(8,9,10)]/2)
parameters_test_reparam[c(14,15,16)] = atanh(parameters_test_reparam[c(14,15,16)])

# parameters_test_reparam[c(2,3,4,5,6,7,8,9,14,15,16,17)] = log(parameters_test_reparam[c(2,3,4,5,6,7,8,9,14,15,16,17)])
# parameters_test_reparam[c(10,11,12,13)] = atan(parameters_test_reparam[c(10,11,12,13)]/2)
# parameters_test_reparam[c(18,19,20,21)] = atanh(parameters_test_reparam[c(18,19,20,21)])

init_param = parameters_test_reparam

optimal = optim(init_param, neg_likelihood_exact, method = "BFGS", control = list(trace=6, REPORT = 1, reltol = 1e-5), n_rows = n_rows, data_sample = simulated_sample, n_cols = n_cols)

estimated_param = rep(optimal$par[1],1+5*ncolor_test)

#estimated_param[c(2,3,4,5,8,9)] = exp(optimal$par[c(2,3,4,5,8,9)])
#estimated_param[c(6,7)] = 2*atan(optimal$par[c(6,7)])
#estimated_param[c(10,11)] = tanh(optimal$par[c(10,11)])

estimated_param[c(2,3,4,5,6,7,11,12,13)] = exp(optimal$par[c(2,3,4,5,6,7,11,12,13)])
estimated_param[c(8,9,10)] = 2*atan(optimal$par[c(8,9,10)])
estimated_param[c(14,15,16)] = tanh(optimal$par[c(14,15,16)])

# estimated_param[c(2,3,4,5,6,7,8,9,14,15,16,17)] = exp(optimal$par[c(2,3,4,5,6,7,8,9,14,15,16,17)])
# estimated_param[c(10,11,12,13)] = 2*atan(optimal$par[c(10,11,12,13)])
# estimated_param[c(18,19,20,21)] = tanh(optimal$par[c(18,19,20,21)])

estimated_param[1]
estimated_param = matrix(estimated_param[2:(5*ncolor_test+1)],nrow=ncolor_test)
estimated_param


estimated_probabilities = find_back_probs(optimal$par, n_rows, simulated_sample, n_cols)

image(matrix(apply(estimated_probabilities,1,which.max),nrow=n_grid), x = 1:n_grid, y = 1:n_grid, xlab = "", ylab = "", col = tim.colors(64))


values = apply(vals, MARGIN= 1, FUN = dabeley, param=estimated_param[1,])
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
values = apply(vals, MARGIN= 1, FUN = dabeley, param=estimated_param[2,])
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
values = apply(vals, MARGIN= 1, FUN = dabeley, param=estimated_param[3,])
values[which(values==Inf)]=0
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
values = apply(vals, MARGIN= 1, FUN = dabeley, param=estimated_param[4,])
values[which(values==Inf)]=0
image.plot(x=X_cor, y = y_cor, z = matrix(values,nrow=100))
