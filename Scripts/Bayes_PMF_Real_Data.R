#This script generates data from a Bayseian PMF distribution and then implements a Gibbs sampler
library(stats)

set.seed(1)

#load training set
load("../Data Structures/train_R_nm.rda")
#rename for cooperation with old code
R_nm = train_R_nm

#recover the number of users and items
N = nrow(R_nm) #number of users
M = ncol(R_nm) #number of items

#remove old name
rm(train_R_nm)

#load cluster_mat
load("../Data Structures/train_I_nm.rda")

#rename for cooperation with old code
I_nm = train_I_nm 

#remove old named data structures
rm(train_I_nm)

#generate a list of which items each user has ratd
user_rated_n = NULL
for(n in 1:N){
  user_rated_n[[n]] = which(I_nm[n,]>0)
}

#generate a list of which users have rated each item
item_rated_m = NULL
for(m in 1:M){
  item_rated_m[[m]] = which(I_nm[,m]>0)
}

#set D
D = 2

#set hyperparameters according to suggestions in paper/MATLAB code
mu_u0 = rep(0,D)
mu_v0 = rep(0,D)
v_0 = D
W_0 = diag(D)
W_0_inv = solve(W_0)
beta_0 = 2
alpha = 2

##############
# Initialize #
##############

#initialize Lambda_U at random
Lambda_U = rWishart(1,v_0,W_0)
Lambda_U = Lambda_U[,,1]
Lambda_U_inv = solve(Lambda_U)


#initiailze mu_u at random
mu_u = rnorm(D)

#initialize U at random
U = matrix(0,nrow = N,ncol = D)
for(i in 1:nrow(U)){
  U[i,] = rnorm(D)
}

#initialize Lambda_V at random
Lambda_V = rWishart(1,v_0,W_0)
Lambda_V = Lambda_V[,,1]
Lambda_V_inv = solve(Lambda_V)

#initiailze mu_u at random
mu_v = rnorm(D)

# initialize V at random
V = matrix(0,nrow = M,ncol = D)
for(i in 1:nrow(V)){
  V[i,] = rnorm(D)
}

#set number of iterations
warmup = 1000
iterations = 1000
n_it = warmup + iterations

#create data structures to store chains
chain_mu_u_id = matrix(0,nrow= iterations, ncol = D)
chain_lambda_u_idd = array(0,dim = c(iterations,D,D))
chain_U_ind = array(0,dim=c(iterations,N,D))
chain_mu_v_id = matrix(0,nrow= iterations, ncol = D)
chain_lambda_v_idd = array(0,dim = c(iterations,D,D))
chain_V_imd = array(0,dim=c(iterations,M,D))

#################
# Gibbs Sampler #
#################

#Lets do the parameters independently first

for(it in 1:n_it){
  
  #####################
  # mu_u and Lambda_u #
  #####################
  
  #find sufficient stats first
  U_mean = apply(U,2,mean)
  S = matrix(0,nrow = ncol(U),ncol = ncol(U))
  for(n in 1:N){
    S = S + U[n,] %*% t(U[n,])
  }
  S = S/N
  
  #find full conditional parameters for mu_u and Lambda_u
  beta_star = beta_0 + N
  v_star = v_0 + N
  avg = (beta_0*mu_u0  + N*U_mean)/(beta_star)
  W_0_star_inv = W_0_inv + N*S + beta_0*N/beta_star * (U_mean - mu_u0) %*% t((U_mean - mu_u0))
  W_0_star = solve(W_0_star_inv)
  
  #draw Lambda_U first
  Lambda_U = rWishart(1,v_star,W_0_star)
  Lambda_U = Lambda_U[,,1]
  Lambda_U_inv = solve(Lambda_U)
  
  #draw mu_u from its full conditional
  mu_u = chol(1/beta_star*Lambda_U_inv) %*% rnorm(D) + avg
  
  #####
  # U #
  #####
  
  for(n in 1:N){
    
    #find sufficient stats
    S = matrix(0,nrow = D, ncol = D)
    
    for(m in user_rated_n[[n]]){
      S = S + V[m,] %*% t(V[m,])
    }
    
    s = rep(0,D)
    for(m in user_rated_n[[n]]){
      s = s + V[m,]*R_nm[n,m]
    }
    #maybe I should rename them becaue S and s are similar..
    
    #find variance and mean for normal
    covar = Lambda_U + alpha*S
    covar_inv = solve(covar)
    avg = covar_inv%*%(alpha*s + Lambda_U %*% mu_u )
    
    #draw U_n
    U[n,] = chol(covar_inv) %*% rnorm(D) + avg
  }
  
  #####################
  # mu_v and Lambda_v #
  #####################
  
  #find sufficient stats first
  V_mean = apply(V,2,mean)
  S = matrix(0,nrow = ncol(V),ncol = ncol(V))
  for(m in 1:M){
    S = S + V[m,] %*% t(V[m,])
  }
  S = S/M
  
  #find full conditional parameters for mu_u and Lambda_u
  beta_star = beta_0 + M
  v_star = v_0 + M
  avg = (beta_0*mu_u0  + M*V_mean)/(beta_star)
  W_0_star_inv = W_0_inv + M*S + beta_0*M/beta_star * (mu_u0 - V_mean) %*% t((mu_u0 - V_mean))
  W_0_star = solve(W_0_star_inv)
  
  #draw Lambda_V first
  Lambda_V = rWishart(1,v_star,W_0_star)
  Lambda_V = Lambda_V[,,1]
  Lambda_V_inv = solve(Lambda_V)
  
  #initiailze mu_v at random
  mu_v = chol(1/beta_star*Lambda_V_inv) %*% rnorm(D) + avg
  
  #####
  # V #
  #####
  
  for(m in 1:M){
    
    #find sufficient stats
    S = matrix(0,nrow = D, ncol = D)
    for(n in item_rated_m[[m]]){
      S = S + U[n,] %*% t(U[n,])
    }
    
    s = rep(0,D)
    for(n in item_rated_m[[m]]){
      s = s + U[n,]*R_nm[n,m]
    }
    #maybe I should rename them becaue S and s are similar..
    
    #find variance and mean for normal
    covar = Lambda_V + alpha*S
    covar_inv = solve(covar)
    avg = covar_inv%*%(alpha*s + Lambda_V %*% mu_v )
    
    #draw V_m
    V[m,] = chol(covar_inv) %*% rnorm(D) + avg
  }
  
  ##########
  # Record #
  ##########
  
  #record the sampled values
  if(it>warmup){
    #record mu_u
    chain_mu_u_id[it - warmup,] = mu_u
    
    chain_lambda_u_idd[it - warmup,,] = Lambda_U
    
    #record U
    chain_U_ind[it - warmup,,] = U
    
    #record mu_v
    chain_mu_v_id[it - warmup,] = mu_v
    
    chain_lambda_v_idd[it - warmup,,] = Lambda_V
    
    # record V
    chain_V_imd[it - warmup,,] = V
  }
  
  print(it)
  
}

######################
# posterior analysis #
######################

#posterior mean for mu_u
post_mu_u = apply(chain_mu_u_id[,],2,mean)
#chainplot
matplot(chain_mu_u_id,type="l")

#posterior mean for lambda_u
post_Lambda_U = matrix(0,nrow = D,ncol = D)
for(i in 1:nrow(post_Lambda_U)){
  for(j in 1:ncol(post_Lambda_U)){
    post_Lambda_U[i,j] = mean(chain_lambda_u_idd[,i,j])
  }
}
matplot(chain_lambda_u_idd[,1,],type="l")

#posterior mean for U
post_U_nd = matrix(0,nrow = N,ncol = D)
for(n in 1:N){
  post_U_nd[n,] = apply(chain_U_ind[,n,],2,mean)
}
#chainplot
matplot(chain_U_ind[,n,],type="l")

# posterior mean for mu_v
post_mu_v = apply(chain_mu_v_id,2,mean)
#chainplot
matplot(chain_mu_v_id,type="l")

#posterior mean for lambda_v
post_Lambda_V = matrix(0,nrow = D,ncol = D)
for(i in 1:nrow(post_Lambda_V)){
  for(j in 1:ncol(post_Lambda_V)){
    post_Lambda_V[i,j] = mean(chain_lambda_v_idd[,i,j])
  }
}

#posterior mean for V
post_V_md = matrix(0,nrow = M,ncol = D)
for(m in 1:M){
  post_V_md[m,] = apply(chain_V_imd[,m,],2,mean)
}
#chainplot
matplot(chain_V_imd[,m,],type="l")

#test how posterior means do on prediction on training set
post_R_nm = matrix(0,nrow = N, ncol = M)
for(n in 1:N){
  for(m in user_rated_n[[n]]){
    post_R_nm[n,m] = sum(post_U_nd[n,] * post_V_md[m,])
    # post_R_nm[n,m] = sum(true_U[n,] * true_V[m,])
  }
}
R_nm = R_nm[which(I_nm>0)]
post_R_nm = post_R_nm[which(I_nm>0)]
plot(R_nm[], post_R_nm[])

#######################
# test set prediction #
#######################

load("../Data Structures/test_R_nm.rda")
#rename for cooperation with old code
R_nm = test_R_nm

#remove old name
rm(test_R_nm)

#load cluster_mat
load("../Data Structures/test_I_nm.rda")

#rename for cooperation with old code
I_nm = test_I_nm 

#remove old named data structures
rm(test_I_nm)

#generate a list of which items each user has ratd
user_rated_n = NULL
for(n in 1:N){
  user_rated_n[[n]] = which(I_nm[n,]>0)
}

#generate a list of which users have rated each item
item_rated_m = NULL
for(m in 1:M){
  item_rated_m[[m]] = which(I_nm[,m]>0)
}

#test how posterior means do on prediction on training set
post_R_nm = matrix(0,nrow = N, ncol = M)
for(n in 1:N){
  for(m in user_rated_n[[n]]){
    post_R_nm[n,m] = sum(post_U_nd[n,] * post_V_md[m,])
  }
}
R_nm = R_nm[which(I_nm>0)]
post_R_nm = post_R_nm[which(I_nm>0)]
plot(R_nm[], post_R_nm[])

MAE = mean(abs(post_R_nm - R_nm))
MAE #this is 3.9435...
#also, for D = 2 we have mu_u => 0 while mu_v explodes
#somethigns not right here....
