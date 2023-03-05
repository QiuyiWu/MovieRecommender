#This script runs a k-means algorithm on real data, then runs a BPMF on each cluster
library(stats)

set.seed(1234)

#load training set
load("../Data Structures/train_R_nm.rda")

#recover the number of users and items
N = nrow(train_R_nm) #number of users
M = ncol(train_R_nm) #number of items

#load the training cluster matrix
load("../Data Structures/train_I_nm.rda")

#run the k-means algorithm for k = 2
K = 2
# out = kmeans(train_I_nm,centers = K, nstart = 25)

#find the cluster ids
# clust_id = out$cluster
clust_id = c(rep(1,500),rep(2,500))

# clust_id = rep(1,886)
# 
# #let's try on one cluster?
# train_R_nm = train_R_nm[which(clust_id==2),]
# train_I_nm = train_I_nm[which(clust_id==2),]
# 
# R_nm = train_R_nm
# I_nm = train_I_nm
# 
# N = nrow(R_nm)
# M = ncol(R_nm)
# 
# #remove old structures
# rm(train_R_nm)
# rm(train_I_nm)
# 
# #generate a list of which items each user has ratd
# user_rated_n = NULL
# for(n in 1:N){
#   user_rated_n[[n]] = which(I_nm[n,]>0)
# }
# 
# #generate a list of which users have rated each item
# item_rated_m = NULL
# for(m in 1:M){
#   item_rated_m[[m]] = which(I_nm[,m]>0)
# }
# 
# #set D
# D = 2
# 
# #set hyperparameters according to suggestions in paper/MATLAB code
# mu_u0 = rep(0,D)
# mu_v0 = rep(0,D)
# v_0 = D
# W_0 = diag(D)
# W_0_inv = solve(W_0)
# beta_0 = 2
# alpha = 2
# 
# ##############
# # Initialize #
# ##############
# 
# #initialize Lambda_U at random
# Lambda_U = rWishart(1,v_0,W_0)
# Lambda_U = Lambda_U[,,1]
# Lambda_U_inv = solve(Lambda_U)
# 
# 
# #initiailze mu_u at random
# mu_u = rnorm(D)
# 
# #initialize U at random
# U = matrix(0,nrow = N,ncol = D)
# for(i in 1:nrow(U)){
#   U[i,] = rnorm(D)
# }
# 
# #initialize Lambda_V at random
# Lambda_V = rWishart(1,v_0,W_0)
# Lambda_V = Lambda_V[,,1]
# Lambda_V_inv = solve(Lambda_V)
# 
# #initiailze mu_u at random
# mu_v = rnorm(D)
# 
# # initialize V at random
# V = matrix(0,nrow = M,ncol = D)
# for(i in 1:nrow(V)){
#   V[i,] = rnorm(D)
# }
# 
# #set number of iterations
# warmup = 1
# iterations = 1000
# n_it = warmup + iterations
# 
# #create data structures to store chains
# chain_mu_u_id = matrix(0,nrow= iterations, ncol = D)
# chain_lambda_u_idd = array(0,dim = c(iterations,D,D))
# chain_U_ind = array(0,dim=c(iterations,N,D))
# chain_mu_v_id = matrix(0,nrow= iterations, ncol = D)
# chain_lambda_v_idd = array(0,dim = c(iterations,D,D))
# chain_V_imd = array(0,dim=c(iterations,M,D))
# 
# #################
# # Gibbs Sampler #
# #################
# 
# #lets try something
# mean_rating = mean(R_nm[which(I_nm>0)])
# 
# for(it in 1:n_it){
#   
#   #####################
#   # mu_u and Lambda_u #
#   #####################
#   
#   #find sufficient stats first
#   U_mean = apply(U,2,mean)
#   S = matrix(0,nrow = ncol(U),ncol = ncol(U))
#   for(n in 1:N){
#     S = S + U[n,] %*% t(U[n,])
#   }
#   S = S/N
#   
#   #find full conditional parameters for mu_u and Lambda_u
#   beta_star = beta_0 + N
#   v_star = v_0 + N
#   avg = (beta_0*mu_u0  + N*U_mean)/(beta_star)
#   W_0_star_inv = W_0_inv + N*S + beta_0*N/beta_star * (U_mean - mu_u0) %*% t((U_mean - mu_u0))
#   W_0_star = solve(W_0_star_inv)
#   
#   #draw Lambda_U first
#   Lambda_U = rWishart(1,v_star,W_0_star)
#   Lambda_U = Lambda_U[,,1]
#   Lambda_U_inv = solve(Lambda_U)
#   
#   #draw mu_u from its full conditional
#   mu_u = chol(1/beta_star*Lambda_U_inv) %*% rnorm(D) + avg
#   
#   #####
#   # U #
#   #####
#   
#   for(n in 1:N){
#     
#     #find sufficient stats
#     S = matrix(0,nrow = D, ncol = D)
#     
#     for(m in user_rated_n[[n]]){
#       S = S + V[m,] %*% t(V[m,])
#     }
#     
#     s = rep(0,D)
#     for(m in user_rated_n[[n]]){
#       s = s + V[m,]*(R_nm[n,m]-mean_rating)
#     }
#     #maybe I should rename them becaue S and s are similar..
#     
#     #find variance and mean for normal
#     covar = Lambda_U + alpha*S
#     covar_inv = solve(covar)
#     avg = covar_inv%*%(alpha*s + Lambda_U %*% mu_u )
#     
#     #draw U_n
#     U[n,] = chol(covar_inv) %*% rnorm(D) + avg
#   }
#   
#   #####################
#   # mu_v and Lambda_v #
#   #####################
#   
#   #find sufficient stats first
#   V_mean = apply(V,2,mean)
#   S = matrix(0,nrow = ncol(V),ncol = ncol(V))
#   for(m in 1:M){
#     S = S + V[m,] %*% t(V[m,])
#   }
#   S = S/M
#   
#   #find full conditional parameters for mu_u and Lambda_u
#   beta_star = beta_0 + M
#   v_star = v_0 + M
#   avg = (beta_0*mu_u0  + M*V_mean)/(beta_star)
#   W_0_star_inv = W_0_inv + M*S + beta_0*M/beta_star * (mu_u0 - V_mean) %*% t((mu_u0 - V_mean))
#   W_0_star = solve(W_0_star_inv)
#   
#   #draw Lambda_V first
#   Lambda_V = rWishart(1,v_star,W_0_star)
#   Lambda_V = Lambda_V[,,1]
#   Lambda_V_inv = solve(Lambda_V)
#   
#   #initiailze mu_v at random
#   mu_v = chol(1/beta_star*Lambda_V_inv) %*% rnorm(D) + avg
#   
#   #####
#   # V #
#   #####
#   
#   for(m in 1:M){
#     
#     #find sufficient stats
#     S = matrix(0,nrow = D, ncol = D)
#     for(n in item_rated_m[[m]]){
#       S = S + U[n,] %*% t(U[n,])
#     }
#     
#     s = rep(0,D)
#     for(n in item_rated_m[[m]]){
#       s = s + U[n,]*(R_nm[n,m]-mean_rating)
#     }
#     #maybe I should rename them becaue S and s are similar..
#     
#     #find variance and mean for normal
#     covar = Lambda_V + alpha*S
#     covar_inv = solve(covar)
#     avg = covar_inv%*%(alpha*s + Lambda_V %*% mu_v )
#     
#     #draw V_m
#     V[m,] = chol(covar_inv) %*% rnorm(D) + avg
#   }
#   
#   ##########
#   # Record #
#   ##########
#   
#   #record the sampled values
#   if(it>warmup){
#     #record mu_u
#     chain_mu_u_id[it - warmup,] = mu_u
#     
#     chain_lambda_u_idd[it - warmup,,] = Lambda_U
#     
#     #record U
#     chain_U_ind[it - warmup,,] = U
#     
#     #record mu_v
#     chain_mu_v_id[it - warmup,] = mu_v
#     
#     chain_lambda_v_idd[it - warmup,,] = Lambda_V
#     
#     # record V
#     chain_V_imd[it - warmup,,] = V
#   }
#   
#   print(it)
#   
# }
# 
# ######################
# # posterior analysis #
# ######################
# 
# #posterior mean for mu_u
# post_mu_u = apply(chain_mu_u_id[,],2,mean)
# #chainplot
# matplot(chain_mu_u_id,type="l")
# 
# #posterior mean for lambda_u
# post_Lambda_U = matrix(0,nrow = D,ncol = D)
# for(i in 1:nrow(post_Lambda_U)){
#   for(j in 1:ncol(post_Lambda_U)){
#     post_Lambda_U[i,j] = mean(chain_lambda_u_idd[,i,j])
#   }
# }
# matplot(chain_lambda_u_idd[,1,],type="l")
# 
# #posterior mean for U
# post_U_nd = matrix(0,nrow = N,ncol = D)
# for(n in 1:N){
#   post_U_nd[n,] = apply(chain_U_ind[,n,],2,mean)
# }
# #chainplot
# matplot(chain_U_ind[,n,],type="l")
# 
# # posterior mean for mu_v
# post_mu_v = apply(chain_mu_v_id,2,mean)
# #chainplot
# matplot(chain_mu_v_id,type="l")
# 
# #posterior mean for lambda_v
# post_Lambda_V = matrix(0,nrow = D,ncol = D)
# for(i in 1:nrow(post_Lambda_V)){
#   for(j in 1:ncol(post_Lambda_V)){
#     post_Lambda_V[i,j] = mean(chain_lambda_v_idd[,i,j])
#   }
# }
# 
# #posterior mean for V
# post_V_md = matrix(0,nrow = M,ncol = D)
# for(m in 1:M){
#   post_V_md[m,] = apply(chain_V_imd[,m,],2,mean)
# }
# #chainplot
# matplot(chain_V_imd[,m,],type="l")
# 
# #test how posterior means do on prediction on training set
# post_R_nm = matrix(0,nrow = N, ncol = M)
# for(n in 1:N){
#   for(m in user_rated_n[[n]]){
#     post_R_nm[n,m] = sum(post_U_nd[n,] * post_V_md[m,]) + mean_rating
#     # post_R_nm[n,m] = sum(true_U[n,] * true_V[m,])
#   }
# }
# R_nm = R_nm[which(I_nm>0)]
# post_R_nm = post_R_nm[which(I_nm>0)]
# plot(R_nm[], post_R_nm[])
# 
# #######################
# # test set prediction #
# #######################
# 
# load("../Data Structures/test_R_nm.rda")
# #rename for cooperation with old code
# R_nm = test_R_nm[which(clust_id==2),]
# 
# #remove old name
# rm(test_R_nm)
# 
# #load cluster_mat
# load("../Data Structures/test_I_nm.rda")
# 
# #rename for cooperation with old code
# I_nm = test_I_nm[which(clust_id==2),]
# 
# #remove old named data structures
# rm(test_I_nm)
# 
# #generate a list of which items each user has ratd
# user_rated_n = NULL
# for(n in 1:N){
#   user_rated_n[[n]] = which(I_nm[n,]>0)
# }
# 
# #generate a list of which users have rated each item
# item_rated_m = NULL
# for(m in 1:M){
#   item_rated_m[[m]] = which(I_nm[,m]>0)
# }
# 
# #test how posterior means do on prediction on training set
# post_R_nm = matrix(0,nrow = N, ncol = M)
# for(n in 1:N){
#   for(m in user_rated_n[[n]]){
#     post_R_nm[n,m] = sum(post_U_nd[n,] * post_V_md[m,]) + mean_rating
#   }
# }
# R_nm = R_nm[which(I_nm>0)]
# post_R_nm = post_R_nm[which(I_nm>0)]
# plot(R_nm[], post_R_nm[])
# 
# MAE = mean(abs(post_R_nm - R_nm))
# MAE 
# rMSE= sqrt(mean((post_R_nm-R_nm)^2))
# rMSE
# #for the first 114 users, we have an MAE of 0.935 and an rMSE of 1.20682
# 
# #If we k-means cluster, then the MAE for the large clsuter is 0.9726351 and rMSE is 1.248449
# #for the small cluster, MAE is 0.8563893 and rMSE is 1.127419
# #the implied overall MAE is 0.9593831

######################
# whole sampler code #
######################

#create a map from global ID to cluster ID
clust_map = matrix(0,nrow = N,ncol = 3)
clust_map[,1] = 1:N
clust_map[,2] = clust_id
#compute cluster ID
for(k in 1:K){
  clust_map[which(clust_id == k),3] = 1:length(which(clust_id==k))
}

#define K different training cluster matrices depending on cluster ID 
train_I_k_nm = NULL
train_I_k_nm[[K +1]] = matrix(0,nrow = 100, ncol = 100)
for(k in 1:K){
  train_I_k_nm[[k]] = train_I_nm[which(clust_id==k),]
}
train_I_k_nm[[K +1]] = NULL

rm(train_I_nm)

#define K different training matrices based on cluster ID and find mean ratings
train_R_k_nm = NULL
mean_rating_k = NULL
train_R_k_nm[[K+1]] = matrix(0,nrow = 100, ncol = 100)
for(k in 1:K){
  train_R_k_nm[[k]] = train_R_nm[which(clust_id==k),]
  mean_rating_k[[k]] = mean(train_R_k_nm[[k]][which(train_I_k_nm[[k]]>0)])
}
train_R_k_nm[[K+1]] = NULL

#remove old data structures
rm(train_R_nm)

#generate a list of which items each user has rated by cluster
user_rated_k_n = NULL
for(k in 1:K){
  user_rated_n = NULL
  train_I_nm = train_I_k_nm[[k]]
  for(n in 1:nrow(train_I_nm)){
    user_rated_n[[n]] = which(train_I_nm[n,]>0)
  }
  user_rated_k_n[[k]] = user_rated_n
}
rm(train_I_nm)
rm(user_rated_n)

#generate a list of which users have rated each item by cluster
item_rated_k_m = NULL
for(k in 1:K){
  item_rated_m = NULL
  train_I_nm = train_I_k_nm[[k]]
  for(m in 1:M){
    item_rated_m[[m]] = which(train_I_nm[,m]>0)
  }
  item_rated_k_m[[k]] = item_rated_m
}
rm(item_rated_m)


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
Lambda_U_kdd = array(0,dim = c(K,D,D))
for(k in 1:K){
  Lambda_U_kdd[k,,] = rWishart(1,v_0,W_0)
}
#precompute inverses
Lambda_U_inv_kdd = Lambda_U_kdd
for(k in 1:K){
  Lambda_U_inv_kdd[k,,] = solve(Lambda_U_kdd[k,,])
}

#initiailze mu_u at random
mu_u_kd = matrix(0,nrow = K, ncol = D)
for(k in 1:K){
  mu_u_kd[k,] = rnorm(D)
}

#initialize U at random
U_k_nd = NULL
for(k in 1:K){
  U_nd = matrix(0,nrow = nrow(train_R_k_nm[[k]]),ncol = D)
  for(n in 1:nrow(U_nd)){
    U_nd[n,] = rnorm(D)
  }
  U_k_nd[[k]] = U_nd
}

#initialize Lambda_V at random
Lambda_V_kdd = array(0,dim = c(K,D,D))
for(k in 1:K){
  Lambda_V_kdd[k,,] = rWishart(1,v_0,W_0)
}
#precompute inverses
Lambda_V_inv_kdd = Lambda_V_kdd
for(k in 1:K){
  Lambda_V_inv_kdd[k,,] = solve(Lambda_V_kdd[k,,])
}

#initiailze mu_v at random
mu_v_kd = matrix(0,nrow = K, ncol = D)
for(k in 1:K){
  mu_v_kd[k,] = rnorm(D)
}

#initialize V at random
V_k_md = NULL
for(k in 1:K){
  V_md = matrix(0,nrow = ncol(train_R_k_nm[[k]]),ncol = D)
  for(n in 1:nrow(V_md)){
    V_md[n,] = rnorm(D)
  }
  V_k_md[[k]] = V_md
}

#set number of iterations
warmup = 1000
iterations = 1000
n_it = warmup + iterations

#create data structures to store chains
chain_mu_u_kid = array(0,dim= c(K,iterations,D))

chain_U_k_ind = NULL
chain_U_k_ind[[K+1]] = array(0,dim = c(iterations, N, D))
for(k in 1:K){
  chain_U_k_ind[[k]] = array(0,dim = c(iterations, nrow(train_R_k_nm[[k]]), D))
}
chain_U_k_ind[[K+1]] = NULL

chain_mu_v_kid = array(0,dim= c(K,iterations,D))

chain_V_k_imd = NULL
chain_V_k_imd[[K+1]] = array(0,dim = c(iterations, N, D))
for(k in 1:K){
  chain_V_k_imd[[k]] = array(0,dim = c(iterations, ncol(train_R_k_nm[[k]]), D))
}
chain_V_k_imd[[K+1]] = NULL


#################
# Gibbs Sampler #
#################
for(it in 1:n_it){

  #cycle through clusters
  for(k in 1:K){

    ###############################
    # set params based on cluster #
    ###############################

    #U params
    Lambda_U = Lambda_U_kdd[k,,]
    Lambda_U_inv = Lambda_U_inv_kdd[k,,]
    mu_u = mu_u_kd[k,]
    U = U_k_nd[[k]]

    #V params
    Lambda_V = Lambda_V_kdd[k,,]
    Lambda_V_inv = Lambda_V_inv_kdd[k,,]
    mu_v = mu_v_kd[k,]
    V = V_k_md[[k]]

    #ratings matrix
    R_nm = train_R_k_nm[[k]]
    N = nrow(R_nm)
    M = ncol(R_nm)

    #which users/items have ratings
    user_rated_n = user_rated_k_n[[k]]
    item_rated_m = item_rated_k_m[[k]]
    
    #find mean rating
    mean_rating = mean_rating_k[[k]]

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
        s = s + V[m,]*(R_nm[n,m]-mean_rating)
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
        s = s + U[n,]*(R_nm[n,m]-mean_rating)
      }
      #maybe I should rename them becaue S and s are similar..

      #find variance and mean for normal
      covar = Lambda_V + alpha*S
      covar_inv = solve(covar)
      avg = covar_inv%*%(alpha*s + Lambda_V %*% mu_v )

      #draw V_m
      V[m,] = chol(covar_inv) %*% rnorm(D) + avg
    }

    #store updated parameter values

    #U params
    Lambda_U_kdd[k,,] = Lambda_U
    Lambda_U_inv_kdd[k,,] = Lambda_U_inv
    mu_u_kd[k,] = mu_u
    U_k_nd[[k]] = U

    #V params
    Lambda_V_kdd[k,,] = Lambda_V
    Lambda_V_inv_kdd[k,,] = Lambda_V_inv
    mu_v_kd[k,] = mu_v
    V_k_md[[k]] = V


  }#end cluster loop



  ##########
  # Record #
  ##########

  #record the sampled values
  if(it>warmup){
    for(k in 1:K){
      #record mu_u
      chain_mu_u_kid[k,it-warmup,] = mu_u_kd[k,]

      #record U
      chain_U_k_ind[[k]][it-warmup,,] = U_k_nd[[k]]

      #record mu_v
      chain_mu_v_kid[k,it-warmup,] = mu_v_kd[k,]

      #record V
      chain_V_k_imd[[k]][it-warmup,,] = V_k_md[[k]]
    }
  }

  print(it)

}

######################
# posterior analysis #
######################

#posterior mean for mu_u
post_mu_u_kd = matrix(0,nrow = K, ncol = D)
for(k in 1:K){
  post_mu_u_kd[k,] = apply(chain_mu_u_kid[k,,],2,mean)
}
#chainplot
matplot(chain_mu_u_kid[k,,],type="l")


#posterior mean for U
post_U_k_nd = NULL
post_U_k_nd[[K+1]] = matrix(0,nrow = N,ncol = D)
for(k in 1:K){
  post_U_k_nd[[k]] = matrix(0,nrow = nrow(train_R_k_nm[[k]]),ncol = D)
  for(n in 1:nrow(train_R_k_nm[[k]])){
    post_U_k_nd[[k]][n,] = apply(chain_U_k_ind[[k]][,n,],2,mean)
  }
}
post_U_k_nd[[K+1]] = NULL
#chainplot
matplot(chain_U_k_ind[[k]][,n,],type="l")

#posterior mean for mu_v
post_mu_v_kd = matrix(0,nrow = K, ncol = D)
for(k in 1:K){
  post_mu_v_kd[k,] = apply(chain_mu_v_kid[k,,],2,mean)
}
#chainplot
matplot(chain_mu_v_kid[k,,],type="l")

#posterior mean for V
post_V_k_md = NULL
post_V_k_md[[K+1]] = matrix(0,nrow = N,ncol = D)
for(k in 1:K){
  post_V_k_md[[k]] = matrix(0,nrow = ncol(train_R_k_nm[[k]]),ncol = D)
  for(m in 1:ncol(train_R_k_nm[[k]])){
    post_V_k_md[[k]][m,] = apply(chain_V_k_imd[[k]][,m,],2,mean)
  }
}
post_V_k_md[[K+1]] = NULL
#chainplot
matplot(chain_V_k_imd[[k]][,n,],type="l")

#load training set
load("../Data Structures/train_R_nm.rda")

#recover the number of users and items
N = nrow(train_R_nm) #number of users
M = ncol(train_R_nm) #number of items

R_nm = train_R_nm

#load the training cluster matrix
load("../Data Structures/train_I_nm.rda")

I_nm = train_I_nm

#test how posterior means do on prediction on training set
post_R_nm = matrix(0,nrow = N, ncol = M)
for(n in 1:N){
  k = clust_id[n]
  id = clust_map[n,3]
  for(m in user_rated_k_n[[k]][[id]]){
    post_R_nm[n,m] = sum(post_U_k_nd[[k]][id,] * post_V_k_md[[k]][m,]) + mean_rating_k[[k]]
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
  k = clust_id[n]
  id = clust_map[n,3]
  for(m in user_rated_n[[n]]){
    post_R_nm[n,m] = sum(post_U_k_nd[[k]][id,] * post_V_k_md[[k]][m,]) + mean_rating_k[[k]]
  }
}
R_nm = R_nm[which(I_nm>0)]
post_R_nm = post_R_nm[which(I_nm>0)]
plot(R_nm[], post_R_nm[])

MAE = mean(abs(post_R_nm - R_nm))
MAE #this is 0.9901067 with 20 iterations, lol
#but apparently 45.04 with 2000 warmup and 1000 iterations???
