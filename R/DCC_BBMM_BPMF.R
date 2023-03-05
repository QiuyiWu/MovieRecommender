#This script implements a gibbs sampler for the beta-binomial clustering method
#generate synthetic data from the model, run the model, see if it works

library(stats)

set.seed(12)

#################
# Generate Data #
#################

#Set the number of clusters
K = 3
#set D
D = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#set hyperparameters
mu_u0 = rep(0,D)
mu_v0 = rep(0,D)
v_0 = D
W_0 = diag(D)
W_0_inv = solve(W_0)
beta_0 = 2
alpha = 2
a = 1
b = 1

#load I_nm
#load training set
load("train_I_nm.rda")
I_nm = train_I_nm
rm(train_I_nm)

N = nrow(I_nm)
M = ncol(I_nm)

#generate a list of which items each user has rated
user_rated_n = NULL
for(n in 1:N){
  user_rated_n[[n]] = which(I_nm[n,]>0)
}

#generate a list of which users have rated each item
item_rated_m = NULL
for(m in 1:M){
  item_rated_m[[m]] = which(I_nm[,m]>0)
}

#generate the BPMF side of the model


#################
# Generate R_nm #
#################

#load training set
load("train_R_nm.rda")
R_nm = train_R_nm
rm(train_R_nm)

#################
# Gibbs sampler #
#################


#initialize phi at random
phi_k = rgamma(K,1,1)
phi_k = phi_k/sum(phi_k)


#initialize z randomly
z_n = sample(1:K,N,replace = TRUE)

#create a map from global ID to cluster ID
clust_map = matrix(0,nrow = N,ncol = 3)
clust_map[,1] = 1:N
clust_map[,2] = z_n
#compute cluster ID
for(k in 1:K){
  clust_map[which(z_n == k),3] = 1:length(which(z_n==k))
}

#find the number of users assigned to each cluster
num_clust_k = rep(0,K)
for(k in 1:K){
  num_clust_k[k] = length(which(z_n == k))
}

#initialize p randomly
p_km = matrix(0,nrow = K, ncol = M)
for(k in 1:K){
  for(m in 1:M){
    p_km[k,m] = rbeta(1,a,b)
  }
}

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
U_knd = array(0,dim=c(K,N,D))
for(k in 1:K){
  for(n in 1:N){
    U_knd[k,n,] = rnorm(D)
  }
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
V_kmd = array(0,dim=c(K,M,D))
for(k in 1:K){
  for(m in 1:M){
    V_kmd[k,m,] = rnorm(D)
  }
}

#set number of iterations
warmup = 3000
iterations = 1000
n_it = warmup + iterations

#make chains
chain_phi_ik = matrix(0,nrow = iterations,ncol = K)
chain_z_in = matrix(0,nrow = iterations, ncol = N)
chain_p_ikm = array(0,dim = c(iterations,K,M))

#create data structures to store chains
chain_mu_u_kid = array(0,dim= c(K,iterations,D))

chain_U_k_ind = NULL
chain_U_k_ind[[K+1]] = array(0,dim = c(iterations, N, D))
for(k in 1:K){
  chain_U_k_ind[[k]] = array(0,dim = c(iterations, N, D))
}
chain_U_k_ind[[K+1]] = NULL

chain_mu_v_kid = array(0,dim= c(K,iterations,D))

chain_V_k_imd = NULL
chain_V_k_imd[[K+1]] = array(0,dim = c(iterations, M, D))
for(k in 1:K){
  chain_V_k_imd[[k]] = array(0,dim = c(iterations, M, D))
}
chain_V_k_imd[[K+1]] = NULL

#demean data
mean_rating = mean(R_nm[I_nm>0])
R_nm[I_nm>0] = R_nm[I_nm>0] - mean_rating


for(it in 1:n_it){
  
  #######
  # phi #
  #######
  
  #sample phi from its full conditional
  for(k in 1:K){
    phi_k[k] = rgamma(1,num_clust_k[k]+alpha,1)
  }
  phi_k = phi_k/sum(phi_k)
  
  #######
  # z_n #
  #######
  
  #for each user
  for(n in 1:N){
    
    #de-increment num_cluster
    num_clust_k[z_n[n]] = num_clust_k[z_n[n]] - 1
    
    #define a probability vector for the full conditional
    probs = rep(0,K)
    
    for(k in 1:K){ #for each cluster
      probs[k] = sum(I_nm[n,]*log(p_km[k,])) + sum((1-I_nm[n,])*log(1-p_km[k,])) + log(phi_k[k])
      for(m in user_rated_n[[n]]){
        probs[k] = probs[k] - 1/2*((R_nm[n,m] - sum(U_knd[k,n,]*V_kmd[k,m,]))/alpha)^2
      }
    }
    
    #convert from log-scale to normal scale
    probs = exp(probs - max(probs))/sum(exp(probs - max(probs)))
    
    #correct for numerical instability?
    probs[probs < 0.00001] = 0.00001
    probs[probs > 0.99999] = 0.99999
    
    #sample from the mulitinomal full conditional
    z_n[n] = sample(1:K,1,prob = probs)
    
    #re-increment num_cluster
    num_clust_k[z_n[n]] = num_clust_k[z_n[n]] + 1
    
  }
  
  clust_map[,2] = z_n
  #compute cluster ID
  for(k in 1:K){
    clust_map[which(z_n == k),3] = 1:length(which(z_n==k))
  }
  #####
  # p #
  #####
  
  #sample p from its full conditional
  for(k in 1:K){ #for each topic
    for(m in 1:M){ #for each item
      users = which(z_n == k)
      u_count = num_clust_k[k]
      y_count = sum(I_nm[users,m])
      p_km[k,m] = rbeta(1,a + y_count,b + u_count - y_count)
    }
  }
  
  #cycle through clusters
  for(k in 1:K){
    
    ###############################
    # set params based on cluster #
    ###############################
    
    #users
    users = which(z_n ==k)
    if(length(users)>0){
      clust_N = length(users)
      
      #U params
      Lambda_U = Lambda_U_kdd[k,,]
      Lambda_U_inv = Lambda_U_inv_kdd[k,,]
      mu_u = mu_u_kd[k,]
      U = U_knd[k,,]
      
      #V params
      Lambda_V = Lambda_V_kdd[k,,]
      Lambda_V_inv = Lambda_V_inv_kdd[k,,]
      mu_v = mu_v_kd[k,]
      V = V_kmd[k,,]
      
      #which users/items have ratings
      clust_user_rated_n =NULL
      for(user in 1:length(users)){
        clust_user_rated_n[[user]] = user_rated_n[[users[user]]]
      }
      
      clust_item_rated_m = NULL
      for(mov in 1:M){
        temp = clust_map[item_rated_m[[mov]],]
        if(is.null(nrow(temp)) == FALSE){
          clust_item_rated_m[[mov]] = temp[temp[,2]==k,1]
        } else {
          if(temp[2] ==k){
            clust_item_rated_m[[mov]] = temp[1]
          } else {
            clust_item_rated_m[[mov]] = integer(0) #this is janky AF 
          }
        }
      }
      
      #####################
      # mu_u and Lambda_u #
      #####################
      
      #find sufficient stats first
      U_mean = apply(U[users,],2,mean)
      S = matrix(0,nrow = ncol(U),ncol = ncol(U))
      for(n in users){
        S = S + U[n,] %*% t(U[n,])
      }
      S = S/clust_N
      
      #find full conditional parameters for mu_u and Lambda_u
      beta_star = beta_0 + clust_N
      v_star = v_0 + clust_N
      avg = (beta_0*mu_u0  + clust_N*U_mean)/(beta_star)
      W_0_star_inv = W_0_inv + clust_N*S + beta_0*clust_N/beta_star * (U_mean - mu_u0) %*% t((U_mean - mu_u0))
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
          s = s + V[m,]*(R_nm[n,m])
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
      W_0_star_inv = W_0_inv + M*S + beta_0*M/beta_star * (mu_v0 - V_mean) %*% t((mu_v0 - V_mean))
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
        for(n in clust_item_rated_m[[m]]){
          S = S + U[n,] %*% t(U[n,])
        }
        
        s = rep(0,D)
        for(n in clust_item_rated_m[[m]]){
          s = s + U[n,]*(R_nm[n,m])
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
      U_knd[k,,] = U
      
      #V params
      Lambda_V_kdd[k,,] = Lambda_V
      Lambda_V_inv_kdd[k,,] = Lambda_V_inv
      mu_v_kd[k,] = mu_v
      V_kmd[k,,] = V
      
      
    }
    
    
  }#end cluster loop
  
  #record the sampled values
  if(it>warmup){
    #record phi
    chain_phi_ik[it-warmup,] = phi_k
    
    #record z
    chain_z_in[it-warmup,] = z_n
    
    #record p
    chain_p_ikm[it-warmup,,] = p_km
    
    #record mu_u
    chain_mu_u_kid[,it-warmup,] = mu_u_kd[,]
    
    #record mu_v
    chain_mu_v_kid[,it-warmup,] = mu_v_kd[,]
    
    
    for(k in 1:K){
      #record U
      chain_U_k_ind[[k]][it-warmup,,] = U_knd[k,,]
      
      #record V
      chain_V_k_imd[[k]][it-warmup,,] = V_kmd[k,,]
    }
  }
  
  print(it)
}

######################
# Posterior Analysis #
######################

#posterior mean for phi
post_phi_k = apply(chain_phi_ik,2,mean)

#z_n
post_z_nk = matrix(0,nrow = N,ncol = K)
for(n in 1:N){
  for(k in 1:K){
    post_z_nk[n,k] = length(which(chain_z_in[,n]==k))/length(chain_z_in[,n])
  }
}

#posterior mean for p
post_p_km = matrix(0,nrow = K, ncol = M)
for(k in 1:K){
  for(m in 1:M){
    post_p_km[k,m] = mean(chain_p_ikm[,k,m])
  }
}

#posterior mean for mu_u
post_mu_u_kd = matrix(0,nrow = K, ncol = D)
for(k in 1:K){
  post_mu_u_kd[k,] = apply(chain_mu_u_kid[k,,],2,mean)
}

#posterior mean for U
post_U_k_nd = NULL
post_U_k_nd[[K+1]] = matrix(0,nrow = N,ncol = D)
for(k in 1:K){
  post_U_k_nd[[k]] = matrix(0,nrow = N,ncol = D)
  for(n in 1:N){
    post_U_k_nd[[k]][n,] = apply(chain_U_k_ind[[k]][,n,],2,mean)
  }
}
post_U_k_nd[[K+1]] = NULL

#posterior mean for mu_v
post_mu_v_kd = matrix(0,nrow = K, ncol = D)
for(k in 1:K){
  post_mu_v_kd[k,] = apply(chain_mu_v_kid[k,,],2,mean)
}

#posterior mean for V
post_V_k_md = NULL
post_V_k_md[[K+1]] = matrix(0,nrow = N,ncol = D)
for(k in 1:K){
  post_V_k_md[[k]] = matrix(0,nrow = M,ncol = D)
  for(m in 1:M){
    post_V_k_md[[k]][m,] = apply(chain_V_k_imd[[k]][,m,],2,mean)
  }
}
post_V_k_md[[K+1]] = NULL

#test how posterior means do on prediction on training set
post_R_nm = matrix(0,nrow = N, ncol = M)
for(n in 1:N){
  for(m in user_rated_n[[n]]){
    for(k in 1:K){
      post_R_nm[n,m] = post_R_nm[n,m] + post_z_nk[n,k]*sum(post_U_k_nd[[k]][n,] * post_V_k_md[[k]][m,])
    }
    post_R_nm[n,m] = post_R_nm[n,m] + mean_rating
  }
}
R_nm = R_nm[which(I_nm>0)]
post_R_nm = post_R_nm[which(I_nm>0)]
plot(R_nm[], post_R_nm[])

#######################
# Test Set Prediction #
#######################

load("test_R_nm.rda")
#rename for cooperation with old code
R_nm = test_R_nm

#remove old name
rm(test_R_nm)

#load cluster_mat
load("test_I_nm.rda")

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
    for(k in 1:K){
      post_R_nm[n,m] = post_R_nm[n,m] + post_z_nk[n,k]*sum(post_U_k_nd[[k]][n,] * post_V_k_md[[k]][m,])
    }
    post_R_nm[n,m] = post_R_nm[n,m] + mean_rating
  }
}
R_nm = R_nm[which(I_nm>0)]
post_R_nm = post_R_nm[which(I_nm>0)]
plot(R_nm[], post_R_nm[])

MAE = mean(abs(post_R_nm - R_nm))
MAE #this is 0.9901067 with 20 iterations, lol
#but apparently 45.04 with 2000 warmup and 1000 iterations???
rMSE = sqrt(mean((post_R_nm - R_nm)^2))
save(MAE,file = paste("MAE_D_",toString(D),".rda",sep=""))
save(rMSE,file = paste("rMSE_D_",toString(D),".rda",sep=""))
save(post_mu_u_kd,file=paste("post_mu_u_D_",toString(D),".rda",sep=""))
save(post_mu_v_kd,file=paste("post_mu_v_D_",toString(D),".rda",sep=""))
save(post_z_nk,file=paste("post_z_",toString(D),".rda",sep=""))