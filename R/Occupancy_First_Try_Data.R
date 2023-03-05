#this script aims to implement an occupancy-type model for clustering users/movies 
#and predicting which users watched and rated which movies.
seed = 1
set.seed(seed)

#######################
# Generate true model #
#######################

#set parameters
K = 2 #number of user clusters
L = 2 #number of movie clusters

#hyperparamters for prior distirubtions
a = 1
b = 1
c = 1 
d = 1
alpha = 1
beta = 1

load("cluster_mat.rda")

I_nm = cluster_mat
rm(cluster_mat)
N = nrow(I_nm)
M = ncol(I_nm)

##############
# Initialize #
##############

#initialize phi randomly
phi_k = rep(0,K)
for(k in 1:K){
  phi_k[k] = rgamma(1,1,1)
}
phi_k = phi_k/sum(phi_k)



#initialize U_n randomly
U_n = sample(1:K,N,replace = TRUE)



#initailize theta_l at random
theta_l = rep(0,L)
for(l in 1:L){
  theta_l[l]= rgamma(1,1,1)
}
theta_l = theta_l/sum(theta_l)


#initialize M_m at random
M_m = sample(1:L,M,replace= TRUE)



#intiailize psi_kl randomly
psi_kl = matrix(0,nrow = K,ncol = L)
for(k in 1:K){
  for(l in 1:L){
    psi_kl[k,l] = rbeta(1,1,1)
  }
}

#initialzie W_nm at random
W_nm = matrix(0,nrow = N,ncol = M)
for(n in 1:N){
  for(m in 1:M){
    W_nm[n,m] = sample(0:1,1)
  }
}

#intiailize p_kl randomly
p_kl = matrix(0,nrow = K,ncol = L)
for(k in 1:K){
  for(l in 1:L){
    p_kl[k,l] = rbeta(1,1,1)
  }
}


#find counts

#find the number of users in each cluster
num_u_k = rep(0,K)
for(k in 1:K){
  num_u_k[k] = length(which(U_n==k))
}

#find the number of movies in each clsuter
num_m_l = rep(0,L)
for(l in 1:L){
  num_m_l[l] = length(which(M_m==l))
}

#find the number of each user-movie combo
num_tot_kl = matrix(0,nrow = K,ncol = L)
num_W_kl = matrix(0,nrow = K,ncol = L)
num_WI_kl = matrix(0,nrow = K,ncol = L)
for(n in 1:N){
  for(m in 1:M){
    num_tot_kl[U_n[n],M_m[m]] = num_tot_kl[U_n[n],M_m[m]] + 1
    num_W_kl[U_n[n],M_m[m]] = num_W_kl[U_n[n],M_m[m]] + W_nm[n,m]
    num_WI_kl[U_n[n],M_m[m]] = num_WI_kl[U_n[n],M_m[m]] + W_nm[n,m]*I_nm[n,m]
  }
}


#set number of iterations
warmup = 3000
iterations = 1000
n_it = warmup + iterations

#set up data structures for chain
chain_phi_ik = matrix(0,nrow= iterations,ncol = K)
chain_u_in = matrix(0,nrow = iterations, ncol = N)
chain_theta_il = matrix(0,nrow = iterations, ncol = L)
chain_m_im = matrix(0,nrow = iterations,ncol = M)
chain_psi_ikl = array(0,dim=c(iterations,K,L))
# chain_W_inm = array(0,dim=c(iterations,N,M))
chain_p_ikl = array(0,dim = c(iterations,K,L))

#################
# Gibbs Sampler #
#################

for(it in 1:n_it){
  
  #######
  # Phi #
  #######
  for(k in 1:K){
    phi_k[k] = rgamma(1,alpha+num_u_k[k],1)
  }
  phi_k = phi_k/sum(phi_k)
  
  #######
  # U_n #
  #######
  for(n in 1:N){
    #record currentassignmnet
    old_k = U_n[n]
    
    #de-incremnt count
    num_u_k[old_k] = num_u_k[old_k] - 1
    
    #make a probablitiy vector for full conditoinal
    #on log-scale
    probs = rep(0,K)
    for(k in 1:K){
      probs[k] = sum(W_nm[n,]*log(psi_kl[k,M_m]) + (1-W_nm[n,])*log(1-psi_kl[k,M_m])) + phi_k[k]
    }
    
    #normalize probs
    probs = exp(probs - max(probs))/sum(exp(probs-max(probs)))
    
    #draw from full conditional
    new_k = sample(1:K,1,replace =TRUE,prob = probs)
    
    #record new value
    U_n[n] = new_k
    
    #re-increment count
    num_u_k[new_k] = num_u_k[new_k] + 1
  }
  
  #########
  # theta #
  #########
  for(l in 1:L){
    theta_l[l]= rgamma(1,beta + num_m_l[l],1)
  }
  theta_l = theta_l/sum(theta_l)

  #######
  # M_m #
  #######
  for(m in 1:M){
    #record currentassignmnet
    old_l = M_m[m]

    #de-incremnt count
    num_m_l[old_l] = num_m_l[old_l] - 1

    #make a probablitiy vector for full conditoinal
    #on log-scale
    probs = rep(0,L)
    for(l in 1:L){
      probs[l] = sum(W_nm[,m]*log(psi_kl[U_n,l]) + (1-W_nm[,m])*log(1-psi_kl[U_n,l])) + theta_l[l]
    }

    #normalize probs
    probs = exp(probs - max(probs))/sum(exp(probs-max(probs)))

    #draw from full conditional
    new_l = sample(1:L,1,replace =TRUE,prob = probs)

    #record new value
    M_m[m] = new_l

    #re-increment count
    num_m_l[new_l] = num_m_l[new_l] + 1
  }
  
  
  ########
  # W_nm #
  ########
  for(n in 1:N){
    for(m in 1:M){
      if(I_nm[n,m]==1){
        W_nm[n,m]=1
      } else{
        k = U_n[n]
        l = M_m[m]
        probs = c(1-psi_kl[k,l],(1-p_kl[k,l])*psi_kl[k,l])
        probs = probs/sum(probs)
        W_nm[n,m] = sample(0:1,1,replace= FALSE,prob = probs)
      }
    }
  }
  
  
  ##########
  # psi_kl #
  ##########
  #recalculate counts (maybe this could be faster?)
  num_tot_kl = matrix(0,nrow = K,ncol = L)
  num_W_kl = matrix(0,nrow = K,ncol = L)
  num_WI_kl = matrix(0,nrow = K,ncol = L)
  for(n in 1:N){
    for(m in 1:M){
      num_tot_kl[U_n[n],M_m[m]] = num_tot_kl[U_n[n],M_m[m]] + 1
      num_W_kl[U_n[n],M_m[m]] = num_W_kl[U_n[n],M_m[m]] + W_nm[n,m]
      num_WI_kl[U_n[n],M_m[m]] = num_WI_kl[U_n[n],M_m[m]] + W_nm[n,m]*I_nm[n,m]
    }
  }
  #sample from full conditional
  for(k in 1:K){
    for(l in 1:L){
      psi_kl[k,l] = rbeta(1,a + num_W_kl[k,l],b + num_tot_kl[k,l] - num_W_kl[k,l])
    }
  }
  
  ########
  # p_kl #
  ########
  for(k in 1:K){
    for(l in 1:L){
      p_kl[k,l] = rbeta(1,c + num_WI_kl[k,l],d + num_W_kl[k,l] - num_WI_kl[k,l])
    }
  }
  
  
  ##########
  # Record #
  ##########
  if(it>warmup){
    #record phi
    chain_phi_ik[it-warmup,] = phi_k
    
    #record U_n
    chain_u_in[it-warmup,] = U_n
    
    #record theta
    chain_theta_il[it-warmup,] = theta_l
    
    #record M_m
    chain_m_im[it-warmup,] = M_m
    
    #record psi_kl
    chain_psi_ikl[it-warmup,,] = psi_kl
    
    # #record W_nm
    # chain_W_inm[it-warmup,,] = W_nm
    
    #record p_kl
    chain_p_ikl[it-warmup,,] = p_kl
  }
  
  print(it)
}

######################
# posterior analysis #
######################

#phi
post_phi_k = apply(chain_phi_ik[,],2,mean)

#U_n
post_U_n = apply(chain_u_in[,],2,mean)

#theta
post_theta_l = apply(chain_theta_il,2,mean)

#M_m
post_M_m = apply(chain_m_im,2,mean)

#psi_kl
post_psi_kl = matrix(0,nrow = K,ncol = L)
for(k in 1:K){
  for(l in 1:L){
    post_psi_kl[k,l] = mean(chain_psi_ikl[,k,l])
  }
}


#W_nm
post_W_nm = matrix(0,nrow = N,ncol =M)
# for(n in 1:N){
#   for(m in 1:M){
#     post_W_nm[n,m] = mean(chain_W_inm[,n,m])
#   }
# }

#p_kl
post_p_kl = matrix(0,nrow = K,ncol = L)
for(k in 1:K){
  for(l in 1:L){
    post_p_kl[k,l] = mean(chain_p_ikl[,k,l])
  }
}


Occ_out = list("post_phi" = post_phi_k,
               "post_U" = post_U_n,
               "post_theta" = post_theta_l,
               "post_M" = post_M_m,
               "post_psi" = post_psi_kl,
               "post_W" = post_W_nm,
               "post_p" = post_p_kl)

save(Occ_out,file = "Occ_out_data.rda")