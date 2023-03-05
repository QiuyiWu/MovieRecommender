#this script aims to implement an occupancy-type model for clustering users/movies 
#and predicting which users watched and rated which movies.
seed = 1
set.seed(seed)

#######################
# Generate true model #
#######################

#set parameters
N = 1000 #number of users
M = 10000 #number o movies
K = 2 #number of user clusters
L = 2 #number of movie clusters

#hyperparamters for prior distirubtions
a = 1
b = 1
c = 1 
d = 1
alpha = 1
beta = 1

#Generate phi
true_phi_k = rep(0,K)
for(k in 1:K){
  true_phi_k[k] = rgamma(1,alpha,1)
}
true_phi_k = true_phi_k/sum(true_phi_k)
#set manually...
true_phi_k[1] = 1/3
true_phi_k[2] = 2/3

#Generate U_n
true_U_n = sample(1:K,N,replace = TRUE,true_phi_k)

#Generate theta
true_theta_l = rep(0,L)
for(l in 1:L){
  true_theta_l[l]= rgamma(1,beta,1)
}
true_theta_l = true_theta_l/sum(true_theta_l)
#set manually
true_theta_l[1] = 1/4
true_theta_l[2] = 3/4

#Generate M_m
true_M_m = sample(1:L,M,replace = TRUE,true_theta_l)

#genreate kappa_kl
true_kappa_kl = matrix(0,nrow = K,ncol = L)
for(k in 1:K){
  for(l in 1:L){
    true_kappa_kl[k,l] = rbeta(1,a,b)
  }
}

#maybe set these manually?
true_kappa_kl[1,1] = 1/5
true_kappa_kl[1,2] = 2/5
true_kappa_kl[2,1] = 3/5
true_kappa_kl[2,2] = 4/5

true_kappa_kl = (0.01412552/mean(true_kappa_kl)*true_kappa_kl)

#Generate W_nm
true_W_nm = matrix(0,nrow = N, ncol = M)
for(n in 1:N){
  for(m in 1:M){
    prob = true_kappa_kl[true_U_n[n],true_M_m[m]]
    true_W_nm[n,m] = sample(0:1,1,replace = TRUE,c(1-prob,prob))
  }
}

###############################
# Set up variational families #
###############################

#initialize gamma_k at a uniform distribution?
gamma_k = rgamma(K,N,K)
# #initialize gamma_k at truth
# gamma_k = N*true_phi_k

#initialize delta_nk at a uniform distribution?
delta_nk = matrix(0,nrow = N,ncol = K)
for(n in 1:N){
  delta_nk[n,] = rgamma(K,1,1)
  delta_nk[n,] = delta_nk[n,]/sum(delta_nk[n,])
}
# #initialize delta_nk at truth
# delta_nk = matrix(0,nrow = N,ncol = K)
# for(n in 1:N){
#   for(k in 1:K){
#     if(true_U_n[n]==k){
#       delta_nk[n,k] = 1
#     }
#   }
# }

#initialize epsilon_l at a uniform
epsilon_l = rgamma(L,M,L)
# #initialize epsilon_l at truth
# epsilon_l = M*true_theta_l

#initialize zeta_ml at random
zeta_ml = matrix(0,nrow = M,ncol = L)
for(m in 1:M){
  zeta_ml[m,] = rgamma(L,1,1)
  zeta_ml[m,] = zeta_ml[m,]/sum(zeta_ml[m,])
}
# #initialize delta_nk at truth
# zeta_ml = matrix(0,nrow = M,ncol = L)
# for(m in 1:M){
#   for(l in 1:L){
#     if(true_M_m[m]==l){
#       zeta_ml[m,l] = 1
#     }
#   }
# }

#initailize eta_kl and iota_kl at random
eta_kl = matrix(0,nrow = K,ncol = L)
iota_kl = matrix(0,nrow = K,ncol = L)
for(k in 1:K){
  for(l in 1:L){
    eta_kl[k,l] = rgamma(1,N*M,2*K*L)
    iota_kl[k,l] = rgamma(1,N*M,2*K*L)
  }
}
# #initialize eta_kl and iota_kl at truth
# eta_kl = true_kappa_kl
# iota_kl = 1- eta_kl
# eta_kl = N*M/(K*L)*eta_kl
# iota_kl = N*M/(K*L)*iota_kl



#############
# CAVI Algo #
#############
n_it = 10

chain_gamma_ik = matrix(0,nrow = n_it,ncol = K)

for(it in 1:n_it){
  
  #######
  # phi #
  #######
  gamma_k = alpha + apply(delta_nk,2,sum)

  #######
  # U_n #
  #######
  for(n in 1:N){
    for(k in 1:K){

      sum = 0
      for(m in 1:M){
        for(l in 1:L){
          sum1 = true_W_nm[n,m]*(digamma(eta_kl[k,l]) - digamma(eta_kl[k,l]+iota_kl[k,l]))
          sum2 = (1-true_W_nm[n,m])*(digamma(iota_kl[k,l]) - digamma(eta_kl[k,l]+iota_kl[k,l]))
          sum = sum + zeta_ml[m,l]*(sum1 + sum2)
        }
      }


      delta_nk[n,k] = sum + digamma(gamma_k[k]) - digamma(sum(gamma_k))
    }
    
    delta_nk[n,] = exp(delta_nk[n,] - max(delta_nk[n,]))/sum(exp(delta_nk[n,]-max(delta_nk[n,])))
  }

  #########
  # theta #
  #########
  epsilon_l = beta + apply(zeta_ml,2,sum)

  #######
  # M_n #
  #######
  for(m in 1:M){
    for(l in 1:L){

      sum = 0
      for(n in 1:N){
        for(k in 1:K){
          sum1 = true_W_nm[n,m]*(digamma(eta_kl[k,l]) - digamma(eta_kl[k,l]+iota_kl[k,l]))
          sum2 = (1-true_W_nm[n,m])*(digamma(iota_kl[k,l]) - digamma(eta_kl[k,l]+iota_kl[k,l]))
          sum = sum + delta_nk[n,k]*(sum1 + sum2)
        }
      }

      zeta_ml[m,l] = sum + digamma(epsilon_l[l]) - digamma(sum(epsilon_l))
    }
    zeta_ml[m,] = exp(zeta_ml[m,] - max(zeta_ml[m,]))/sum(exp(zeta_ml[m,]-max(zeta_ml[m,])))
  }
  
  ##########
  # psi_kl #
  ##########
  
  for(k in 1:K){
    for(l in 1:L){
      
      eta_sum = 0
      iota_sum = 0
      for(n in 1:N){
        for(m in 1:M){
          eta_sum = eta_sum + true_W_nm[n,m]*delta_nk[n,k]*zeta_ml[m,l]
          iota_sum = iota_sum + (1-true_W_nm[n,m])*delta_nk[n,k]*zeta_ml[m,l]
        }
      }
      
      eta_kl[k,l] = eta_sum + a
      iota_kl[k,l] = iota_sum + b
      
    }
  }
  
  
  chain_gamma_ik[it,] = gamma_k
  print(it)
  
}


#look at how good inference is!
matplot(chain_gamma_ik,type="l")
gamma_k
epsilon_l
eta_kl/(eta_kl + iota_kl)