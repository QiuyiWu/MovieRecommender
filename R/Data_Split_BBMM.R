
set.seed(1)

#load a ratings matrix
load("../Data Structures/rat_mat.rda")
#recover the number of users and items
N = nrow(rat_mat) #number of users
M = ncol(rat_mat) #number of items

#divide data into training and test set
#for each user, put half of that user's item in the test set?
train_R_nm = matrix(0,nrow = N, ncol =M)
test_R_nm = matrix(0,nrow = N, ncol =M)
items = 1:M
for(n in 1:N){
  train_items = sample(x = items,size = round(length(items)/2), replace = FALSE)
  test_items = 1:length(items)
  test_items = test_items[-train_items]
  train_R_nm[n,train_items] = rat_mat[n,train_items]
  test_R_nm[n,test_items] = rat_mat[n,test_items]
}


#make cluster matrices for the training set
train_I_nm = train_R_nm
train_I_nm[which(train_I_nm>0)] = 1

rm(train_R_nm)
rm(test_R_nm)

#This script implements a gibbs sampler for the beta-binomial clustering method
#generate synthetic data from the model, run the model, see if it works

#######################
# Set Hyperparamaters #
#######################

#Set the number of clusters
K = 2
#set alpha
alpha = 1
#set a
a = 1
#set b
b = 1
#maybe these should be modified away from a uniform distribution?

#############
# load Data #
#############

#rename cluster_mat
y_nm = train_I_nm
# y_nm = cluster_mat
#remove for memory
rm(train_I_nm)
# rm(cluster_mat)


#recover the number of users
N = nrow(y_nm)
#recover th enumber o fitems
M = ncol(y_nm)

#################
# Gibbs Sampler #
#################

#for now, we do the parts individually

#initialize phi at random
phi_k = rgamma(K,1,1)
phi_k = phi_k/sum(phi_k)

#initialize z randomly
z_u = sample(1:K,N,replace = TRUE)

#initialize p randomly
p_km = matrix(0,nrow = K, ncol = M)
for(k in 1:K){
  for(m in 1:M){
    p_km[k,m] = rbeta(1,a,b)
  }
}

#find the number of users assigned to each cluster
num_clust_k = rep(0,K)
for(k in 1:K){
  num_clust_k[k] = length(which(z_u == k))
}

#set number of iterations
warmup = 1000
iterations = 1000
n_it = warmup + iterations

#make data structures to store chains
chain_phi_ik = matrix(0,nrow = iterations,ncol = K)
chain_z_iu = matrix(0,nrow = iterations, ncol = N)
chain_p_ikm = array(0,dim = c(iterations,K,M))

###########
# Sampler #
###########

for(it in 1:n_it){
  
  #######
  # phi #
  #######
  
  #sample phi from its full conditional
  for(k in 1:K){
    phi_k[k] = rgamma(1,num_clust_k[k]+alpha,1)
  }
  phi_k = phi_k/sum(phi_k)
  
  #####
  # z #
  #####
  
  #sample phi from its full conditional
  for(n in 1:N){ #for each user
    
    #de-increment num_cluster
    num_clust_k[z_u[n]] = num_clust_k[z_u[n]] - 1
    
    #define a probability vector for the full conditional
    probs = rep(0,K)
    
    for(k in 1:K){ #for each cluster
      probs[k] = sum(y_nm[n,]*log(p_km[k,])) + sum((1-y_nm[n,])*log(1-p_km[k,])) + log(phi_k[k])
    }
    
    #convert from log-scale to normal scale
    probs = exp(probs - max(probs))/sum(exp(probs - max(probs)))
    
    #correct for numerical instability?
    probs[probs < 0.00001] = 0.00001
    probs[probs > 0.99999] = 0.99999
    
    #sample from the mulitinomal full conditional
    z_u[n] = sample(1:K,1,prob = probs)
    
    #re-increment num_cluster
    num_clust_k[z_u[n]] = num_clust_k[z_u[n]] + 1
  }
  
  #####
  # p #
  #####
  
  #sample p from its full conditional
  for(k in 1:K){ #for each topic
    for(m in 1:M){ #for each item
      users = which(z_u == k)
      u_count = length(users)
      y_count = sum(y_nm[users,m])
      p_km[k,m] = rbeta(1,a + y_count,b + u_count - y_count)
    }
  }
  
  ##########
  # Record #
  ##########
  
  #record the sampled values
  if(it>warmup){
    #record phi
    chain_phi_ik[it-warmup,] = phi_k
    
    #record z
    chain_z_iu[it-warmup,] = z_u
    
    #record p
    chain_p_ikm[it-warmup,,] = p_km
    
  }
  
  #print iterate
  print(it)
  
}

BBMM_out = list(Chain_phi = chain_phi_ik,
                Chain_z = chain_z_iu,
                Chain_p = chain_p_ikm)
save(BBMM_out, file = "../Data Structures/BBMM_out.rda")
