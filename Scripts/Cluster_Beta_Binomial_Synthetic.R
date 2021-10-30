#This script implements a gibbs sampler for the beta-binomial clustering method
#generate synthetic data from the model, run the model, see if it works

set.seed(1)

#################
# Generate Data #
#################

#Set the number of clusters
K = 2
#set the number of users
N = 1000
#Set th enumber o fitems
M = 2000

#set alpha
alpha = 1

#generate the true phi from a Dir(alpha) dist
true_phi_k = rgamma(K,alpha,1)
true_phi_k = true_phi_k/sum(true_phi_k)

#generate the true cluster assignments z from the true phi
true_z_u = sample(1:K, N, replace = TRUE, prob = true_phi_k)

#generate the true probabilities
a = 0.01412
b = 1 - a
true_p_km = matrix(0,nrow = K, ncol = M)
for(k in 1:K){
  for(m in 1:M){
    true_p_km[k,m] = rbeta(1,a,b)
  }
}

#generate the true y_nm
y_nm = matrix(0,nrow = N, ncol =M)
for(n in 1:N){
  cluster = true_z_u[n]
  for(m in 1:M){
    probs = c(1 - true_p_km[cluster,m],true_p_km[cluster,m])
    y_nm[n,m] = sample(0:1,1,prob = probs)
  }
}

#################
# Gibbs Sampler #
#################

#for now, we do the parts individually

#initialize phi at random
phi_k = rgamma(K,1,1)
phi_k = phi_k/sum(phi_k)
#initialize phi at truth
# phi_k = true_phi_k

#initialize z randomly
z_u = sample(1:K,N,replace = TRUE)
# #initialize z at truth
# z_u = true_z_u

#initialize p randomly
p_km = matrix(0,nrow = K, ncol = M)
for(k in 1:K){
  for(m in 1:M){
    p_km[k,m] = rbeta(1,a,b)
  }
}
#initialize p at truth
# p_km = true_p_km

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






######################
# Posterior Analysis #
######################

#posterior mean for phi
post_phi_k = apply(chain_phi_ik,2,mean)
#chain plot
matplot(chain_phi_ik,type="l")
abline(h = true_phi_k[1])
abline(h = true_phi_k[2])
#posterior mean vs truth
plot(true_phi_k,post_phi_k)


#posterior mean for z
post_z_uk = matrix(0,nrow = N,ncol = K)
for(n in 1:N){
  for(k in 1:K){
    post_z_uk[n,k] = length(which(chain_z_iu[,n]==k))/length(chain_z_iu[,n])
  }
}
plot(true_z_u,post_z_uk[,1])

#posterior mean for p
post_p_km = matrix(0,nrow = K, ncol = M)
for(k in 1:K){
  for(m in 1:M){
    post_p_km[k,m] = mean(chain_p_ikm[,k,m])
  }
}
plot(true_p_km[1,],post_p_km[2,])
points(true_p_km[2,],post_p_km[1,])
