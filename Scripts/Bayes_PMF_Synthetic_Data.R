#This script generates data from a Bayseian PMF distribution and then implements a Gibbs sampler
library(stats)

#generate data
set.seed(1)

#load cluster_mat
load("../Data Structures/cluster_mat.rda")

#rename to indicator functions
I_nm = cluster_mat[1:100,] #subset for time...

# recover number of users N and number of items M
N= nrow(I_nm)
M = ncol(I_nm)

#remove cluster_mat for memory rasons
rm(cluster_mat)


#set D
D = 5

#set hyperparameters
mu_u0 = mu_v0 = rep(0,D)
v_0 = D
W_0 = diag(D)
W_0_inv = solve(W_0)
beta_0 = 1
alpha = 2


##############
# generate U #
##############

#draw Lambda_U
true_Lambda_U = rWishart(1,v_0,W_0)
true_Lambda_U = true_Lambda_U[,,1]
true_Lambda_U_inv = solve(true_Lambda_U)
#draw mu_u
true_mu_u = chol(1/beta_0*true_Lambda_U_inv) %*% rnorm(D) + mu_u0

#draw U_i for each user u_i
true_U = matrix(0,nrow = N,ncol = D)
for(i in 1:nrow(true_U)){
  true_U[i,] = chol(true_Lambda_U_inv) %*% rnorm(D) + true_mu_u
}

##############
# generate V #
##############

#draw Lambda_V
true_Lambda_V = rWishart(1,v_0,W_0)
true_Lambda_V = true_Lambda_V[,,1]
true_Lambda_V_inv = solve(true_Lambda_V)
#draw mu_v
true_mu_v = chol(1/beta_0*true_Lambda_V_inv) %*% rnorm(D) + mu_u0

#draw V_i for each item m_i
true_V = matrix(0,nrow = M,ncol = D)
for(i in 1:nrow(true_V)){
  true_V[i,] = chol(true_Lambda_V_inv) %*% rnorm(D) + true_mu_v
}

##############
# generate R #
##############

#make R 
R_nm = matrix(0,nrow = N, ncol = M)
count = 0
for(i in 1:N){
  for(j in 1:M){
    if( (I_nm[i,j]>0) == TRUE){
      R_nm[i,j] = 1/alpha * rnorm(1) + sum(true_U[i,] * true_V[j,])
    }
  }
}


##############
# Initialize #
##############

# #initialize Lambda_U at random
# Lambda_U = rWishart(1,v_0,W_0)
# Lambda_U = Lambda_U[,,1]
# Lambda_U_inv = solve(Lambda_U)
# #initialize to truth
Lambda_U = true_Lambda_U
Lambda_U_inv = solve(Lambda_U)

# #initiailze mu_u at random
# mu_u = rnorm(D)
# #initialize to truth
mu_u = true_mu_u

# #initialize U at random
# U = matrix(0,nrow = N,ncol = D)
# for(i in 1:nrow(U)){
#   U[i,] = rnorm(D) 
# }
#initialize to truth
U = true_U

#initialize V at random
#initailize to truth
V = true_V

#set number of iterations
warmup = 1
iterations = 1000
n_it = warmup + iterations

#create data structures to store chains
chain_mu_u_id = matrix(0,nrow= iterations, ncol = D)
chain_U_ind = array(0,dim=c(iterations,N,D))

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
    S = S + U[i,] %*% t(U[i,])
  }
  S = S/N

  #find full conditional parameters for mu_u and Lambda_u
  beta_star = beta_0 + N
  v_star = v_0 + N
  avg = (beta_0*mu_u0  + N*U_mean)/(beta_star)
  W_0_star_inv = W_0_inv + N*S + beta_0*N/beta_star * (mu_u0 - U_mean) %*% t((mu_u0 - U_mean))
  W_0_star = solve(W_0_star_inv)

  #draw Lambda_U first
  Lambda_U = rWishart(1,v_star,W_0)
  Lambda_U = Lambda_U[,,1]
  Lambda_U_inv = solve(Lambda_U)

  #initiailze mu_u at random
  mu_u = chol(1/beta_star*Lambda_U_inv) %*% rnorm(D) + avg

  #####
  # U #
  #####
  
  for(n in 1:N){
    
    #find sufficient stats
    S = matrix(0,nrow = D, ncol = D)
    for(m in 1:M){
      if((I_nm[n,m]>0) == TRUE){
        S = S + V[m,] %*% t(V[m,])
      }
    }
    
    s = rep(0,D)
    for(m in 1:M){
      if((I_nm[n,m]>0) == TRUE){
        s = s + V[m,]*R_nm[n,m]
      }
    }
    #maybe I should rename them becaue S and s are similar..
    
    #find variance and mean for normal
    covar = Lambda_U + alpha*S
    covar_inv = solve(covar)
    avg = covar_inv%*%(alpha*s + Lambda_U %*% mu_u )
    
    #draw U_n
    U[n,] = chol(covar_inv) %*% rnorm(D) + avg
  }

  ##########
  # Record #
  ##########

  #record the sampled values
  if(it>warmup){
    #record mu_u
    chain_mu_u_id[it - warmup,] = mu_u
    
    #record U
    chain_U_ind[it - warmup,,] = U

  }

  print(it)

}

#posterior analysis

#posterior mean for mu_u
post_mu_u = apply(chain_mu_u_id,2,mean)
#chainplot
matplot(chain_mu_u_id,type="l")
for(d in 1:D){
  abline(h = true_mu_u[d])
}
#posterior mean vs truth
plot(true_mu_u,post_mu_u)

#posterior mean for U
post_U_nd = matrix(0,nrow = N,ncol = D)
for(n in 1:N){
  post_U_nd[n,] = apply(chain_U_ind[,n,],2,mean)
}
#chainplot
matplot(chain_U_ind[,n,],type="l")

plot(true_U[,],post_U_nd[,])

