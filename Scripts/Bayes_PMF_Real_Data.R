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
