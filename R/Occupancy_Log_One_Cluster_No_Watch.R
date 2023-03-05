#this script attempts to implement a logistic regressoin using Polya-Gamma variables
#to predict the probability a movie is rated for one user
library(pgdraw)

#set random seed 
seed = 1
set.seed(seed)

#########################
# Generate "true" model #
#########################

#set hyperparameters
N = 100
M = 100 #the number of movies
p = 2 #the number of covariates
Lambda = diag(p)
sigma2_alpha = 1
sigma = 1

#generate true gammma
true_gamma = chol(Lambda) %*% rnorm(p,0,1)

#generate true alpha
true_alpha = sqrt(sigma2_alpha)*rnorm(1,0,1)

#generate covariates
y_mp = matrix(0,nrow = M, ncol = p)
for(m in 1:M){
  y_mp[m,] = chol(Lambda) %*% rnorm(p,0,sigma2_alpha) + 1
}

#generate true psi
true_psi_m = rep(0,M)
for(m in 1:M){
  true_psi_m[m] = true_alpha + y_mp[m,] %*% true_gamma  
  # + rnorm(1,0,sigma)
}

#generate true W_m
true_W_nm = matrix(0,nrow = N, ncol = M)
for(m in 1:M){
  prob_s = exp(true_psi_m[m])/(1 + exp(true_psi_m[m]))
  prob_f = 1/(1 + exp(true_psi_m[m]))
  true_W_nm[,m] = sample(0:1,size = N, replace = TRUE,prob = c(prob_f,prob_s))
}


##############
# Initialize #
##############

#initailize alpha at truth
alpha = true_alpha
# #initailize alaph randomly
# alpha = rnorm(1,0,1)

#initialize gamma at truth
gamma = true_gamma
# #initialize gamma randomly
# gamma = chol(Lambda)%*%rnorm(p,0,1)

#initialize polya-gamma variables
w_nm = matrix(0,nrow = N, ncol = M)
for(n in 1:N){
  for(m in 1:M){
    w_nm[m] = pgdraw(1,alpha + y_mp[m,]%*%gamma)
  }
}


#iitaiilze W_m at truth
W_nm = true_W_nm

#initialize kappa
kappa_nm = matrix(0,nrow = N, ncol = M)
for(n in 1:N){
  for(m in 1:M){
    kappa_nm[n,m] = W_nm[n,m] - 1/2
  }
}


#set number of iterations
warmup = 1000
iterations = 1000
n_it = warmup + iterations

#set up chains
chain_alpha_i = rep(0,iterations)
chain_gamma_ip = matrix(0,nrow = iterations, ncol = p)

#################
# Gibbs Sampler #
#################

for(it in 1:n_it){
  
  #######
  # w_m #
  #######
  for(m in 1:M){
    psi = alpha + y_mp[m,]%*%gamma
    for(n in 1:N){
      w_nm[n,m] = pgdraw(1,psi)
    }
  }
  
  # #########
  # # alpha #
  # #########
  # V_alpha = 1/(sum(w_nm) + 1/sigma2_alpha)
  # m_alpha = 0
  # for(n in 1:N){
  #   for(m in 1:M){
  #     m_alpha = m_alpha + kappa_nm[n,m] - w_nm[n,m]*(y_mp[m,]%*%gamma)
  #   }
  # }
  # m_alpha = V_alpha*m_alpha
  # alpha = sqrt(V_alpha)*rnorm(1,0,1) + m_alpha
  
  #########
  # gamma #
  #########
  #find V_gamma
  V_gamma = matrix(0,nrow = p, ncol = p)
  for(n in 1:N){
    V_gamma = V_gamma + t(y_mp) %*% diag(w_nm[n,]) %*% y_mp
  }
  V_gamma = V_gamma + solve(Lambda)
  V_gamma = solve(V_gamma)
  
  #find m_gamma
  m_gamma = rep(0,p)
  for(n in 1:N){
    m_gamma = m_gamma + t(y_mp) %*% kappa_nm[n,] - t(y_mp) %*% diag(w_nm[n,]) %*% rep(alpha,M)
  }
  m_gamma = V_gamma %*% m_gamma
  gamma = chol(V_gamma)%*%rnorm(p,0,1) + m_gamma
  
  
  if(it > warmup){
    chain_alpha_i[it-warmup] = alpha
    chain_gamma_ip[it-warmup,] = gamma
  }
  print(it)
}

#####################
# Posterior analyis #
#####################

#alpha
post_alpha = mean(chain_alpha_i)
post_alpha
plot(chain_alpha_i,type="l")

#gamma
post_gamma_p = rep(0,p)
for(i in 1:p){
  post_gamma_p[i] = mean(chain_gamma_ip[,i])
}
post_gamma_p

matplot(chain_gamma_ip,type="l")


#psi
post_psi_m = rep(0,M)
for(m in 1:M){
  post_psi_m[m] = post_alpha + y_mp[m,]%*%post_gamma_p
}
plot(true_psi_m,post_psi_m)

lines(-10:10,-10:10)
