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
N = 1000
M = 1000 #the number of movies
p = 2 #the number of covariates
Pi = diag(p)
Lambda = diag(p)
sigma2_alpha = 1
sigma2_delta = 1
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

#generate true pi
true_pi = chol(Pi) %*% rnorm(p,0,1)

#generate true alpha
true_delta = sqrt(sigma2_delta)*rnorm(1,0,1)

#generate covariates
z_mp = matrix(0,nrow = M, ncol = p)
for(m in 1:M){
  z_mp[m,] = chol(Pi) %*% rnorm(p,0,sigma)
}

#generate true_theta_m
true_theta_m = rep(0,M)
for(m in 1:M){
  true_theta_m[m] = true_delta + z_mp[m,] %*% true_pi
}

I_nm = matrix(0,nrow = N, ncol = M)
for(m in 1:M){
  prob_s = exp(true_theta_m[m])/(1 + exp(true_theta_m[m]))
  prob_f = 1/(1 + exp(true_theta_m[m]))
  probs = c(prob_f,prob_s)
  for(n in 1:N){
    if(true_W_nm[n,m] > 0){
      I_nm[n,m] = sample(0:1,size = 1,prob = probs)
    }
  }
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
kappa_nm = W_nm - 1/2

#initalize lambda
lambda_nm = I_nm - W_nm/2


#initialize delta at truth
delta = true_delta
# #initialize delta at random
# delta = rnorm(1,0,1)

#initialize pi at truth
pi = true_pi
# #initialize pi at random
# pi = chol(Pi)%*%rnorm(p,0,1)

#initialize polya-gamma variables
v_nm = matrix(0,nrow = N, ncol = M)
for(m in 1:M){
  theta = delta + z_mp[m,]%*%pi
  for(n in 1:N){
    v_m[m] = pgdraw(1,theta)
  }
}

#set number of iterations
warmup = 1000
iterations = 1000
n_it = warmup + iterations

#set up chains
chain_alpha_i = rep(0,iterations)
chain_gamma_ip = matrix(0,nrow = iterations, ncol = p)
# chain_w_inm = array(0,dim = c(iterations,N,M))
chain_delta_i = rep(0,iterations)
chain_pi_ip = matrix(0,nrow = iterations,ncol = p)


#################
# Gibbs Sampler #
#################

for(it in 1:n_it){
  
  # #######
  # # w_m #
  # #######
  # for(m in 1:M){
  #   psi = alpha + y_mp[m,]%*%gamma
  #   for(n in 1:N){
  #     w_nm[n,m] = pgdraw(1,psi)
  #   }
  # }
  # 
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
  # 
  # #########
  # # gamma #
  # #########
  # #find V_gamma
  # V_gamma = matrix(0,nrow = p, ncol = p)
  # for(n in 1:N){
  #   V_gamma = V_gamma + t(y_mp) %*% diag(w_nm[n,]) %*% y_mp
  # }
  # V_gamma = V_gamma + solve(Lambda)
  # V_gamma = solve(V_gamma)
  # 
  # #find m_gamma
  # m_gamma = rep(0,p)
  # for(n in 1:N){
  #   m_gamma = m_gamma + t(y_mp) %*% kappa_nm[n,] - t(y_mp) %*% diag(w_nm[n,]) %*% rep(alpha,M)
  # }
  # m_gamma = V_gamma %*% m_gamma
  # gamma = chol(V_gamma)%*%rnorm(p,0,1) + m_gamma
  
  ########
  # W_nm #
  ########
  for(m in 1:M){
    prob_s = 1/(1+exp(delta+z_mp[m,]%*%pi))*exp(alpha+y_mp[m,]%*%gamma)/(1 + exp(alpha+y_mp[m,]%*%gamma))
    prob_f = 1/(1+exp(delta+z_mp[m,]%*%pi))*1/(1 + exp(alpha+y_mp[m,]%*%gamma))
    probs = c(prob_f,prob_s)
    probs = probs/sum(probs)
    for(n in 1:N){
      if(I_nm[n,m] > 0){
        W_nm[n,m] = 1
      } else {
        W_nm[n,m] = sample(0:1,size = 1,prob = probs)
      }
    }
  }
  #update lambda_nm
  lambda_nm = I_nm - W_nm/2

  
  ########
  # v_nm #
  ########
  for(m in 1:M){
    theta = delta + z_mp[m,]%*%pi
    for(n in 1:N){
      if(W_nm[n,m]>0){
        v_nm[n,m] = pgdraw(1,theta)
      } else {
        v_nm[n,m] = 0
      }
    }
  }

  
  #########
  # delta #
  #########
  V_delta = 1/(sum(v_nm[W_nm>0]) + 1/sigma2_delta)
  m_delta = 0
  for(n in 1:N){
    for(m in 1:M){
      if(W_nm[n,m]>0){
        m_delta = m_delta + lambda_nm[n,m] - v_nm[n,m]*(z_mp[m,]%*%pi)
      }
    }
  }
  m_delta = V_delta*m_delta
  delta = sqrt(V_delta)*rnorm(1,0,1) + m_delta

  ######
  # pi #
  ######
  #Find v_pi
  V_pi = matrix(0,nrow = p, ncol = p)
  for(n in 1:N){
    ind = which(W_nm[n,]>0)
    V_pi = V_pi + t(z_mp[ind,]) %*% diag(v_nm[n,ind]) %*% z_mp[ind,]
  }
  V_pi = V_pi + solve(Pi)
  V_pi = solve(V_pi)

  #find m_pi
  m_pi = rep(0,p)
  for(n in 1:N){
    ind = which(W_nm[n,]>0)
    m_pi = m_pi + t(z_mp[ind,]) %*% lambda_nm[n,ind] - t(z_mp[ind,]) %*% diag(v_nm[n,ind]) %*% rep(delta,length(ind))
  }
  m_pi = V_pi %*% m_pi
  pi = chol(V_pi)%*%rnorm(p,0,1) + m_pi
  
  ##########
  # Record #
  ##########
  if(it > warmup){
    chain_alpha_i[it-warmup] = alpha
    chain_gamma_ip[it-warmup,] = gamma
    # chain_w_inm[it-warmup,,] = W_nm
    chain_delta_i[it-warmup] = delta
    chain_pi_ip[it-warmup,] = pi
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

#W_m
post_W_nm = matrix(0,nrow = N, ncol = M)
for(n in 1:N){
  for(m in 1:M){
    post_W_nm[n,m] = mean(chain_w_inm[,n,m])
  }
}
plot(exp(true_psi_m)/(1+exp(true_psi_m)),post_W_nm[4,])
plot(exp(true_psi_m)/(1+exp(true_psi_m)),apply(post_W_nm[,],2,mean))

#delta
post_delta = mean(chain_delta_i)
post_delta
plot(chain_delta_i,type="l")

#pi
post_pi_p = rep(0,p)
for(i in 1:p){
  post_pi_p[i] = mean(chain_pi_ip[,i])
}
post_pi_p
matplot(chain_pi_ip,type="l")

#theta
post_theta_m = rep(0,M)
for(m in 1:M){
  post_theta_m[m] = post_delta + z_mp[m,]%*%post_pi_p
}
plot(true_theta_m,post_theta_m)
plot(exp(true_theta_m)/(1+exp(true_theta_m)),exp(post_theta_m)/(1+exp(post_theta_m)))
lines(0:1,0:1)