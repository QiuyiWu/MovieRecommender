#This script attempts to implement a GP to predict ratings for one user
#The GP operates over the tag genome space

#set random seed
set.seed(1)

#load a ratings matrix
load("../Data Structures/rat_mat.rda")
#recover the number of users and items
num_user = nrow(rat_mat)
num_item = ncol(rat_mat)

#divide data into training and test set
train_items = sample(1:num_item,round(num_item/2))
test_items = 1:num_item
test_items = test_items[-train_items]
train_mat = rat_mat[,train_items]
test_mat = rat_mat[,test_items]

#investigate the trianing set
num_train = NULL
for(i in 1:nrow(train_mat)){
  num_train[i] = length(which(train_mat[i,]>0))
}
min(num_train)
# all users have at least 5 items in the training set 
# we may need to adjust this process to guarantee training set size on a per-user basis

#load the distance matrix
load("../Data Structures/dist_mat.rda")

#set a user 
user = 1

#find the movies which the user has rated in the training and test sets
test_items = which(test_mat[user,]>0)
train_items = which(train_mat[user,]>0)
all_items = c(test_items,train_items)

#subset ratings
test_response = test_mat[user,test_items]
train_response = train_mat[user,train_items]

#subset distance matrix
test_dist = dist_mat[test_items,test_items]
train_dist = dist_mat[train_items,train_items]
all_dist = dist_mat[all_items,all_items]

rm(dist_mat)

#define hyperparameter values
l = 0.01*1:100
tau = 0.001*1:100
sigma = 0.01*90:110

#construct a matrix where each combination is one row
hyper = matrix(0,nrow = length(l)*length(tau)*length(sigma),ncol = 3)
for(i in 1:length(l)){
  for(j in 1:length(tau)){
    for(k in 1:length(sigma)){
      ind = length(tau)*length(sigma)*(i-1) + length(sigma)*(j-1)+k
      hyper[ind,1] = l[i]
      hyper[ind,2] = tau[j]
      hyper[ind,3] = sigma[k]
    }
  }
}
colnames(hyper) = c("l","tau","sigma")

#find average on train response
avg = mean(train_response)

#do a grid search over hyperparameters
ll = rep(0,nrow(hyper)) #pre-allocate memory to record ll
# MAE = rep(0,nrow(hyper)) #pre-allocate memory to record MAE
for(i in 1:nrow(hyper)){
  #set hyperparameters
  l = hyper[i,1]
  tau = hyper[i,2]
  sigma = hyper[i,3]
  
  #find covariance kernels
  # test_kernel = sigma*exp(-test_dist^2/(2*l)) + tau*diag(length(test_items))
  train_kernel = sigma*exp(-train_dist^2/(2*l)) + tau*diag(length(train_items))
  # all_kernel = sigma*exp(-all_dist^2/(2*l)) + tau*diag(length(all_items))

  #find log-likelhood
  ll[i] = -0.5*det(train_kernel) - 0.5*(train_response - avg) %*% solve(train_kernel) %*%(train_response - avg)
  
  # #find mean on training set
  # avg = mean(train_response)
  # 
  # #find predictive conditional distribution on test set
  # A = test_kernel
  # B = all_kernel[1:length(test_items),(length(test_items)+1):ncol(all_kernel)]
  # C = train_kernel
  # 
  # pred_mean = avg + B %*% solve(C) %*%(train_response - avg)
  # 
  # MAE[i] = mean(abs(test_response-pred_mean))
  
  print(i/nrow(hyper))
}

# min(MAE)
# hyper[which(MAE == min(MAE)),]
# mean(abs(avg-test_response))

#set kernel hyperparaemters
hypers = hyper[which(ll == max(ll)),]
l = hypers[1]
tau = hypers[2]
sigma = hypers[3]

#find kernels
test_kernel = sigma*exp(-test_dist^2/(2*l)) + tau*diag(length(test_items))
train_kernel = sigma*exp(-train_dist^2/(2*l)) + tau*diag(length(train_items))
all_kernel = sigma*exp(-all_dist^2/(2*l)) + tau*diag(length(all_items))

#find predictive conditional distribution on test set
A = test_kernel
B = all_kernel[1:length(test_items),(length(test_items)+1):ncol(all_kernel)]
C = train_kernel

#find predicted values
pred_mean = avg + B %*% solve(C) %*%(train_response - avg)

#plot predicted vs actual
plot(test_response,pred_mean)
hypers

mean(abs(test_response-pred_mean))
mean(abs(test_response-avg))

