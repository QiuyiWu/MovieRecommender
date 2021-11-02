#this script clusters users according to which items they rate
#then it implements memory-based CF according to cluster
library(stats)
set.seed(1)

###########################
# Collaborative Filtering #
###########################

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

#define a similarity function - here cosine similarity
sim = function(a,b){
  return(sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2))))
}

#load the cluster matrix
load("../Data Structures/cluster_mat.rda")

#run the k-means algorithm for k = 2
K = 2
out = kmeans(cluster_mat,centers = K, nstart = 25)

#find the cluster ids
clust_id = out$cluster

#look at cluster numbers 
table(clust_id)

#create a map from global ID to cluster ID
clust_map = matrix(0,nrow = num_user,ncol = 3)
clust_map[,1] = 1:num_user
clust_map[,2] = clust_id
#compute cluster ID
for(k in 1:K){
  clust_map[which(clust_id == k),3] = 1:length(which(clust_id==k))
}

#define K different training matrices depending on cluster ID
train_mats = NULL
train_mats[[K +1]] = matrix(0,nrow = 100, ncol = 100)
for(k in 1:K){
  train_mats[[k]] = rat_mat[which(clust_id==k),]
}
train_mats[[K +1]] = NULL

#define K different similarity matrices depending on cluster ID
sim_mats = NULL
sim_mats[[K +1]] = matrix(0,nrow = 100, ncol = 100)
for(k in 1:K){
  #find similarities between each user
  sim_mat = matrix(0,nrow=nrow(train_mats[[k]]),ncol = nrow(train_mats[[k]]))
  for(i in 1:nrow(sim_mat)){
    for(j in 1:ncol(sim_mat)){
      sim_mat[i,j] = sim(train_mats[[k]][i,],train_mats[[k]][j,])
    }
  }
  
  #saveoutput
  sim_mats[[k]] = sim_mat
}
sim_mats[[K +1]] = NULL

#set diag equal to 0
for(k in 1:K){
  diag(sim_mats[[k]]) = 0
}

#set number of nearest neighbors to find
n = 25

#find k nearest neighbors of each user, in the cluster
nn = matrix(0,nrow = num_user,ncol = n)
for(k in 1:K){
  for(i in 1:nrow(sim_mats[[k]])){
    temp = order(sim_mats[[k]][i,],decreasing = TRUE)
    nn[which(clust_map[,2] == k & clust_map[,3] == i),] = temp[1:n]
  }
}

#convert from cluster-specific ID to global ID
for(i in 1:nrow(nn)){
  k = clust_id[i]
  for(j in 1:ncol(nn)){
    nn[i,j] = which(clust_map[,2] == k & clust_map[,3] == nn[i,j])
  }
}

#now we proceed with generic Memory-based Collaborative Filtering

#make predictions on test set
pred = matrix(0,nrow = nrow(test_mat),ncol = ncol(test_mat))
num_pred = 0
cum_err = 0
num_empty = 0
for(i in 1:nrow(test_mat)){
  for(j in 1:nrow(test_mat)){
    
    #if there is a rating for this user-item pair
    if(test_mat[i,j] > 0){
      #iterate number of predictions counter
      num_pred = num_pred + 1
      
      #find the ratings prediction from similar users that have rated the movie
      prediction = (test_mat[nn[i,],j])
      prediction = prediction[which(prediction > 0)]
      if(length(prediction>0)){ #count only those users that have rated the movie
        pred[i,j] = mean(prediction)
      } else { #if no similar users rated the movie, use the average item rating
        num_empty = num_empty + 1
        pred[i,j] = mean(test_mat[which(test_mat[,j]>0),j])
      }
      
      #find the cumulative absolute error
      cum_err = cum_err + abs(test_mat[i,j] - pred[i,j])
    }
    
  }
}
#find the mean absolute error
MAE = cum_err/num_pred
percent_empty = num_empty/num_pred
