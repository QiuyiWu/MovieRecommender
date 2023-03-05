#This script has a simple implementation of a memory-based collaborative filtering method

#set random seed
set.seed(1)

#load a ratings matrix
load("../data/rat_mat.rda")
#recover the number of users and items
num_user = nrow(rat_mat)
num_item = ncol(rat_mat)

#load training and test sets
load("../data/train_R_nm.rda")
load("../data/test_R_nm.rda")
#rename so it works with old code
train_mat = train_R_nm
test_mat = test_R_nm
#remove old data structures
rm(train_R_nm)
rm(test_R_nm)

#define a similarity function - here cosine similarity
sim = function(a,b){
  return(sum(a*b)/(sqrt(sum(a^2))*sqrt(sum(b^2))))
}

#find similarities between each user
sim_mat = matrix(0,nrow=num_user,ncol = num_user)
for(i in 1:nrow(sim_mat)){
  for(j in 1:nrow(sim_mat)){
    sim_mat[i,j] = sim(train_mat[i,],train_mat[j,])
  }
}

#set the diagonals of sim_mat to 0
#set similarity between a vector and itself to 0
#makes finding maximums easier
diag(sim_mat) = 0

#set number of nearest neighbors to find
k = 25

#find k nearest neighbors of each user
nn = matrix(0,nrow = num_user,ncol = k)
for(i in 1:nrow(nn)){
  temp = order(sim_mat[i,],decreasing = TRUE)
  nn[i,] = temp[1:k]
}

#make predictions on test set
pred = matrix(0,nrow = nrow(test_mat),ncol = ncol(test_mat))
num_pred = 0
cum_err = 0
num_empty = 0
#find global average for impuration
glob_avg = mean(train_mat[train_mat>0])
#find user average for imputation
user_avg = rep(0,num_user)
for(n in 1:num_user){
  user_avg[n] = mean(train_mat[n,train_mat[n,]>0])
}
#find item average for imputaiton
item_avg = rep(0,num_item)
for(m in 1:num_item){
  item_avg[m] = mean(train_mat[train_mat[,m]>0,m])
}
#some items are NaN - not present in training set - 1959!
sum(is.nan(item_avg))

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
      } else { #if no similar users rated the movie, use the average rating from the training set
        num_empty = num_empty + 1
        # pred[i,j] = item_avg[j] #item mean
        # pred[i,j] = user_avg[i] #user mean
        pred[i,j] = glob_avg
      }
      
      #find the cumulative absolute error
      cum_err = cum_err + abs(test_mat[i,j] - pred[i,j])
    }
    
  }
}
#find the mean absolute error
MAE = cum_err/num_pred
percent_empty = num_empty/num_pred

#if we dont do imputation for empty predictions then MAE is 3.167...bad method
#impute for item mean: NaN because some items don't appear in the training set?
#if we imput for user mean then we have 0.9543072
#if we use the global mean then the MAE is 0.953949 - probably recommend this?

save(MAE,file = "../data/CF_MAE.rda")
