#This script has a simple implementation of a memory-based collaborative filtering method

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
