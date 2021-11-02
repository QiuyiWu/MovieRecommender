#This script has an implementation of a memory-based collaborative filtering method
#which uses isomap to find distances between users 
#instead of similarity measures

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

#find distances between users in unreduced space in a format that the isomap function accepts
space_dist = dist(train_mat)
#apply isomap, setting number of dimensions to 3
library(vegan)
L = isomap(space_dist, ndim = 50, k=25)
plot(L)
#find the coordinates of points in reduced space
coord = L$points 

#find the distance between all users in reduced space
dist_mat = matrix(0,nrow=num_user,ncol = num_user)
for(i in 1:nrow(dist_mat)){
  for(j in 1:nrow(dist_mat)){
    dist_mat[i,j] = sum((coord[i,] - coord[j,])^2)
  }
}
#set diagonal equal to the largest value for computaitonal reasons
diag(dist_mat) = max(dist_mat)


#set number of nearest neighbors to find
k = 25

#find k nearest neighbors of each user
nn = matrix(0,nrow = num_user,ncol = k)
for(i in 1:nrow(nn)){
  temp = order(dist_mat[i,],decreasing = FALSE)
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
