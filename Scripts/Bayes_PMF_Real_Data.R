#This script generates data from a Bayseian PMF distribution and then implements a Gibbs sampler
library(stats)

set.seed(1)

#load a ratings matrix
load("../Data Structures/rat_mat.rda")
#recover the number of users and items
N = nrow(rat_mat) #number of users
M = ncol(rat_mat) #number of items

#divide data into training and test set
#for each user, put half of that user's item in the test set?

#load cluster_mat
load("../Data Structures/cluster_mat.rda")

#rename to indicator functions
I_nm = cluster_mat[,] #subset for time...

#remove cluster mat
rm(cluster_mat)

#generate a list of which items each user has rated
user_rated_n = NULL
for(n in 1:N){
  user_rated_n[[n]] = which(I_nm[n,]>0)
}

#generate a list of which users have rated each item
item_rated_m = NULL
for(m in 1:M){
  item_rated_m[[m]] = which(I_nm[,m]>0)
}


train_items = sample(1:num_item,round(num_item/2))
test_items = 1:num_item
test_items = test_items[-train_items]
train_mat = rat_mat[,train_items]
test_mat = rat_mat[,test_items]

rm(rat_mat)