#This script generates a training set and a test set based off of
# a rating matrix and a clustering matrix
#This is more effiicent than running this code in every script

set.seed(1)

#load a ratings matrix
load("../Data Structures/rat_mat.rda")
#recover the number of users and items
N = nrow(rat_mat) #number of users
M = ncol(rat_mat) #number of items

#divide data into training and test set
#for each user, put half of that user's item in the test set?
train_R_nm = matrix(0,nrow = N, ncol =M)
test_R_nm = matrix(0,nrow = N, ncol =M)
items = 1:M
for(n in 1:N){
  train_items = sample(x = items,size = round(length(items)/2), replace = FALSE)
  test_items = 1:length(items)
  test_items = test_items[-train_items]
  train_R_nm[n,train_items] = rat_mat[n,train_items]
  test_R_nm[n,test_items] = rat_mat[n,test_items]
}


#make cluster matrices for the training set
train_I_nm = train_R_nm
train_I_nm[which(train_I_nm>0)] = 1

#make cluster matrix for the test set
test_I_nm = test_R_nm
test_I_nm[which(test_I_nm>0)] = 1

#save files
save(train_R_nm,file = "../Data Structures/train_mat.rda")
save(test_R_nm,file = "../Data Structures/test_mat.rda")
save(train_I_nm, file = "../Data Structures/train_I_nm.rda")
save(test_I_nm, file = "../Data Structures/test_I_nm.rda")