library(vegan)#isomap library

#This script has an implementation of a content-based filtering method which uses
#OLS regressions, one independent for each user, on a space of isomap-reduced tag genome space
#as covariates and with ratings as predictions.  

#set random seed
set.seed(1)

#load a ratings matrix
load("../data/rat_mat.rda")
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

#load the tag matrix
load("../data/tag_mat.rda")

#load the iso object
load("../data/iso_points.rda")


#train OLS linear models on the training set with iso coords as covariates
#then evaluate on the test set and record Absolute error as well as number of points
cum_err = 0 #record cumulative absolute error
num_pred = 0 #record number of predictions made
avg_err = 0 #record the cumulative error from using user average as the guess
for(i in 1:nrow(rat_mat)){
  #find the movies in the training set
  train_items = which(train_mat[i,] > 0)
  #find the rating data
  y = train_mat[i,train_items]
  #find the covariate data
  x = data.frame(X[train_items,])
  #run a OLS linear model
  model = lm(y ~ x$Dim1 + x$Dim2 + x$Dim3)
  
  avg = mean(y)
  
  #find the test items for the user
  test_items = which(test_mat[i,]>0)
  #find the ratings
  y = test_mat[i,test_items]
  #predict on new values
  x = data.frame(X[test_items,])
  
  #update cumulative absolute error and number of predictions
  cum_err = cum_err + sum(abs(y - predict(model, x)))
  avg_err = avg_err + sum(abs(y - avg))
  num_pred = num_pred + length(y)
  
}

MAE = cum_err/num_pred
MAE
avg_err/num_pred