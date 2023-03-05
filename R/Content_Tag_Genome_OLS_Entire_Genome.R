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

#make dataframe for training set
response = NULL
users = NULL
x = NULL
for(i in 1:nrow(train_mat)){
  #find which items this user has rated
  train_items = which(train_mat[i,]>0)
  
  #record ratings of rated items
  response = c(response,train_mat[i,train_items])
  
  #make predictor variable for users
  users = c(users,rep(i,length(train_items)))
  
  #make predictor variables for isomap
  x = rbind(x,tag_mat[train_items,])
}
df = data.frame(response,x)

model = lm(response ~., data = df)
summary(model)

#make dataframe for test set
response = NULL
users = NULL
x = NULL
for(i in 1:nrow(test_mat)){
  #find which items this user has rated
  test_items = which(test_mat[i,]>0)
  
  #record ratings of rated items
  response = c(response,test_mat[i,test_items])
  
  #make predictor variable for users
  users = c(users,rep(i,length(test_items)))
  
  #make predictor variables for isomap
  x = rbind(x,tag_mat[test_items,])
}
df = data.frame(response,x)

test_pred = predict(model,df)

mean(abs(response-test_pred))
