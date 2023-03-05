#this script includes a baseline measure for recommender system accuracy
#the predicted rating is the mean of all ratings,
#plus the user difference from teh mean
#plus the item difference

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

#find the overall mean
avg = mean(train_mat[which(train_mat>0)])

#find the user mean
user_mean = NULL
for(i in 1:nrow(train_mat)){
  user_mean[i] = mean(train_mat[i,which(train_mat[i,]>0)]) - avg
}

#find the item mean
item_mean = NULL
for(i in 1:ncol(train_mat)){
  item_mean[i] = mean(train_mat[which(train_mat[,i]>0),i]) - avg
}


#make predictions on test set
num_pred = 0
cum_err = 0
avg_err = 0
user_err = 0 
item_err = 0
for(i in 1:nrow(test_mat)){
  for(j in 1:ncol(test_mat)){
    
    #if there is a rating for this user-item pair
    if(test_mat[i,j] > 0){
      #iterate number of predictions counter
      num_pred = num_pred + 1
      
      #find the ratings prediction from average, user, item effects
      prediction = avg + user_mean[i] + item_mean[j]
      
      #find the cumulative absolute error
      cum_err = cum_err + abs(test_mat[i,j] - prediction)
      
      #find the ratings prediction from average
      prediction = avg 
      
      #find the cumulative absolute error
      avg_err = avg_err + abs(test_mat[i,j] - prediction)
      
      #find the ratings prediction from average plus user effects
      prediction = avg + user_mean[i]
      
      #find the cumulative absolute error
      user_err = user_err + abs(test_mat[i,j] - prediction)
      
      #find the ratings prediction from average plus item effects
      prediction = avg + item_mean[j]
      
      #find the cumulative absolute error
      item_err = item_err + abs(test_mat[i,j] - prediction)
    }
    
  }
}
#find the mean absolute error
MAE = cum_err/num_pred
MAE
avg_err/num_pred
user_err/num_pred
item_err/num_pred