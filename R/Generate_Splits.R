#set my directory - only good on my computer
setwd("C:/Users/psmit/Desktop/Research/Recommeder Systems/MovieRecommender/Scripts")

#load data
ratings = read.csv("../../Data/ratings.csv", head = T)
genomescores <- read.csv("../../Data/genome-scores.csv", header = T)

# ratings: use subset of movies we should be limiting ourselves 
# to the set of movieIDs which appear in the genome-scores.csv file.
# ratings has 59047 unique movies; genomescores has 13816 unique movies

#define some data strucutres
movies = unique(genomescores$movieId)
M = length(movies)
row_ind = rep(0,nrow(ratings))
for(m in 1:nrow(ratings)){
  if(ratings[m,2] %in% movies){
    row_ind[m] = 1
  }
  if( m %%10000 == 0){
    print(m)
  }
}

# subset of ratings from genomescores: 13816 unique movies, total 24.67M movies
ratings = ratings[row_ind,]

#remove things and clear memory
rm(row_ind)
rm(genomescores)
gc()

#find unique users
users = unique(ratings$userId)

set.seed(1)
#randomize order of users
users = sample(users,length(users))

#select users until we have 500,000 ratings
num_rat = 500000
cum_rat = 0
cum_ind = NULL
for(i in 1:length(users)){
  #find the user id
  user_ind = users[i]
  
  #find which rows in teh ratings matrix belong to this user and record
  ind = which(ratings$userId == user_ind)
  cum_ind = c(cum_ind,ind)
  
  #find the cumulative number of ratings and print
  cum_rat = cum_rat + length(ind)
  print(cum_rat)
  
  #break if we have enough
  if(cum_rat>num_rat){
    break
  }
}

#record results in smaller dataset
data = ratings[cum_ind,]

#remove ratins
rm(ratings)
gc()

#find a timestap to subset the data at
quantile(data$timestamp,0.75)

#save list of users
write.csv(users[1:i],file = "../data/Users.csv",row.names = FALSE)
