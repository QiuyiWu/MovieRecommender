#load ratings
#this assumes a certain directory structure where the movielens data is stored in a Data
#directory which is in the same file as the github directory
#also this assumes that we are in the Scripts directory
ratings = read.csv("../../Data/ratings.csv")

#subset to user IDS
userIDs = ratings$userId

#find rows corresponding to a subset of users
users = 1:1000
user_rows = NULL
for(i in 1:length(users)){
  user_rows = c(user_rows,which(userIDs == users[i]))
}
sub_ratings = ratings[user_rows,]

#remove unnecessary data structures 
rm(ratings)
rm(userIDs)
rm(user_rows)

#For ease of comparison, we want to find movies which have tag genome data
#load genome scores
#this assumes a certain directory structure where the movielens data is stored in a Data
#directory which is in the same file as the github directory
#also this assumes that we are in the Scripts directory
genome = read.csv("../../Data/genome-scores.csv")
#find all unique movies with tag genome data
tag_movies = unique(genome$movieId)

#manage memory
rm(genome)

#find all unique movies in the ratings data
unique_movies = sort(unique(sub_ratings$movieId))
#find the movies both in the ratings data and the tag data
desired_movies = NULL
for(i in 1:length(unique_movies)){
  if(unique_movies[i] %in% tag_movies){
    desired_movies = c(desired_movies,unique_movies[i])
  }
}

#find the ratings matrix for the movies in the tag genome dataset
rat_mat = matrix(0,nrow = length(users), ncol = length(desired_movies))
for(i in 1:nrow(rat_mat)){
  user_rows = which(sub_ratings$userId == i)
  movies = sub_ratings$movieId[user_rows]
  movies = movies[movies %in% desired_movies]
  for(j in 1:length(movies)){
    rat_mat[i,which(desired_movies == movies[j])] = sub_ratings$rating[which(sub_ratings$movieId[user_rows] == movies[j])]
  }
}

rownames(rat_mat) = users
colnames(rat_mat) = desired_movies

#save the ratings matrix for future use
save(rat_mat,file="../Data Structures/rat_mat.rda")
