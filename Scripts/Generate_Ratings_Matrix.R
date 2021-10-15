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

#generate a ratings matrix for the specified users and all the movies
unique_movies = sort(unique(sub_ratings$movieId))
rat_mat = matrix(0,nrow = length(users), ncol = length(unique_movies))
for(i in 1:nrow(rat_mat)){
  user_rows = which(sub_ratings$userId == i)
  movies = sub_ratings$movieId[user_rows]
  for(j in 1:length(movies)){
    rat_mat[i,which(unique_movies == movies[j])] = sub_ratings$rating[which(sub_ratings$movieId[user_rows] == movies[j])]
  }
}

#save the ratings matrix for future use
save(rat_mat,file="../Data Structures/rat_mat.rda")
