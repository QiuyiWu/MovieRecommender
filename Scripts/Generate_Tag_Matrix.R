# This script generates the tag matrix
#each row is a movie, each column is a tag
# the entries are the relevance scores assigned to that movie by that tag

#load genome scores
#this assumes a certain directory structure where the movielens data is stored in a Data
#directory which is in the same file as the github directory
#also this assumes that we are in the Scripts directory
genome = read.csv("../../Data/genome-scores.csv")

#load the ratings matrix to find the movies we want for this batch of users
load("../Data Structures/rat_mat.rda")
movies = strtoi(colnames(rat_mat))

#remove unnecessary data structures
rm(rat_mat)

#find number of movies and tags
num_mov = length(movies)
tags = unique(genome$tagId)
num_tag = length(unique(genome$tagId))

#generate the tag matrix
#rows are movies, columns are tags
tag_mat = matrix(0,nrow = num_mov, ncol = num_tag)
for(i in 1:num_mov){
  tag_mat[i,] = genome$relevance[which(genome$movieId == movies[i])]
  print(i/num_mov)
}

#rename rows to correspond to movies
rownames(tag_mat) = movies

#save the tag matrix for future use
save(tag_mat,file="../Data Structures/tag_mat.rda")
