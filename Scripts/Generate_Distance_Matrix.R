#This script takes the tag matrix and finds the Euclidean distance between all movies in 
#isomap reduced space

#load the isomap reduced space
load("../Data Structures/iso_points.rda")

#compute the distance between all movies
#Standard Euclidean distance in R^n
dist_mat = matrix(0,nrow=nrow(X),ncol = nrow(X))
for(i in 1:nrow(dist_mat)){
  for(j in 1:ncol(dist_mat)){
    dist_mat[i,j] = sqrt(sum((X[i,]-X[j,])^2))
  }
  print(i/nrow(dist_mat))
}

#save the distance matrix
save(dist_mat,file="../Data Structures/dist_mat.rda")
