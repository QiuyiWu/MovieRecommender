#This script attempts to cluster users using a k-mean algorithm
#this is equivalent to clustering users to minimize Euclidean distance to cluster centroids

library(stats)

#load the cluster matrix
load("../data/cluster_mat.rda")

#find within sum of squares for a varying number of clusters
clusters = 2:25 #the number of clusters to attempt
WSS = rep(0,length(clusters)) #The within sum of squares for each cluster
for(i in clusters){
  out = kmeans(cluster_mat,centers = i, nstart = 25)
  WSS[i-1] = sum(out$withinss)
  print(i/max(clusters))
}

#plot the WCSS curve
plot(clusters,WSS,type="l",
     xlab = "Number of Clusters",
     ylab = "WCSS",
     main = "k-means")
#not sure if there is an elbow - there might be aroudn 7or so?

#run at elbow
out = kmeans(cluster_mat,centers = 5, nstart = 10)

#find the cluster ids
clust_id = out$cluster

#look at cluster numbers 
table(clust_id)

#Find the user vectors in each cluster
clusters = NULL
for(i in 1:max(clust_id)){
  clusters[[i]] = cluster_mat[which(clust_id == i),]
}

#recover movie names
mov = read.csv("../../Data/movies.csv")

#find the top movies in each cluster
movies = NULL
for(i in 1:length(clusters)){
  #find how many times the movies in each cluster were rated
  subset = apply(clusters[[i]],2,sum)
  #sort this in decreasing fashion
  temp = sort(subset,decreasing = TRUE)
  #recover movie id's as strings and then integers
  mov_id_str = names(temp)
  mov_id = strtoi(mov_id_str)
  #subset to top 5 movies in each cluster
  mov_id = mov_id[1:5]
  
  names = NULL
  for(j in 1:length(mov_id)){
    names = c(names,mov[which(mov$movieId == mov_id[j]),2])
  }
  
  movies[[i]] = names
}
