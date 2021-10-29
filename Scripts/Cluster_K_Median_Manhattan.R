#This script attempts to cluster users using a k-mediod algorithm
#this is equivalent to clustering users to minimize manhattan distance to a cluster centroid

library(cluster)
library(factoextra)

#load the cluster matrix
load("../Data Structures/cluster_mat.rda")

#find within sum of squares for a varying number of clusters
clusters = 2*2:25 #the number of clusters to attempt
WCMD = rep(0,length(clusters)) #The within cluster manhattan distance
for(k in clusters){
  #run the k-medians algorithm
  out = pam(cluster_mat,k,metric = "manhattan")
  
  #extract cluster medians and cluster assignemnts
  clust_mediod = out$medoids
  clust_id = out$clustering
  
  #calculate within-cluster manhattan distance
  total = 0
  for(j in 1:k){
    subset = cluster_mat[which(clust_id == j),]
    if(is.null(nrow(subset)) == TRUE){
        total = total + sum(abs(subset - clust_mediod[j,]))
    } else {
      for(n in 1:nrow(subset)){
        total = total + sum(abs(subset[n,] - clust_mediod[j,]))
      }
    }
  }
  
  #record within-cluster manhattan distance
  WCMD[k-1] = total
  print(k/max(clusters))
}


#plot the WSS curve
plot(clusters,WCMD[clusters-1],type="l")
#there might be an elbow around 10-15?

#run at elbow
out = pam(cluster_mat,15,metric = "manhattan")
#recover cluster ids
clust_id = out$clustering
table(clust_id)
#only 5 non-trivial clusters...
#might be the case that there are a lot of movies which are just off on their own?

#Find the user vectors in each cluster
clusters = NULL
for(i in 1:max(clust_id)){
  clusters[[i]] = cluster_mat[which(clust_id == i),]
}

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
  mov_id = mov_id[1:10]
  
  names = NULL
  for(j in 1:length(mov_id)){
    names = c(names,mov[which(mov$movieId == mov_id[j]),2])
  }
  
  movies[[i]] = names
}
