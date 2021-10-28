#This script attempts to cluster users using a k-mean algorithm
#this is equivalent to clustering users to minimize Euclidean distance to cluster centroids

library(stats)

#load the cluster matrix
load("../Data Structures/cluster_mat.rda")

#find within sum of squares for a varying number of clusters
clusters = 2:25 #the number of clusters to attempt
WSS = rep(0,length(clusters)) #The within sum of squares for each cluster
for(i in clusters){
  out = kmeans(cluster_mat,centers = i, nstart = 25)
  WSS[i-1] = sum(out$withinss)
  print(i/max(clusters))
}

#plot the WSS curve
plot(clusters,WSS,type="l")
