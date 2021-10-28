#This script attempts to cluster users using a k-mediod algorithm
#this is equivalent to clustering users to minimize manhattan distance to a cluster centroid

library(cluster)
library(factoextra)

#load the cluster matrix
load("../Data Structures/cluster_mat.rda")

#find within sum of squares for a varying number of clusters
clusters = 2:25 #the number of clusters to attempt
WCMD = rep(0,length(clusters)) #The within cluster manhattan distance
for(i in clusters){
  #set the number of clusters
  k = clusters[i]
  
  #run the k-medians algorithm
  out = pam(cluster_mat,k,metric = "manhattan")
  
  #extract cluster medians and cluster assignemnts
  clust_mediod = out$medoids
  clust_id = out$clustering
  
  #calculate within-cluster manhattan distance
  total = 0
  for(j in 1:k){
    subset = cluster_mat[which(clust_id == j),]
    for(n in 1:nrow(subset)){
      total = total + sum(abs(subset[n,] - clust_mediod[j,]))
    }
  }
  
  #record within-cluster manhattan distance
  WCMD[i-1] = total
  print(i/max(clusters))
}


#plot the WSS curve
plot(clusters,WCMD,type="l")
