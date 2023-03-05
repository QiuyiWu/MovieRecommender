#This script turns a ratings matrix into a clustering matrix
#It replaces the ratings with 1 
#a user-item pair is rated if the entry is 1, 0 otherwise.

#load the ratings matrix
load("../data/rat_mat.rda")

#generate new data scture and repleace each instance of a rating with a 1
cluster_mat = rat_mat
cluster_mat[which(cluster_mat>0)] = 1

#save the clustering matrix
save(cluster_mat,file = "../data/cluster_mat.rda")