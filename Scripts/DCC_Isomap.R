library(vegan)#isomap library

#This script exists to use the isomap algorithm on the tag genome matrix
#this can't be run on my laptop so I have to run it on the Duke Computing Cluster

#load the tag matrix
load("tag_mat.rda")

#prepare the tag data for isomap
tag_dist = dist(tag_mat)
save(tag_dist,file = "tag_dist.rda")
#do isomap on the tag data
L = isomap(tag_dist, ndim=3, k=5)
X = L$points
#save
save(L,file = "iso.rda")
save(X,file = "iso_points.rda")