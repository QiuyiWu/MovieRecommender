library(doParallel)
library(gsrs)

setwd("~/Documents/Research/Computational Advertising/ml-25m")
ratings <- read.csv("ratings.csv", header = T)
genomescores <- read.csv("genome-scores.csv", header = T)

# ratings: use subset of movies we should be limiting ourselves 
# to the set of movieIDs which appear in the genome-scores.csv file.
# ratings has 59047 unique movies; genomescores has 13816 unique movies

# subset of ratings from genomescores: 13816 unique movies, total 24.67M movies
ratings <- ratings[which(ratings$movieId %in% unique(genomescores$movieId)),]

# take a look at first 1000 obs from ratings
ratings2 <- ratings[1:1000,1:3]

registerDoParallel(cores=detectCores()-1) #7 cores
getDoParWorkers()

# Initialization Parameters
K= 60
B= 10
C= 10
lambda = 2
max_iter = 1 # usually more than 10;
tol_1=1e-1
tol_2=1e-1
# Train Test Split
N=dim(ratings3)[1]
test_rate = 0.3
train.row = sample(1:nrow(ratings3), floor((1 - test_rate) * N))
test.row = sample(1:nrow(ratings3), floor( test_rate * N))
train.data=ratings3[train.row,1:3]
test.data=ratings3[test.row,1:3]
# Call gssvd function
start_time <- Sys.time()
a = gssvd(train=train.data, test=test.data, B=B, C=C, K=K,
          lambda=lambda, max_iter=max_iter, verbose=1)
end_time <- Sys.time()
end_time - start_time

stopImplicitCluster()
# Output the result
a$RMSE_Test
head(a$P)
head(a$Q)
head(a$S)
head(a$T)

