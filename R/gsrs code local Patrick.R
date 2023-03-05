library(gsrs)

setwd("C:/Users/psmit/Desktop/Research/Recommeder Systems/MovieRecommender/Scripts")
ratings = read.csv("../../Data/ratings.csv", head = T)
genomescores <- read.csv("../../Data/genome-scores.csv", header = T)

# ratings: use subset of movies we should be limiting ourselves 
# to the set of movieIDs which appear in the genome-scores.csv file.
# ratings has 59047 unique movies; genomescores has 13816 unique movies

# subset of ratings from genomescores: 13816 unique movies, total 24.67M movies
movies =  unique(genomescores$movieId)
ratings <- ratings[which(ratings$movieId %in% movies),]

# take a look at first 1000 obs from ratings
ratings2 <- ratings[1:10000,1:3]

# Initialization Parameters
K= 60
B= 10
C= 10
lambda = 2
max_iter = 1 # usually more than 10;
tol_1=1e-1
tol_2=1e-1
# Train Test Split
N=dim(ratings2)[1]
test_rate = 0.3
train.row = sample(1:N, floor((1 - test_rate) * N))
test.row = sample(1:N, floor( test_rate * N))
train.data=ratings2[train.row,1:3]
test.data=ratings2[test.row,1:3]
# Call gssvd function
start_time <- Sys.time()
a = gssvd(train=train.data, test=test.data, B=B, C=C, K=K,
          tol_1 = tol_1, tol_2 = tol_2,
          lambda=lambda, max_iter=max_iter, verbose=1)
end_time <- Sys.time()
end_time - start_time

# Output the result
a$RMSE_Test
head(a$P)
head(a$Q)
head(a$S)
head(a$T)

