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


# pick user at random, then got all ratings from those users and then stop till ratings number reach 10,000
ratings2 <- NULL
i=oldsampleSet=0
while (i <10000) {
  sampleSet <- sample(unique(ratings$userId), 1)
  ratings2 <- rbind(ratings2, ratings[which(ratings$userId %in% sampleSet), ])
  if (sampleSet == oldsampleSet){
    break
  }
  i = nrow(ratings2)
  oldsampleSet <- sampleSet
}
dim(ratings2)
save(ratings2, file = "sampled_ratings2.RData")

