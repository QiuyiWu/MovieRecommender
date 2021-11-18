load("../Data Structures/MAE_vec.rda")

#MAE for BPMF as D changes for 2k burnin its and 1k its
plot(2:42,MAE_vec,type = "l", col = "blue")

#plot MAE for BPMF as D changes with 3k burnin its and 1k its
load("../Data Structures/4k_it_MAE_vec.rda")
lines(2:(length(MAE_vec)+1),MAE_vec, col = "red")
#load and plot rMSE
load("../Data Structures/4k_it_rMSE_vec.rda")

#plot MAE for memory-based collaborative filtering
load("../Data Structures/CF_MAE.rda")
abline(h = MAE)

load("../Data Structures/kMAE_vec.rda")

