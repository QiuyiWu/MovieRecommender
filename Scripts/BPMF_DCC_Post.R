load("../Data Structures/MAE_vec.rda")

#MAE for BPMF as D changes
plot(2:42,MAE_vec,type = "l")
#plot MAE for memory-based collaborative filtering
load("../Data Structures/CF_MAE.rda")
abline(h = MAE)