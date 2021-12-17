#plot MAE with Mean vs without Mean for 3k burnin its and 1k its
load("../Data Structures/4k_it_MAE_vec.rda")
plot(2:(length(MAE_vec)+1),MAE_vec, type= "l",col = "red",
     xlab = "D",ylab = "MAE")
load("../Data Structures/4k_it_demean_MAE_vec.rda")
lines(2:(length(MAE_vec)+1),MAE_vec, col = "blue")
load("../Data Structures/4k_it_demed_MAE_vec.rda")
lines(2:(length(MAE_vec)+1),MAE_vec, col = "green")
legend(39, 3.5, legend=c("Standard", "de-Mean","de-Med"),
       col=c("red", "blue","green"), lty=c(1,1,1), cex=0.8)
#zoom in?
load("../Data Structures/4k_it_demean_MAE_vec.rda")
plot(2:(length(MAE_vec)+1),MAE_vec, type = "l",
     col = "blue", ylim = c(0.85,1),
     xlab = "D",ylab = "MAE")
load("../Data Structures/4k_it_demed_MAE_vec.rda")
lines(2:(length(MAE_vec)+1),MAE_vec, col = "green")
legend(15, 0.95, legend=c("de-Mean","de-Med"),
       col=c("blue","green"), lty=c(1,1), cex=0.8)

#plot MAE for BPMF as D changes with 3k burnin its and 1k its
load("../Data Structures/4k_it_demean_MAE_vec.rda")
plot(2:(length(MAE_vec)+1),MAE_vec, type= "l",col = "red",
     xlab = "D",ylab = "MAE",
     xlim=c(2,46),ylim=c(0.85,0.93))
load("../Data Structures/4k_it_demean_kmean_MAE_vec.rda")
lines(2:(length(MAE_vec)+1),MAE_vec, col = "blue")
load("../Data Structures/4k_it_demean_kmean_item_const_MAE_vec.rda")
lines(2:(length(MAE_vec)+1),MAE_vec, col = "green")
legend(30, 0.92, legend=c("Standard", "Independent k-means", "Linked k-means"),
       col=c("red", "blue","green"), lty=c(1,1,1), cex=0.8)

#load same datastructures but demenaed and with k_means
load("../Data Structures/4k_it_demean_kmean_MAE_vec.rda")

load("../Data Structures/4k_it_demean_kmean_rMSE_vec.rda")

#load and plot regular vs kmean for dmeaned data
load("../Data Structures/4k_it_demean_MAE_vec.rda")
plot(2:(length(MAE_vec)+1),MAE_vec,col="red",type="l")
load("../Data Structures/4k_it_demean_kmean_MAE_vec.rda")
lines(2:(length(MAE_vec)+1),MAE_vec,col="blue",type="l")

#plot MAE for memory-based collaborative filtering
load("../Data Structures/CF_MAE.rda")
abline(h = MAE)

############
# BBMMBPMF #
############

#first plot MAE from regular demeaned BPMF
load("../Data Structures/4k_it_demean_MAE_vec.rda")
plot(2:(length(MAE_vec)+1),MAE_vec, type = "l",
     col = "red", ylim = c(0.85,1),
     xlab = "D",ylab = "MAE")
#plot MAE from demeaned BBMMBPMF
load("../Data Structures/4k_it_demean_BBMMBPMF_MAE_vec.rda")
lines(2:(length(MAE_vec)+1),MAE_vec, col = "blue")
#plot MAE from kmean BPMF
load("../Data Structures/4k_it_demean_kmean_MAE_vec.rda")
lines(2:(length(MAE_vec)+1),MAE_vec, col = "green")
#plot MAE from kmean BPMF with fake cluster
load("../Data Structures/4k_it_demean_K_2_rand_BPMF_MAE_vec.rda")
lines(2:(length(MAE_vec)+1),MAE_vec, col = "black")
legend(30, 0.92, legend=c("demeaned BPMF", "demeaned BBMM BPMF", "demeaned kBPMF","demeaned kBPMF fake clusters"),
       col=c("red", "blue","green","black"), lty=c(1,1,1,1), cex=0.8)

#examine z
load("../Data Structures/post_z.rda")
matplot(t(z_mat),type="l")


#examine BBMM split on different random seeds
load("../Data Structures/phi_mat.rda")
for(j in 1:ncol(phi_mat)){
        phi_mat[,j] = sort(phi_mat[,j])
}
matplot(t(phi_mat),type="l")
