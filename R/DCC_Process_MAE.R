D = 56


MAE_vec = rep(0,D)
for(d in 1:D){
  if(d!=0){
    load(paste("K_MAE_D_",toString(d+1),".rda",sep=""))
    MAE_vec[d] = MAE 
  }
}

save(MAE_vec,file = "4k_it_demean_K_2_rand_BPMF_MAE_vec.rda")
# 
# rMSE_vec = rep(0,D)
# for(d in 1:D){
#   if(d!=7){
#   load(paste("rMSE_D_",toString(d+1),".rda",sep=""))
#   rMSE_vec[d] = rMSE
#   }
# }
# 
# save(rMSE_vec,file = "4k_it_demean_BBMMBPMF_rMSE_vec.rda")
