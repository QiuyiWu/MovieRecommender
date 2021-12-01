D = 60


MAE_vec = rep(0,D)
for(d in 1:D){
  if(d!=2){
    load(paste("MAE_D_",toString(d+1),".rda",sep=""))
    MAE_vec[d] = MAE 
  }
}

save(MAE_vec,file = "4k_it_demean_MAE_vec.rda")

rMSE_vec = rep(0,D)
for(d in 1:D){
  if(d!=2){
  load(paste("rMSE_D_",toString(d+1),".rda",sep=""))
  rMSE_vec[d] = rMSE
  }
}

save(rMSE_vec,file = "4k_it_demean_rMSE_vec.rda")
