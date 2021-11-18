D = 46

MAE_vec = rep(0,D)
for(d in 1:D){
  load(paste("MAE_D_",toString(d+1),".rda",sep=""))
  MAE_vec[d] = MAE
}

save(MAE_vec,file = "4k_it_MAE_vec.rda")

rMSE_vec = rep(0,D)
for(d in 1:D){
  load(paste("rMSE_D_",toString(d+1),".rda",sep=""))
  rMSE_vec[d] = rMSE
}

save(rMSE_vec,file = "4k_it_rMSE_vec.rda")
