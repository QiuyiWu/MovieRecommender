D = 35

kMAE_vec = rep(0,D)
for(d in 1:D){
  if(d != 4){
    load(paste("K_MAE_D_",toString(d+1),".rda",sep=""))
    kMAE_vec[d] = MAE
  }
}

save(kMAE_vec,file = "kMAE_vec.rda")