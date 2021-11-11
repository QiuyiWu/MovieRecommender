D = 41

MAE_vec = rep(0,D)
for(d in 1:D){
  load(paste("MAE_D_",toString(d+1),".rda",sep=""))
  MAE_vec[d] = MAE
}

save(MAE_vec,file = "MAE_vec.rda")