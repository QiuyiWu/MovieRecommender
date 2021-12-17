D = 35
K = 2

z_mat = matrix(0,nrow = K,ncol = D)
for(d in 1:D){
  if(d!=7){
    load(paste("post_z_",toString(d+1),".rda",sep=""))
    z_mat[,d] = apply(post_z_nk,2,mean)
  }
}

save(z_mat,file = "post_z.rda")