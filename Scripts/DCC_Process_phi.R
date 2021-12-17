num = 200

phi_mat = matrix(0,nrow = 2,ncol = num)

for(i in 1:num){
  load(paste("phi_post_",toString(i),".rda",sep=""))
  phi_mat[,i] = phi_post
  print(i)
}


save(phi_mat,file = "phi_mat.rda")