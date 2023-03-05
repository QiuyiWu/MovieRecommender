load("../Data Structures/Occ_out.rda")

post_phi_k = Occ_out$post_phi
post_U_n = Occ_out$post_U
post_theta_l = Occ_out$post_theta
post_M_m = Occ_out$post_M
post_psi_kl = Occ_out$post_psi
post_W_nm = Occ_out$post_W
post_p_kl = Occ_out$post_p


plot(true_phi_k,post_phi_k[1:K],
     xlim= c(0,1),ylim = c(0,1))
lines(0:1,0:1)

plot(true_U_n,post_U_n)

plot(true_theta_l,post_theta_l)

plot(true_M_m,post_M_m)

plot(true_psi_kl[1:K,],post_psi_kl[1:K,],
     xlim= c(0,1),ylim = c(0,1))

plot(true_W_nm,post_W_nm)

plot(true_p_kl[,],post_p_kl[1:K,],
     xlim = c(0,1),ylim = c(0,1))
