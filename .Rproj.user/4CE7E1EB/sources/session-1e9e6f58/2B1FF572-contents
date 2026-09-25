beta_increment_ord <- function(pi_T_hat, pi_C_hat,
                               rho_matrix,
                               cluster_size,
                               number_cluster,
                               y_by_cluster,
                               treat_cluster) {
  
  # K = number of indicator components used (i.e., K = (#categories - 1))
  K <- length(pi_T_hat)
  if (length(pi_C_hat) != K) stop("pi_T_hat and pi_C_hat must have same length K.")
  if (!all(dim(rho_matrix) == c(K, K))) stop("rho_matrix must be K x K.")
  if (length(treat_cluster) != number_cluster) stop("treat_cluster length must equal number_cluster.")
  
  # ---- build arm-specific pieces ----
  a_mat_T <- make_a_mat(pi_T_hat)
  a_mat_C <- make_a_mat(pi_C_hat)
  
  Q_mat_T <- make_Q_mat_fast(pi_T_hat)
  Q_mat_C <- make_Q_mat_fast(pi_C_hat)
  
  d_mat_T <- make_d_mat(pi_T_hat, treat_last_col = TRUE)    # K x (K+1)
  d_mat_C <- make_d_mat(pi_C_hat, treat_last_col = FALSE)   # last col 0
  
  # ---- pieces used in dv_* (cluster-level) ----
  inv_part_T <- solve(Q_mat_T + (cluster_size - 1) * rho_matrix)
  inv_part_C <- solve(Q_mat_C + (cluster_size - 1) * rho_matrix)
  
  dv_T <- t(d_mat_T) %*% a_mat_T %*% inv_part_T %*% a_mat_T   # (K+1) x K
  dv_C <- t(d_mat_C) %*% a_mat_C %*% inv_part_C %*% a_mat_C   # (K+1) x K
  
  dvd_T <- dv_T %*% d_mat_T   # (K+1) x (K+1)
  dvd_C <- dv_C %*% d_mat_C   # (K+1) x (K+1)
  
  inv_sum_DVD <- 2 / cluster_size / number_cluster * solve(dvd_T + dvd_C)
  
  # ---- big-block pieces (within cluster) ----
  oneJ <- matrix(rep(1, cluster_size), ncol = 1)
  
  D_T <- kronecker(oneJ, d_mat_T)          # (J*K) x (K+1)
  D_C <- kronecker(oneJ, d_mat_C)          # (J*K) x (K+1)
  
  inv_A_T <- kronecker(diag(cluster_size), a_mat_T)  # (J*K) x (J*K)
  inv_A_C <- kronecker(diag(cluster_size), a_mat_C)
  
  U <- (1 / cluster_size) * matrix(1, nrow = cluster_size, ncol = cluster_size)
  U_vert <- diag(cluster_size) - U
  
  inv_R_T <- kronecker(U_vert, solve(Q_mat_T - rho_matrix)) +
    kronecker(U,      solve(Q_mat_T + (cluster_size - 1) * rho_matrix))
  
  inv_R_C <- kronecker(U_vert, solve(Q_mat_C - rho_matrix)) +
    kronecker(U,      solve(Q_mat_C + (cluster_size - 1) * rho_matrix))
  
  inv_V_T <- inv_A_T %*% inv_R_T %*% inv_A_T
  inv_V_C <- inv_A_C %*% inv_R_C %*% inv_A_C
  
  # ---- data aggregation across clusters ----
  Y_mat <- do.call(cbind, y_by_cluster)  # (J*K) x (#clusters)
  if (!all(dim(Y_mat) == c(cluster_size * K, number_cluster))) {
    stop("Each element of y_by_cluster must be length cluster_size*K, and list length = number_cluster.")
  }
  
  Y_C <- Y_mat[, treat_cluster == 0, drop = FALSE] # Each column is one cluster’s stacked vector of binary indicators.
  Y_T <- Y_mat[, treat_cluster == 1, drop = FALSE]
  
  sum_Y_T <- if (ncol(Y_T) > 0) rowSums(Y_T) else rep(0, nrow(Y_mat)) # calculate the sum of stacked indicator vectors across all treatment clusters.
  sum_Y_C <- if (ncol(Y_C) > 0) rowSums(Y_C) else rep(0, nrow(Y_mat))
  
  DVY_T <- t(D_T) %*% inv_V_T %*% sum_Y_T
  DVY_C <- t(D_C) %*% inv_V_C %*% sum_Y_C
  sum_DVY <- DVY_T + DVY_C
  
  # ---- mean part ----
  dvpi_T <- dv_T %*% pi_T_hat
  dvpi_C <- dv_C %*% pi_C_hat
  sum_DVpi_hat <- cluster_size * number_cluster / 2 * (dvpi_T + dvpi_C)
  
  part2 <- sum_DVY - sum_DVpi_hat
  increment <- inv_sum_DVD %*% part2
  
  list(
    increment = increment,
    part2 = part2,
    d_mat_T = d_mat_T, d_mat_C = d_mat_C,
    Q_mat_T = Q_mat_T, Q_mat_C = Q_mat_C
  )
}

#### update beta by ISLR
beta_increment = function(pi_T_hat,pi_C_hat,
                          rho_matrix,
                          cluster_size ,
                          number_cluster,
                          y01_by_cluster){
  pi_T1 = pi_T_hat[1]
  pi_T2 = pi_T_hat[2]
  pi_C1 = pi_C_hat[1]
  pi_C2 = pi_C_hat[2]
  
  cdf_T2 = pi_T1 + pi_T2
  cdf_C2 = pi_C1 + pi_C2
  
  a_T1 = (pi_T1*(1-pi_T1))^(-1/2); a_T2 = (pi_T2*(1-pi_T2))^(-1/2)
  a_C1 = (pi_C1*(1-pi_C1))^(-1/2); a_C2 = (pi_C2*(1-pi_C2))^(-1/2)
  a_mat_T = matrix(c( a_T1,0,0, a_T2),nrow=2)
  a_mat_C = matrix(c( a_C1,0,0, a_C2),nrow=2)
  
  r_T = -pi_T1* pi_T2 *a_T1*a_T2 #-pi_T1* pi_T2 / sqrt(pi_T1*(1-pi_T1)*pi_T2*(1-pi_T2))
  r_C = -pi_C1* pi_C2 *a_C1*a_C2
  
  Q_mat_T = matrix(c( 1,r_T ,r_T , 1),nrow=2)
  Q_mat_C = matrix(c( 1,r_C ,r_C , 1),nrow=2)  
  
  d_mat_T = matrix(c(pi_T1*(1-pi_T1),-pi_T1*(1-pi_T1),
                     0, cdf_T2*(1-cdf_T2),
                     pi_T1*(1-pi_T1),  cdf_T2*(1-cdf_T2)-pi_T1*(1-pi_T1)),
                   ncol=3) # 2*3
  
  d_mat_C = matrix(c(pi_C1*(1-pi_C1),-pi_C1*(1-pi_C1),
                     0, cdf_C2*(1-cdf_C2),
                     0,  0) ,
                   ncol=3)
  inv_part_T = solve(Q_mat_T+(cluster_size-1)*rho_matrix)
  inv_part_C = solve(Q_mat_C+(cluster_size-1)*rho_matrix)
  
  dv_T = t(d_mat_T) %*% a_mat_T %*% inv_part_T %*%a_mat_T
  dv_C = t(d_mat_C) %*% a_mat_C %*% inv_part_C %*%a_mat_C
  
  dvd_T = dv_T %*% d_mat_T
  dvd_C = dv_C %*% d_mat_C
  
  inv_sum_DVD = 2/ cluster_size/number_cluster * solve(dvd_T  + dvd_C)
  
  D_T = kronecker(matrix(rep(1,cluster_size ),ncol=1), d_mat_T )
  inv_A_T = kronecker( diag(cluster_size), a_mat_T)
  
  D_C = kronecker(matrix(rep(1,cluster_size ),ncol=1), d_mat_C )
  inv_A_C = kronecker( diag(cluster_size), a_mat_C)
  
  U = 1/cluster_size * matrix(1,ncol=cluster_size,nrow=cluster_size)
  U_vert = diag(cluster_size) - U
  # ## check
  # all(U %*% U == U)
  # all(U %*% U_vert < 1e-9)
  # all(U_vert %*% U_vert - U_vert< 1e-9)
  # all(U+ U_vert-diag(cluster_size) <1e-9)
  
  inv_R_T = kronecker( U_vert, solve(Q_mat_T-rho_matrix))  +
    kronecker( U, solve(Q_mat_T + (cluster_size-1)*rho_matrix))
  
  inv_R_C = kronecker( U_vert, solve(Q_mat_C-rho_matrix))  +
    kronecker( U, solve(Q_mat_C + (cluster_size-1)*rho_matrix))
  
  inv_V_T = inv_A_T %*% inv_R_T %*% inv_A_T
  inv_V_C = inv_A_C %*% inv_R_C %*% inv_A_C
  
  Y_mat = do.call(cbind, y01_by_cluster) 
  
  Y_C = Y_mat[,treat_cluster ==0]
  
  Y_T = Y_mat[,treat_cluster ==1]
  
  sum_Y_T = rowSums(Y_T)           
  sum_Y_C = rowSums(Y_C)  
  
  ## since D_T and inv_V_T are the same for all the treatment clusters.
  DVY_T = t(D_T) %*% inv_V_T %*% sum_Y_T 
  DVY_C = t(D_C) %*% inv_V_C %*% sum_Y_C
  
  sum_DVY = DVY_T + DVY_C
  
  dvpi_T = dv_T %*% pi_T_hat 
  dvpi_C = dv_C %*% pi_C_hat 
  sum_DVpi_hat = cluster_size * number_cluster/2 * (dvpi_T+dvpi_C)
  
  part2 = sum_DVY-sum_DVpi_hat
  
  increment = inv_sum_DVD %*% part2 
  
  return(list(increment=increment,
              part2 = part2))
}

# rho_matrix
# cluster_size 
# number_cluster
# 
# y01_by_cluster <- lapply(split(sim_data$Y, sim_data$cluster), function(yf) {
#   y <- as.character(yf)
#   c(rbind(as.integer(y == "0"), as.integer(y == "1")))
# })
# 
# delta_beta = beta_new
# 
# 
# 
# beta_init = - gee_fit0$beta
# 
# beta_trace = beta_init
# beta_old = beta_init
# beta_new = beta_init
# ### update 
# 
# 
# 
# 
# #### get pi_hat
# alpha <- unname(beta_new[grepl("^Inter:", names(beta_new))])  # length K-1
# beta_treat <- unname(beta_new["treat1"])
# treat_vec <- model.matrix(~ treat, data = sim_data)[, "treat1"]
# # cumulative probs
# eta <- outer(rep(1, nrow(sim_data)), alpha) + treat_vec %o% rep(beta_treat, length(alpha))
# cum <- plogis(eta)
# # category probs (K=3 here), pi_hat
# pi_hat <- cbind(cum[,1], cum[,2]-cum[,1], 1 - cum[,2])
# distinct(as.data.frame(treat_vec))
# pi_hat_CT = distinct(as.data.frame(pi_hat))[,-3]
# pi_hat_CT
# pi_C_hat = matrix(as.numeric(pi_hat_CT[1,]),ncol=1)
# pi_T_hat = matrix(as.numeric(pi_hat_CT[2,]),ncol=1)
# 
# #### update beta by ISLR
# beta_increment = function(pi_T_hat,pi_C_hat,
#                           rho_matrix,
#                           cluster_size ,
#                           number_cluster,
#                           y01_by_cluster){
#   pi_T1 = pi_T_hat[1]
#   pi_T2 = pi_T_hat[2]
#   pi_C1 = pi_C_hat[1]
#   pi_C2 = pi_C_hat[2]
# 
#   cdf_T2 = pi_T1 + pi_T2
#   cdf_C2 = pi_C1 + pi_C2
#   
#   a_T1 = (pi_T1*(1-pi_T1))^(-1/2); a_T2 = (pi_T2*(1-pi_T2))^(-1/2)
#   a_C1 = (pi_C1*(1-pi_C1))^(-1/2); a_C2 = (pi_C2*(1-pi_C2))^(-1/2)
#   a_mat_T = matrix(c( a_T1,0,0, a_T2),nrow=2)
#   a_mat_C = matrix(c( a_C1,0,0, a_C2),nrow=2)
#   
#   r_T = -pi_T1* pi_T2 *a_T1*a_T2 #-pi_T1* pi_T2 / sqrt(pi_T1*(1-pi_T1)*pi_T2*(1-pi_T2))
#   r_C = -pi_C1* pi_C2 *a_C1*a_C2
#   
#   Q_mat_T = matrix(c( 1,r_T ,r_T , 1),nrow=2)
#   Q_mat_C = matrix(c( 1,r_C ,r_C , 1),nrow=2)  
#   
#   d_mat_T = matrix(c(pi_T1*(1-pi_T1),-pi_T1*(1-pi_T1),
#                      0, cdf_T2*(1-cdf_T2),
#                      pi_T1*(1-pi_T1),  cdf_T2*(1-cdf_T2)-pi_T1*(1-pi_T1)),
#                    ncol=3) # 2*3
#   
#   d_mat_C = matrix(c(pi_C1*(1-pi_C1),-pi_C1*(1-pi_C1),
#                      0, cdf_C2*(1-cdf_C2),
#                      0,  0) ,
#                    ncol=3)
#   inv_part_T = solve(Q_mat_T+(cluster_size-1)*rho_matrix)
#   inv_part_C = solve(Q_mat_C+(cluster_size-1)*rho_matrix)
#   
#   dv_T = t(d_mat_T) %*% a_mat_T %*% inv_part_T %*%a_mat_T
#   dv_C = t(d_mat_C) %*% a_mat_C %*% inv_part_C %*%a_mat_C
#   
#   dvd_T = dv_T %*% d_mat_T
#   dvd_C = dv_C %*% d_mat_C
#   
#   inv_sum_DVD = 2/ cluster_size/number_cluster * solve(dvd_T  + dvd_C)
#   
#   D_T = kronecker(matrix(rep(1,cluster_size ),ncol=1), d_mat_T )
#   inv_A_T = kronecker( diag(cluster_size), a_mat_T)
#   
#   D_C = kronecker(matrix(rep(1,cluster_size ),ncol=1), d_mat_C )
#   inv_A_C = kronecker( diag(cluster_size), a_mat_C)
#   
#   U = 1/cluster_size * matrix(1,ncol=cluster_size,nrow=cluster_size)
#   U_vert = diag(cluster_size) - U
#   # ## check
#   # all(U %*% U == U)
#   # all(U %*% U_vert < 1e-9)
#   # all(U_vert %*% U_vert - U_vert< 1e-9)
#   # all(U+ U_vert-diag(cluster_size) <1e-9)
#   
#   inv_R_T = kronecker( U_vert, solve(Q_mat_T-rho_matrix))  +
#     kronecker( U, solve(Q_mat_T + (cluster_size-1)*rho_matrix))
#   
#   inv_R_C = kronecker( U_vert, solve(Q_mat_C-rho_matrix))  +
#     kronecker( U, solve(Q_mat_C + (cluster_size-1)*rho_matrix))
#   
#   inv_V_T = inv_A_T %*% inv_R_T %*% inv_A_T
#   inv_V_C = inv_A_C %*% inv_R_C %*% inv_A_C
#   
#   Y_mat = do.call(cbind, y01_by_cluster) 
#   
#   Y_C = Y_mat[,treat_cluster ==0]
#   
#   Y_T = Y_mat[,treat_cluster ==1]
#   
#   sum_Y_T = rowSums(Y_T)           
#   sum_Y_C = rowSums(Y_C)  
#   
#   ## since D_T and inv_V_T are the same for all the treatment clusters.
#   DVY_T = t(D_T) %*% inv_V_T %*% sum_Y_T 
#   DVY_C = t(D_C) %*% inv_V_C %*% sum_Y_C
#   
#   sum_DVY = DVY_T + DVY_C
#   
#   
#   dvpi_T = dv_T %*% pi_T_hat 
#   dvpi_C = dv_C %*% pi_C_hat 
#   sum_DVpi_hat = cluster_size * number_cluster/2 * (dvpi_T+dvpi_C)
#   
#   part2 = sum_DVY-sum_DVpi_hat
#   
#   increment = inv_sum_DVD %*% part2 
#   
#   return(list(increment=increment,
#               part2 = part2))
# }
# 
# beta_update = beta_increment(pi_T_hat=pi_T_hat,
#                            pi_C_hat=pi_C_hat,
#                            rho_matrix=rho_matrix,
#                            cluster_size =cluster_size ,
#                            number_cluster=number_cluster)
# score = as.numeric(beta_update$part2) 
# increment = as.numeric(beta_update$increment) 
# beta_old = beta_new
# beta_new = beta_old +as.numeric(increment) 
# 
# beta_trace = rbind(beta_trace,beta_new)
# beta_new 
# 
# 
# 
# beta_init = - gee_fit0$beta
# # beta <- beta_init
# beta_trace <- rbind(beta)   # store the initial point
# score_all <- NULL
# beta_old = beta_init
# beta_new = beta_init
# 
# for (i in 1:10){
#   
#   #### get pi_hat
#   alpha <- unname(beta_new[grepl("^Inter:", names(beta_new))])  # length K-1
#   beta_treat <- unname(beta_new["treat1"])
#   treat_vec <- model.matrix(~ treat, data = sim_data)[, "treat1"]
#   # cumulative probs
#   eta <- outer(rep(1, nrow(sim_data)), alpha) + treat_vec %o% rep(beta_treat, length(alpha))
#   print(distinct(as.data.frame(eta)))
#   cum <- plogis(eta)
#   print(distinct(as.data.frame(cum)))
#   # category probs (K=3 here), pi_hat
#   pi_hat <- cbind(cum[,1], cum[,2]-cum[,1], 1 - cum[,2])
#   print(distinct(as.data.frame(pi_hat)))
#   distinct(as.data.frame(treat_vec))
#   pi_hat_CT = distinct(as.data.frame(pi_hat))[,-3]
#   pi_hat_CT
#   pi_C_hat = matrix(as.numeric(pi_hat_CT[1,]),ncol=1)
#   pi_T_hat = matrix(as.numeric(pi_hat_CT[2,]),ncol=1)
#   
#   beta_update = beta_increment(pi_T_hat=pi_T_hat,
#                                pi_C_hat=pi_C_hat,
#                                rho_matrix=rho_matrix,
#                                cluster_size =cluster_size ,
#                                number_cluster=number_cluster)
#   score_now = as.numeric(beta_update$part2) 
#   increment = as.numeric(beta_update$increment) 
#   
#   # # Backtracking to ensure valid probs & (optionally) decreasing score
#   # upd <- update_with_backtracking(beta_old = beta_new,
#   #                                 increment = increment,
#   #                                 sim_data = sim_data,
#   #                                 score_old = score_now,
#   #                                 score_fn = NULL)  # or = score_fn
#   # beta_new <- upd$beta_new
#   beta_new <- beta_old + increment
#   beta_trace = rbind(beta_trace,beta_new)
#   print(i)
#   print(beta_old)
#   print(beta_new)
#   #print(beta_trace)
#   print(increment)
#   print(score_now)
#   score_all = rbind(score_all,score_now)
#   
#   rel_change_beta = max(abs((beta_new - beta_old)/(beta_old + .Machine$double.eps)))
#   print(rel_change_beta)
#   print(rel_change_beta<1e-5)
#   beta_old <- beta_new
# 
#   if (rel_change_beta<1e-5) break
#   # if (all(abs(score_now)<1e-4)){
#   #   break
#   # }
# }
# 
# beta_new
# beta_trace
# score_all
# 
# beta_new
# #### get pi_hat
# alpha <- unname(beta_new[grepl("^Inter:", names(beta_new))])  # length K-1
# beta_treat <- unname(beta_new["treat1"])
# treat_vec <- model.matrix(~ treat, data = sim_data)[, "treat1"]
# # cumulative probs
# eta <- outer(rep(1, nrow(sim_data)), alpha) + treat_vec %o% rep(beta_treat, length(alpha))
# print(distinct(as.data.frame(eta)))
# cum <- plogis(eta)
# print(distinct(as.data.frame(cum)))
# # category probs (K=3 here), pi_hat
# pi_hat <- cbind(cum[,1], cum[,2]-cum[,1], 1 - cum[,2])
# print(distinct(as.data.frame(pi_hat)))
# distinct(as.data.frame(treat_vec))
# pi_hat_CT = distinct(as.data.frame(pi_hat))[,-3]
# pi_hat_CT
# pi_C_hat = matrix(as.numeric(pi_hat_CT[1,]),ncol=1)
# pi_T_hat = matrix(as.numeric(pi_hat_CT[2,]),ncol=1)
# 
# sandwich_var = function(pi_T_hat,
#                         pi_C_hat,
#                         rho_matrix,
#                         cluster_size ,
#                         number_cluster){
#   pi_T1 = pi_T_hat[1]
#   pi_T2 = pi_T_hat[2]
#   pi_C1 = pi_C_hat[1]
#   pi_C2 = pi_C_hat[2]
#   
#   cdf_T2 = pi_T1 + pi_T2
#   cdf_C2 = pi_C1 + pi_C2
#   
#   a_T1 = (pi_T1*(1-pi_T1))^(-1/2); a_T2 = (pi_T2*(1-pi_T2))^(-1/2)
#   a_C1 = (pi_C1*(1-pi_C1))^(-1/2); a_C2 = (pi_C2*(1-pi_C2))^(-1/2)
#   a_mat_T = matrix(c( a_T1,0,0, a_T2),nrow=2)
#   a_mat_C = matrix(c( a_C1,0,0, a_C2),nrow=2)
#   
#   r_T = -pi_T1* pi_T2 *a_T1*a_T2 #-pi_T1* pi_T2 / sqrt(pi_T1*(1-pi_T1)*pi_T2*(1-pi_T2))
#   r_C = -pi_C1* pi_C2 *a_C1*a_C2
#   
#   Q_mat_T = matrix(c( 1,r_T ,r_T , 1),nrow=2)
#   Q_mat_C = matrix(c( 1,r_C ,r_C , 1),nrow=2)  
#   
#   d_mat_T = matrix(c(pi_T1*(1-pi_T1),-pi_T1*(1-pi_T1),
#                      0, cdf_T2*(1-cdf_T2),
#                      pi_T1*(1-pi_T1),  cdf_T2*(1-cdf_T2)-pi_T1*(1-pi_T1)),
#                    ncol=3) # 2*3
#   
#   d_mat_C = matrix(c(pi_C1*(1-pi_C1),-pi_C1*(1-pi_C1),
#                      0, cdf_C2*(1-cdf_C2),
#                      0,  0) ,
#                    ncol=3)
#   inv_part_T = solve(Q_mat_T+(cluster_size-1)*rho_matrix)
#   inv_part_C = solve(Q_mat_C+(cluster_size-1)*rho_matrix)
#   
#   dv_T = t(d_mat_T) %*% a_mat_T %*% inv_part_T %*%a_mat_T
#   dv_C = t(d_mat_C) %*% a_mat_C %*% inv_part_C %*%a_mat_C
#   
#   dvd_T = dv_T %*% d_mat_T
#   dvd_C = dv_C %*% d_mat_C
#   
#   inv_sum_DVD = 2/ cluster_size/number_cluster * solve(dvd_T  + dvd_C)
#   
#   D_T = kronecker(matrix(rep(1,cluster_size ),ncol=1), d_mat_T )
#   inv_A_T = kronecker( diag(cluster_size), a_mat_T)
#   
#   D_C = kronecker(matrix(rep(1,cluster_size ),ncol=1), d_mat_C )
#   inv_A_C = kronecker( diag(cluster_size), a_mat_C)
#   
#   U = 1/cluster_size * matrix(1,ncol=cluster_size,nrow=cluster_size)
#   U_vert = diag(cluster_size) - U
#   # ## check
#   # all(U %*% U == U)
#   # all(U %*% U_vert < 1e-9)
#   # all(U_vert %*% U_vert - U_vert< 1e-9)
#   # all(U+ U_vert-diag(cluster_size) <1e-9)
#   
#   inv_R_T = kronecker( U_vert, solve(Q_mat_T-rho_matrix))  +
#     kronecker( U, solve(Q_mat_T + (cluster_size-1)*rho_matrix))
#   
#   inv_R_C = kronecker( U_vert, solve(Q_mat_C-rho_matrix))  +
#     kronecker( U, solve(Q_mat_C + (cluster_size-1)*rho_matrix))
#   
#   inv_V_T = inv_A_T %*% inv_R_T %*% inv_A_T
#   inv_V_C = inv_A_C %*% inv_R_C %*% inv_A_C
#   
#   Y_mat = do.call(cbind, y01_by_cluster) # each column is a cluster
#   
#   Y_C = Y_mat[,treat_cluster ==0]
#   
#   Y_T = Y_mat[,treat_cluster ==1]
#   
#   # num_trt = number_cluster/2 #sum(treat_cluster==1) # number_cluster/2
#   # num_ct = number_cluster/2
#   
#   
#   m <- 2 * cluster_size
#   nT <- sum(treat_cluster == 1)
#   nC <- sum(treat_cluster == 0)
#   
#   # residual matrices: columns are r_i = Y_i - pi_hat
#   R_T <- Y_T - matrix(as.numeric(pi_T_hat), nrow = m, ncol = nT)
#   R_C <- Y_C - matrix(as.numeric(pi_C_hat), nrow = m, ncol = nC)
#   
#   # Sums of r_i r_i^T:
#   S_T <- R_T %*% t(R_T)   # m x m
#   S_C <- R_C %*% t(R_C)   # m x m
# 
#   
#   DV_T = t(D_T) %*% inv_V_T 
#   DV_C = t(D_C) %*% inv_V_C
#   
#   meat_T = DV_T %*% S_T  %*% t(DV_T)
#   meat_C = DV_C %*% S_C  %*% t(DV_C)
#   meat = meat_T +meat_C
#   
#   sandwich_var_beta = #number_cluster *
#     inv_sum_DVD %*% meat %*% inv_sum_DVD
#   
# }
# sqrt(diag(sandwich_var_beta) )
# pnorm(beta_new/sqrt(diag(sandwich_var_beta) ), lower.tail = F)
# 
# 
# # test stats
# 
# 
# 
# # valid_beta <- function(beta, sim_data, eps = 1e-12) {
# #   alpha <- unname(beta[grepl("^Inter:", names(beta))])
# #   beta_treat <- unname(beta["treat1"])
# #   treat_vec <- model.matrix(~ treat, data = sim_data)[, "treat1"]
# #   
# #   eta <- outer(rep(1, nrow(sim_data)), alpha) + treat_vec %o% rep(beta_treat, length(alpha))
# #   cum <- plogis(eta)                    # 2 columns for K=3
# #   # monotone CDFs:
# #   if (any(!is.finite(cum))) return(FALSE)
# #   if (any(cum[,2] <= cum[,1])) return(FALSE)
# #   
# #   # implied category probs
# #   pi_hat <- cbind(cum[,1], cum[,2]-cum[,1], 1 - cum[,2])
# #   if (any(!is.finite(pi_hat))) return(FALSE)
# #   if (any(pi_hat < eps | pi_hat > 1-eps)) return(FALSE)
# #   
# #   TRUE
# # }
# # 
# # update_with_backtracking <- function(beta_old, increment, sim_data,
# #                                      max_halves = 20, score_old = NULL, score_fn = NULL) {
# #   s <- 1
# #   for (h in 0:max_halves) {
# #     beta_try <- beta_old + s * increment
# #     ok <- valid_beta(beta_try, sim_data)
# #     # Optional: require score to decrease in infinity-norm
# #     if (ok && !is.null(score_fn)) {
# #       sc_try <- score_fn(beta_try)
# #       ok <- max(abs(sc_try)) <= max(abs(score_old))
# #     }
# #     if (ok) return(list(beta_new = beta_try, step = s))
# #     s <- s / 2
# #   }
# #   stop("Backtracking failed: step too small or model invalid")
# # }
# 
# 
# ### true probabilities.
# # > p0_C;p1_C; p2_C
# # [1] 0.762
# # [1] 0.125
# # [1] 0.113
# # > p0_T;p1_T;p2_T
# # [1] 0.69
# # [1] 0.155
# # [1] 0.155
# 
# 
# beta_increment(pi_T_hat=c(0.69,0.155),
#                pi_C_hat=c(0.762,0.125),
#                rho_matrix=rho_matrix,
#                cluster_size =cluster_size ,
#                number_cluster=number_cluster)
#   
