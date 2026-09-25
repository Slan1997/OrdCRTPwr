sandwich_var_ord <- function(pi_T_hat,
                             pi_C_hat,
                             rho_matrix,
                             cluster_size,
                             number_cluster,
                             y_by_cluster,
                             treat_cluster) {
  
  # coerce to plain vectors (just make sure)
  pi_T_hat <- as.numeric(pi_T_hat)
  pi_C_hat <- as.numeric(pi_C_hat)
  
  K <- length(pi_T_hat)
  if (length(pi_C_hat) != K) stop("pi_T_hat and pi_C_hat must have same length K.")
  if (!all(dim(rho_matrix) == c(K, K))) stop("rho_matrix must be K x K.")
  if (length(y_by_cluster) != number_cluster) stop("y_by_cluster length must equal number_cluster.")
  if (length(treat_cluster) != number_cluster) stop("treat_cluster length must equal number_cluster.")
  
  # ---- build model-based pieces ----
  a_mat_T <- make_a_mat(pi_T_hat)          # K x K
  a_mat_C <- make_a_mat(pi_C_hat)
  
  Q_mat_T <- make_Q_mat_fast(pi_T_hat)     # K x K
  Q_mat_C <- make_Q_mat_fast(pi_C_hat)
  
  d_mat_T <- make_d_mat(pi_T_hat, treat_last_col = TRUE)   # K x (K+1)
  d_mat_C <- make_d_mat(pi_C_hat, treat_last_col = FALSE)  # K x (K+1), last col 0
  
  inv_part_T <- solve(Q_mat_T + (cluster_size - 1) * rho_matrix)
  inv_part_C <- solve(Q_mat_C + (cluster_size - 1) * rho_matrix)
  
  dv_T  <- t(d_mat_T) %*% a_mat_T %*% inv_part_T %*% a_mat_T    # (K+1) x K
  dv_C  <- t(d_mat_C) %*% a_mat_C %*% inv_part_C %*% a_mat_C    # (K+1) x K
  dvd_T <- dv_T %*% d_mat_T                                      # (K+1) x (K+1)
  dvd_C <- dv_C %*% d_mat_C                                      # (K+1) x (K+1)
  
  inv_sum_DVD <- 2 / cluster_size / number_cluster * solve(dvd_T + dvd_C)
  # this is the "bread" (model-based var part)
  
  # ---- within-cluster precision (inv_V) ----
  oneJ <- matrix(1, nrow = cluster_size, ncol = 1)
  
  D_T <- kronecker(oneJ, d_mat_T)                # (J*K) x (K+1)
  D_C <- kronecker(oneJ, d_mat_C)                # (J*K) x (K+1)
  
  inv_A_T <- kronecker(diag(cluster_size), a_mat_T)  # (J*K) x (J*K)
  inv_A_C <- kronecker(diag(cluster_size), a_mat_C)
  
  U      <- (1 / cluster_size) * matrix(1, cluster_size, cluster_size)
  U_vert <- diag(cluster_size) - U
  
  inv_R_T <- kronecker(U_vert, solve(Q_mat_T - rho_matrix)) +
    kronecker(U,      solve(Q_mat_T + (cluster_size - 1) * rho_matrix))
  
  inv_R_C <- kronecker(U_vert, solve(Q_mat_C - rho_matrix)) +
    kronecker(U,      solve(Q_mat_C + (cluster_size - 1) * rho_matrix))
  
  inv_V_T <- inv_A_T %*% inv_R_T %*% inv_A_T      # (J*K) x (J*K)
  inv_V_C <- inv_A_C %*% inv_R_C %*% inv_A_C
  
  # ---- assemble data matrix ----
  Y_mat <- do.call(cbind, y_by_cluster)  # (J*K) x (#clusters)
  if (!all(dim(Y_mat) == c(cluster_size * K, number_cluster))) {
    stop("Each y_by_cluster element must be length cluster_size*K and list length = number_cluster.")
  }
  
  Y_T <- Y_mat[, treat_cluster == 1, drop = FALSE]
  Y_C <- Y_mat[, treat_cluster == 0, drop = FALSE]
  
  nT <- ncol(Y_T)
  nC <- ncol(Y_C)
  m  <- cluster_size * K   # length of stacked vector per cluster
  
  # residual matrices: columns are r_i = Y_i - mu_hat
  # mu_hat is the stacked mean vector for a cluster:
  # (1_J ⊗ I_K) pi_hat  == rep(pi_hat, each=J) with subject-major stacking
  mu_T <- rep(pi_T_hat, times = cluster_size)   # length m (subject-major)
  mu_C <- rep(pi_C_hat, times = cluster_size)
  
  R_T <- if (nT > 0) (Y_T - matrix(mu_T, nrow = m, ncol = nT)) else matrix(0, m, 0)
  R_C <- if (nC > 0) (Y_C - matrix(mu_C, nrow = m, ncol = nC)) else matrix(0, m, 0)
  
  S_T <- if (nT > 0) (R_T %*% t(R_T)) else matrix(0, m, m)
  S_C <- if (nC > 0) (R_C %*% t(R_C)) else matrix(0, m, m)
  
  DV_T <- t(D_T) %*% inv_V_T    # (K+1) x m
  DV_C <- t(D_C) %*% inv_V_C
  
  meat_T <- DV_T %*% S_T %*% t(DV_T)   # (K+1) x (K+1)
  meat_C <- DV_C %*% S_C %*% t(DV_C)
  meat <- meat_T + meat_C
  
  sandwich_var_beta <- inv_sum_DVD %*% meat %*% inv_sum_DVD
  
  list(
    var_beta = inv_sum_DVD,
    sandwich_var_beta = sandwich_var_beta,
    meat = meat
  )
}

sandwich_var = function(pi_T_hat,
                        pi_C_hat,
                        rho_matrix,
                        cluster_size ,
                        number_cluster){
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

  Y_mat = do.call(cbind, y01_by_cluster) # each column is a cluster

  Y_C = Y_mat[,treat_cluster ==0]

  Y_T = Y_mat[,treat_cluster ==1]

  # num_trt = number_cluster/2 #sum(treat_cluster==1) # number_cluster/2
  # num_ct = number_cluster/2


  m <- 2 * cluster_size
  nT <- sum(treat_cluster == 1)
  nC <- sum(treat_cluster == 0)

  # residual matrices: columns are r_i = Y_i - pi_hat
  R_T <- Y_T - matrix(as.numeric(pi_T_hat), nrow = m, ncol = nT)
  R_C <- Y_C - matrix(as.numeric(pi_C_hat), nrow = m, ncol = nC)

  # Sums of r_i r_i^T:
  S_T <- R_T %*% t(R_T)   # m x m
  S_C <- R_C %*% t(R_C)   # m x m


  DV_T = t(D_T) %*% inv_V_T
  DV_C = t(D_C) %*% inv_V_C

  meat_T = DV_T %*% S_T  %*% t(DV_T)
  meat_C = DV_C %*% S_C  %*% t(DV_C)
  meat = meat_T +meat_C
 
  sandwich_var_beta = #number_cluster *
    inv_sum_DVD %*% meat %*% inv_sum_DVD
  
  return(list(var_beta = inv_sum_DVD,
              sandwich_var_beta = sandwich_var_beta))

}
#sqrt(diag(sandwich_var_beta) )
#pnorm(beta_new/sqrt(diag(sandwich_var_beta) ), lower.tail = F)
