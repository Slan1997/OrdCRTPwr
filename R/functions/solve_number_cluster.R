## contains:
# solve_number_cluster()
# solve_clusters_binary()


expit <- function(x) 1/(1+exp(-x))
logit <- function(p) log(p/(1-p))

pi_T_from_PO <- function(pi_C_hat, theta) {
  K <- length(pi_C_hat)                 # K = #thresholds = #cats - 1
  vartheta_C <- cumsum(pi_C_hat)[1:K]   # cum probs up to k=1..K
  vartheta_T <- expit(logit(vartheta_C) + theta)
  
  pi_T_hat <- numeric(K)
  pi_T_hat[1] <- vartheta_T[1]
  if (K >= 2) pi_T_hat[2:K] <- diff(vartheta_T)
  
  # sanity checks
  if (any(pi_T_hat < -1e-12)) warning("pi_T_hat has negative entries; check inputs/ordering.")
  if (abs(sum(pi_T_hat) - vartheta_T[K]) > 1e-8) warning("pi_T_hat sum mismatch; check computations.")
  pi_T_hat
}

solve_number_cluster <- function(pi_C_hat,
                                 theta,
                                 cluster_size,
                                 rho_matrix,
                                 alpha = 0.05,
                                 target_power = 0.80,
                                 n_cat_index = NULL,
                                 round_even = F) {
  # pi_C_hat: length K (#cats - 1), i.e., omitting last category
  # theta: log OR
  # cluster_size: J
  # rho_matrix: K x K (your estimated rho for the binary-indicator vector)
  # n_cat_index: index of theta in beta vector; default assumes last element (K+1)
  #             If your theta is at position 'n_cat' in your code, pass that in.
  
  K <- length(pi_C_hat)
  
  if (is.null(n_cat_index)) n_cat_index <- K + 1  # common: beta=(gamma1..gammaK, theta)
  
  # 1) Get treatment distribution (first K probs)
  pi_T_hat <- pi_T_from_PO(pi_C_hat, theta)
  
  
  # print(pi_C_hat)
  # print(pi_T_hat)
  # 
  # 2) Build pieces (using YOUR helper functions)
  a_mat_T <- make_a_mat(pi_T_hat)                 # K x K
  a_mat_C <- make_a_mat(pi_C_hat)                 # K x K
  
  Q_mat_T <- make_Q_mat_fast(pi_T_hat)            # K x K
  Q_mat_C <- make_Q_mat_fast(pi_C_hat)            # K x K
  
  d_mat_T <- make_d_mat(pi_T_hat, treat_last_col = TRUE)    # K x (K+1)
  d_mat_C <- make_d_mat(pi_C_hat, treat_last_col = FALSE)   # K x (K+1)
  
  inv_part_T <- solve(Q_mat_T + (cluster_size - 1) * rho_matrix)
  inv_part_C <- solve(Q_mat_C + (cluster_size - 1) * rho_matrix)
  
  dvd_T <- (t(d_mat_T) %*% a_mat_T %*% inv_part_T %*% a_mat_T) %*% d_mat_T
  dvd_C <- (t(d_mat_C) %*% a_mat_C %*% inv_part_C %*% a_mat_C) %*% d_mat_C
  
  # 3) Var(theta_hat) = (2/(J*N)) * [solve(dvd_T + dvd_C)]_{theta,theta}
  M <- solve(dvd_T + dvd_C)
  v0 <- (2 / cluster_size) * M[n_cat_index, n_cat_index]   # so Var = v0 / N
  
  if (v0 <= 0) stop("Computed v0 <= 0. Check rho_matrix, Q/A matrices, and indexing.")
  
  # 4) Back-solve N from power equation:
  # power = Phi(|theta|/sqrt(v0/N) - z_{1-alpha/2})
  # set power = target_power -> |theta|/sqrt(v0/N) = z_{1-alpha/2} + z_{target_power}
  z_alpha_2 <- qnorm(1 - alpha/2)
  z_pow     <- qnorm(target_power)
  z_sum     <- z_alpha_2 + z_pow
  
  N_req <- v0 * (z_sum^2) / (theta^2)
  
  # round up to integer (and even if balanced 1:1)
  N_int <- ceiling(N_req)
  if (round_even && (N_int %% 2 == 1)) N_int <- N_int + 1
  
  # optional: return also achieved power for the rounded N
  var_theta_hat <- v0 / N_int
  achieved_power <- pnorm(abs(theta)/sqrt(var_theta_hat) - z_alpha_2)
  
  list(
    N_required_continuous = N_req,
    N_required_integer = N_int,
    achieved_power = achieved_power,
    v0 = v0,
    pi_T_hat = pi_T_hat,
    pi_C_hat = pi_C_hat
  )
}


solve_clusters_binary <- function(P0, theta_R, rho, cluster_size,
                                  alpha = 0.05, target_power = 0.80,
                                  allocate_prop_C = 0.5,
                                  round_even = F) {
  OR <- exp(theta_R)
  P1 <- OR * P0 / (1 - P0 + OR * P0)
  # print(P0)
  # print(P1)
  DE <- 1 + (cluster_size - 1) * rho
  var_term <- 1/(allocate_prop_C * P0 * (1 - P0)) +
    1/((1 - allocate_prop_C) * P1 * (1 - P1))
  
  z_alpha_2 <- qnorm(1 - alpha/2)
  z_pow     <- qnorm(target_power)
  z_sum     <- z_alpha_2 + z_pow
  
  m_cont <- (DE * var_term / (cluster_size * theta_R^2)) * (z_sum^2)
  
  m_int <- ceiling(m_cont)
  if (round_even && (m_int %% 2 == 1)) m_int <- m_int + 1
  
  final_term1 <- sqrt(cluster_size * m_int * theta_R^2 / (DE * var_term))
  achieved_power <- pnorm(final_term1 - z_alpha_2)
  
  list(
    N_required_continuous = m_cont,
    N_required_integer = m_int,
    achieved_power = achieved_power,
    v0 = (DE * var_term) / cluster_size,   # so Var(theta_hat)= v0 / m
    pi_T_hat = P1 ,                         # treatment prob
    pi_C_hat = P0                      # control prob
  )
}
