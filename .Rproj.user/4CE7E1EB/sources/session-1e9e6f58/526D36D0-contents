# helpers
make_a_mat <- function(pi_vec) {
  K <- length(pi_vec)
  a_vec <- (pi_vec * (1 - pi_vec))^(-1/2)
  diag(a_vec, nrow = K, ncol = K)
}

# make_Q_mat <- function(pi_vec) {
#   K <- length(pi_vec)
#   a_vec <- (pi_vec * (1 - pi_vec))^(-1/2)
#   
#   Q <- diag(1, K)
#   
#   if (K > 1) {
#     for (k in 1:(K - 1)) {
#       for (l in (k + 1):K) {
#         r_kl <- -pi_vec[k] * pi_vec[l] * a_vec[k] * a_vec[l]
#         Q[k, l] <- r_kl
#         Q[l, k] <- r_kl
#       }
#     }
#   }
#   Q
# }

make_Q_mat_fast <- function(pi_vec) {
  a <- (pi_vec * (1 - pi_vec))^(-1/2)
  
  r <- -outer(pi_vec, pi_vec) * outer(a, a)
  
  diag(r) <- 1
  r
}

make_d_mat <- function(pi_vec, treat_last_col = TRUE) {
  K <- length(pi_vec)
  
  vartheta <- cumsum(pi_vec)
  v <- vartheta * (1 - vartheta)
  
  d <- matrix(0, nrow = K, ncol = K + 1)
  
  for (k in 1:K) {
    d[k, k] <- v[k]
    
    if (k >= 2)
      d[k, k - 1] <- -v[k - 1]
    
    if (treat_last_col) {
      prev <- if (k >= 2) v[k - 1] else 0
      d[k, K + 1] <- v[k] - prev
    }
  }
  d
}