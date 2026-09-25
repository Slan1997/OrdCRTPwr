# Convert cumulative-logit coefficients into estimated category probabilities for control and treatment groups
get_pi_from_coef <- function(coef, n_cat) {
  # coef = (gamma_1,...,gamma_{K-1}, beta_treat)
  # n_cat = total number of ordinal categories (K)
  
  #inv_logit <- function(x) 1 / (1 + exp(-x)) # already defined in ordinal_probabilities.R
  
  K <- n_cat
  if (K < 3) stop("Need at least 3 categories.")
  
  n_thresh <- K - 1
  
  if (length(coef) != n_thresh + 1) {
    stop("coef must have length (K-1 thresholds + 1 treat effect).")
  }
  
  gammas <- sort(as.numeric(coef[1:n_thresh]))
  beta_trt <- as.numeric(coef[n_thresh + 1])
  
  treat01 <- c(0, 1)
  
  # (K-1) x 2 cumulative probs
  cum_le <- sapply(
    treat01,
    function(g) inv_logit(gammas + beta_trt * g)
  )
  
  # K x 2 category probs
  pi_mat <- rbind(
    cum_le[1, , drop = FALSE],
    if (K > 2) diff(cum_le) else NULL,
    1 - cum_le[n_thresh, , drop = FALSE]
  )
  
  # # Numerical guard
  # pi_mat <- pmin(pmax(pi_mat, 1e-12), 1 - 1e-12)
  
  # First K-1 categories (drop last)
  pi_use <- pi_mat[1:n_thresh, , drop = FALSE]
  
  list(
    pi_mat   = pi_mat,       # K x 2
    pi_use   = pi_use,       # (K-1) x 2
    pi_C_hat = pi_use[, 1],  # control
    pi_T_hat = pi_use[, 2],  # treatment
    cum_le   = cum_le
  )
}
