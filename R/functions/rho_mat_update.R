# ============================================================
# Estimate the exchangeable within-cluster correlation matrix
# for ordinal category indicators using standardized residuals.
#
# For a K-category ordinal outcome, the function estimates the
# (K - 1) x (K - 1) correlation matrix for the first K - 1
# category indicators.
#
# Method-of-moments estimator:
#
#   rho_hat =
#     sum_i sum_{t > s} e_is e_it^T
#     --------------------------------
#     [sum_i T_i(T_i - 1) / 2] - P
#
# where e_it is the standardized residual vector for subject t
# in cluster i, and P is the number of marginal mean parameters.
#
# Assumptions:
#   - sim_data has columns: cluster, treat, Y
#   - Y is coded as 0, 1, ..., K - 1
#   - treat is coded as 0/1 or factor with levels 0 and 1
#   - beta = c(gamma_0, ..., gamma_{K-2}, beta_treat)
#   - cumulative logit model:
#       logit Pr(Y <= j | treat) = gamma_j + beta_treat * treat
# ============================================================

### inv_logit has been defined in ordinal_probabilities.R
# inv_logit <- function(x, scale=1){
#   1 / (1 + exp(-x/scale))
# }


rho_mat_update <- function(beta,
                           sim_data,
                           useP = TRUE,
                           P = NULL,
                           K = NULL,
                           symmetrize = TRUE) {
  
  # ----------------------------
  # Basic input checks
  # ----------------------------
  required_cols <- c("cluster", "treat", "Y")
  missing_cols <- setdiff(required_cols, names(sim_data))
  
  if (length(missing_cols) > 0) {
    stop(
      "sim_data must contain columns: ",
      paste(required_cols, collapse = ", "),
      ". Missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (!is.numeric(beta) || length(beta) < 3) {
    stop("beta must be a numeric vector with at least 3 elements.")
  }
  
  # Determine number of ordinal categories
  if (is.null(K)) {
    K <- length(sort(unique(sim_data$Y)))
  }
  
  if (!is.numeric(K) || length(K) != 1 || K < 3) {
    stop("K must be a scalar integer >= 3.")
  }
  
  K <- as.integer(K)
  n_thresh <- K - 1
  
  # beta should contain K - 1 thresholds plus one treatment effect
  if (length(beta) != n_thresh + 1) {
    stop(sprintf(
      "beta length should be K = (K - 1 thresholds) + 1 treatment effect. Got %d, expected %d.",
      length(beta), n_thresh + 1
    ))
  }
  
  if (is.null(P)) {
    P <- length(beta)
  }
  
  if (!is.numeric(P) || length(P) != 1 || !is.finite(P) || P < 0) {
    stop("P must be a finite nonnegative numeric scalar.")
  }
  
  # ----------------------------
  # Parse model parameters
  # ----------------------------
  gammas_est <- sort(as.numeric(beta[1:n_thresh]))
  beta_treat_est <- as.numeric(beta[n_thresh + 1])
  
  # Make sure treat is matched consistently to 0/1
  sim_data <- sim_data %>%
    dplyr::mutate(
      treat = as.factor(as.character(treat)),
      Y_num = as.integer(as.character(Y))
    )
  
  if (!all(c("0", "1") %in% levels(sim_data$treat))) {
    stop("treat must contain levels coded as 0 and 1.")
  }
  
  if (any(is.na(sim_data$Y_num))) {
    stop("Y must be coded as numeric categories 0, 1, ..., K - 1.")
  }
  
  if (!all(sim_data$Y_num %in% 0:(K - 1))) {
    stop("Y must be coded as 0, 1, ..., K - 1.")
  }
  
  # ----------------------------
  # Estimated category probabilities by treatment group
  # ----------------------------
  treat01 <- c(0, 1)
  
  # cum_le: (K - 1) x 2 matrix
  # column 1 = control, column 2 = treatment
  cum_le <- sapply(
    treat01,
    function(g) inv_logit(gammas_est + beta_treat_est * g)
  )
  
  # Convert cumulative probabilities into category probabilities
  # pi_mat: K x 2 matrix
  pi_mat <- rbind(
    cum_le[1, , drop = FALSE],
    if (K > 2) diff(cum_le) else NULL,
    1 - cum_le[n_thresh, , drop = FALSE]
  )
  
  # Use only the first K - 1 category indicators
  pi_use <- pi_mat[1:n_thresh, , drop = FALSE]
  
  # Numerical guard against invalid probabilities
  if (any(!is.finite(pi_use)) || any(pi_use <= 0) || any(pi_use >= 1)) {
    stop("Estimated probabilities contain non-finite values or values outside (0, 1).")
  }
  
  # Treatment-specific probability lookup table
  p_i_C_and_T_est <- tibble::tibble(
    treat = factor(c("0", "1"), levels = levels(sim_data$treat))
  )
  
  for (k in 1:n_thresh) {
    p_i_C_and_T_est[[paste0("pi_", k)]] <- pi_use[k, ]
  }
  
  # ----------------------------
  # Create category indicators and standardized residuals
  # ----------------------------
  sim_data_bin_est <- sim_data %>%
    dplyr::transmute(
      cluster,
      treat,
      Y_num
    )
  
  for (k in 1:n_thresh) {
    sim_data_bin_est[[paste0("bin_Y", k)]] <-
      as.integer(sim_data_bin_est$Y_num == (k - 1))
  }
  
  sim_data_bin_est <- sim_data_bin_est %>%
    dplyr::left_join(p_i_C_and_T_est, by = "treat")
  
  e_mat <- sim_data_bin_est
  
  for (k in 1:n_thresh) {
    pi_name <- paste0("pi_", k)
    bin_name <- paste0("bin_Y", k)
    e_name <- paste0("e", k)
    
    e_mat[[e_name]] <- (e_mat[[bin_name]] - e_mat[[pi_name]]) /
      sqrt(e_mat[[pi_name]] * (1 - e_mat[[pi_name]]))
  }
  
  e_cols <- paste0("e", 1:n_thresh)
  
  # ----------------------------
  # Method-of-moments correlation estimator
  # ----------------------------
  rho_parts <- e_mat %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      Ti = dplyr::n(),
      
      # Numerator contribution:
      #   sum_{t > s} e_is e_it^T
      num = list({
        E <- as.matrix(dplyr::pick(dplyr::all_of(e_cols)))
        Ti_current <- nrow(E)
        
        S <- matrix(0, nrow = n_thresh, ncol = n_thresh)
        
        if (Ti_current >= 2) {
          for (s in 1:(Ti_current - 1)) {
            for (t in (s + 1):Ti_current) {
              S <- S + tcrossprod(E[s, ], E[t, ])
            }
          }
        }
        
        S
      }),
      
      # Number of unordered within-cluster pairs:
      #   T_i(T_i - 1) / 2
      n_pairs = dplyr::n() * (dplyr::n() - 1) / 2,
      .groups = "drop"
    )
  
  # Sum numerator matrices across clusters
  num_total <- Reduce("+", rho_parts$num)
  
  # Denominator:
  #   sum_i T_i(T_i - 1) / 2 - P
  denom_total <- sum(rho_parts$n_pairs) - if (useP) P else 0
  
  # Fallback if finite-sample correction gives invalid denominator
  if (!is.finite(denom_total) || denom_total <= 0) {
    denom_total <- sum(rho_parts$n_pairs)
  }
  
  rho_est <- num_total / denom_total
  
  # For an exchangeable working correlation matrix, enforce symmetry
  if (symmetrize) {
    rho_est <- (rho_est + t(rho_est)) / 2
  }
  
  dimnames(rho_est) <- list(
    paste0("Y=", 0:(n_thresh - 1)),
    paste0("Y=", 0:(n_thresh - 1))
  )
  
  rho_est
}