# ============================================================
# 00_create_superpop_scenarios.R
#
# Create the superpopulation scenario grid for the 3-category
# ordinal CRT simulation study.
#
# Output:
#   results/3cat/superpop/scenarios_for_superpop_3cat.csv
#
# This file defines:
#   - ordinal distribution shapes
#   - cutpoints used in the latent proportional odds model
#   - treatment effect settings
#   - latent ICC settings
#   - large superpopulation design used to estimate rho matrices
#     and scalar ICC summaries
# ============================================================


# ----------------------------
# Packages
# ----------------------------
library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(bridgedist)


# ----------------------------
# Source helper functions
# ----------------------------
source("R/functions/bridge_helpers.R")
source("R/functions/ordinal_probabilities.R")


# ----------------------------
# Output directory
# ----------------------------
out_dir <- "results/3cat/superpop"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# ----------------------------
# Global settings
# ----------------------------
epsilon_scale <- 1

# Large superpopulation settings used to estimate dependence parameters.
# These are intentionally much larger than the empirical-power simulation
# trial sizes so that the estimated rho/ICC values are treated as "truth".
number_clusters_superpop <- 5000
cluster_size_superpop    <- 100

# Conditional treatment effects.
# Main manuscript figures use OR = 1.5; OR = 2 may be kept for supplementary
# or sensitivity analyses if needed.
OR_values <- c(1.5, 2)
beta_X_values <- log(OR_values)

# Latent ICC values.
latent_ICC_values <- c(0.01, 0.02, 0.05)


# ----------------------------
# Three-category distribution settings
# ----------------------------
# Categories are coded as 0, 1, 2 in the simulation code.
#
# target_pi0, target_pi1, target_pi2 describe the approximate marginal
# distribution in the control arm when there is no treatment effect and
# under the bridge random-intercept model.
#
# gamma0 and gamma1 are approximate cutpoints used in the latent model.
# These were chosen to produce the desired ordinal distribution shapes.
gamma_mat_init <- tibble(
  ordinal_dist_type = c(
    "bell",
    "common",
    "rare",
    "U-shape",
    "high01low2",
    "even"
  ),
  ordinal_dist_label = c(
    "Bell",
    "Common",
    "Rare",
    "U-shape",
    "Rare-top",
    "Even"
  ),
  gamma0 = c(
    -2.9,
    -2.0,
    1.3,
    -0.094,
    -0.094,
    -0.8
  ),
  gamma1 = c(
    3.32,
    -1.1,
    2.3,
    0.5,
    3.32,
    0.6
  ),
  target_pi0 = c(
    0.10,
    0.10,
    0.80,
    0.45,
    0.45,
    0.33
  ),
  target_pi1 = c(
    0.80,
    0.10,
    0.10,
    0.10,
    0.45,
    0.33
  ),
  target_pi2 = c(
    0.10,
    0.80,
    0.10,
    0.45,
    0.10,
    0.34
  )
)


# ----------------------------
# Helper: compute marginal probabilities by bridge integration
# ----------------------------
get_pi_vec_from_gamma <- function(gamma, beta1, latent_ICC,
                                  epsilon_scale = 1) {
  
  cut_points <- gamma
  phi0 <- get_phi_bridge_from_ICC(ICC = latent_ICC)
  
  n_cut <- length(cut_points)
  n_cat <- n_cut + 1L
  
  p_C_vec <- numeric(n_cat)
  p_T_vec <- numeric(n_cat)
  
  for (y in 0:n_cut) {
    idx <- y + 1L
    
    p_C_vec[idx] <- ordinal_logit_marginal_prob_bridge(
      y = y,
      x = 0,
      beta1 = beta1,
      cutpoints = cut_points,
      epsilon_scale = epsilon_scale,
      b_phi = phi0
    )
    
    p_T_vec[idx] <- ordinal_logit_marginal_prob_bridge(
      y = y,
      x = 1,
      beta1 = beta1,
      cutpoints = cut_points,
      epsilon_scale = epsilon_scale,
      b_phi = phi0
    )
  }
  
  list(
    pi_C = p_C_vec,
    pi_T = p_T_vec,
    pi_pooled = (p_C_vec + p_T_vec) / 2,
    phi = phi0,
    marginal_log_OR = phi0 * beta1
  )
}


# ----------------------------
# Create scenario grid
# ----------------------------
scenarios_for_superpop <- tidyr::expand_grid(
  number_clusters = number_clusters_superpop,
  cluster_sizes   = cluster_size_superpop,
  beta_X_list     = beta_X_values,
  latent_ICC      = latent_ICC_values,
  ordinal_dist_type = gamma_mat_init$ordinal_dist_type
) %>%
  left_join(gamma_mat_init, by = "ordinal_dist_type") %>%
  mutate(
    OR = exp(beta_X_list),
    phi = get_phi_bridge_from_ICC(latent_ICC),
    marginal_log_OR = phi * beta_X_list,
    N = number_clusters * cluster_sizes,
    m = N / 2
  ) %>%
  arrange(
    beta_X_list,
    latent_ICC,
    match(
      ordinal_dist_type,
      c("bell", "common", "rare", "U-shape", "high01low2", "even")
    )
  ) %>%
  mutate(scena = row_number()) %>%
  relocate(scena, .before = 1)


# ----------------------------
# Add implied marginal probabilities
# ----------------------------
pi_info <- lapply(seq_len(nrow(scenarios_for_superpop)), function(i) {
  
  gamma_i <- c(
    scenarios_for_superpop$gamma0[i],
    scenarios_for_superpop$gamma1[i]
  )
  
  out <- get_pi_vec_from_gamma(
    gamma = gamma_i,
    beta1 = scenarios_for_superpop$beta_X_list[i],
    latent_ICC = scenarios_for_superpop$latent_ICC[i],
    epsilon_scale = epsilon_scale
  )
  
  tibble(
    scena = scenarios_for_superpop$scena[i],
    
    pi_0_C = out$pi_C[1],
    pi_1_C = out$pi_C[2],
    pi_2_C = out$pi_C[3],
    
    pi_0_T = out$pi_T[1],
    pi_1_T = out$pi_T[2],
    pi_2_T = out$pi_T[3],
    
    pi_0 = out$pi_pooled[1],
    pi_1 = out$pi_pooled[2],
    pi_2 = out$pi_pooled[3]
  )
}) %>%
  bind_rows()


scenarios_for_superpop <- scenarios_for_superpop %>%
  left_join(pi_info, by = "scena")


# ----------------------------
# Basic checks
# ----------------------------
check_probs <- scenarios_for_superpop %>%
  mutate(
    sum_pi_C = pi_0_C + pi_1_C + pi_2_C,
    sum_pi_T = pi_0_T + pi_1_T + pi_2_T,
    sum_pi_pooled = pi_0 + pi_1 + pi_2
  )

stopifnot(all(abs(check_probs$sum_pi_C - 1) < 1e-8))
stopifnot(all(abs(check_probs$sum_pi_T - 1) < 1e-8))
stopifnot(all(abs(check_probs$sum_pi_pooled - 1) < 1e-8))

stopifnot(all(scenarios_for_superpop$latent_ICC > 0))
stopifnot(all(scenarios_for_superpop$latent_ICC < 1))
stopifnot(all(scenarios_for_superpop$cluster_sizes > 0))
stopifnot(all(scenarios_for_superpop$number_clusters > 0))


# ----------------------------
# Save
# ----------------------------
out_file <- file.path(out_dir, "scenarios_for_superpop_3cat.csv")

write_csv(scenarios_for_superpop, out_file)

message("Saved 3-category superpopulation scenarios to: ", out_file)
message("Number of scenarios: ", nrow(scenarios_for_superpop))


# ----------------------------
# Optional quick view
# ----------------------------
print(
  scenarios_for_superpop %>%
    select(
      scena,
      ordinal_dist_type,
      OR,
      latent_ICC,
      gamma0,
      gamma1,
      target_pi0,
      target_pi1,
      target_pi2,
      pi_0_C,
      pi_1_C,
      pi_2_C,
      pi_0_T,
      pi_1_T,
      pi_2_T
    )
)