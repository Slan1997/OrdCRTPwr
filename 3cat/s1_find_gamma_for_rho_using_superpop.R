### confirmed. the pi in Corr(Yis,Yit) has to be marginal!
## correct sign of beta1 => ordinal dist type should follow the expectations.
### the result of this file will later be combined into "~/Composite_Endpoints/proj2/results/rho_mat_from_superpop0929.csv"
library(dplyr)
library(readr)
library(geepack)
library(ordinal)
library(bridgedist)
source("~/Composite_Endpoints/proj2/functions_clean0903.R")

epsilon_scale = 1; epsilon_sd=pi/sqrt(3)
n_cat <- 3

## find the target gamma 
gamma_mat_init = tibble(ordinal_dist_type = c("even","rare","common","bell","high01low2","U-shape"),
                        gamma0 = c(-.8, 1.3,   -2, -2.9, -0.094, -0.094), 
                        gamma1 = c( .6, 2.3, -1.1, 3.32,   3.32, .5),
                        target_pi0 = c(0.33,0.8,0.1,0.1,0.45, 0.45),
                        target_pi1 = c(0.33,0.1,0.1,0.8,0.45, 0.1),
                        target_pi2 = c(0.34,0.1,0.8,0.1,0.1, 0.45))

# use a big cluster size for super population? or just use scenario specific?
scenarios_for_super_pop = tibble(expand.grid(number_clusters = 5000, #c(10, 20, 30, 50, 100), # remove 4
                                             cluster_sizes = 100, 
                                             beta_X_list   = c(log(1.5),log(2)),
                                             latent_ICC = c( 0.01, 0.02, 0.05),   
                                             ordinal_dist_type = c("even","rare","common","bell","high01low2","U-shape")))  %>%
  full_join(gamma_mat_init,by="ordinal_dist_type") %>%
  mutate(N = number_clusters*cluster_sizes,
         m = N/2,
         scena = row_number()) %>%
  relocate(scena, .before = 1)
scenarios_for_super_pop
num_scena <- nrow(scenarios_for_super_pop)

get_pi_vec_from_gamma <- function(gamma, beta1, latent_ICC, epsilon_scale = 1) {
  cut_points <- gamma
  phi0 <- get_phi_bridge_from_ICC(ICC = latent_ICC)
  
  m <- length(cut_points)          # number of cutpoints
  K <- m + 1L                      # number of categories (0..m)
  
  p_C_vec <- numeric(K)
  p_T_vec <- numeric(K)
  
  for (k in 0:m) {
    idx <- k + 1L                  # R index (1..K)
    
    p_C_vec[idx] <- ordinal_logit_marginal_prob_bridge(
      y = k, x = 0, beta1 = beta1,
      cutpoints = cut_points,
      epsilon_scale = epsilon_scale,
      b_phi = phi0
    )
    
    p_T_vec[idx] <- ordinal_logit_marginal_prob_bridge(
      y = k, x = 1, beta1 = beta1,
      cutpoints = cut_points,
      epsilon_scale = epsilon_scale,
      b_phi = phi0
    )
  }
  
  pi_vec <- (p_C_vec + p_T_vec) / 2
  return(pi_vec)
}

# get_pi_vec_from_gamma <- function(gamma, beta1, latent_ICC, epsilon_scale = 1) {
#   cut_points <- gamma
#   phi0 <- get_phi_bridge_from_ICC(ICC = latent_ICC)
#   
#   # Control
#   p0_C <- ordinal_logit_marginal_prob_bridge(0, x = 0, beta1 = beta1,
#                                              cutpoints = cut_points,
#                                              epsilon_scale = epsilon_scale,
#                                              b_phi = phi0)
#   p1_C <- ordinal_logit_marginal_prob_bridge(1, x = 0, beta1 = beta1,
#                                              cutpoints = cut_points,
#                                              epsilon_scale = epsilon_scale,
#                                              b_phi = phi0)
#   p2_C <- ordinal_logit_marginal_prob_bridge(2, x = 0, beta1 = beta1,
#                                              cutpoints = cut_points,
#                                              epsilon_scale = epsilon_scale,
#                                              b_phi = phi0)
#   # Treatment
#   p0_T <- ordinal_logit_marginal_prob_bridge(0, x = 1, beta1 = beta1,
#                                              cutpoints = cut_points,
#                                              epsilon_scale = epsilon_scale,
#                                              b_phi = phi0)
#   p1_T <- ordinal_logit_marginal_prob_bridge(1, x = 1, beta1 = beta1,
#                                              cutpoints = cut_points,
#                                              epsilon_scale = epsilon_scale,
#                                              b_phi = phi0)
#   p2_T <- ordinal_logit_marginal_prob_bridge(2, x = 1, beta1 = beta1,
#                                              cutpoints = cut_points,
#                                              epsilon_scale = epsilon_scale,
#                                              b_phi = phi0)
#   
#   # average over arms
#   pi_vec <- c(p0_C + p0_T, p1_C + p1_T, p2_C + p2_T) / 2
#   return(pi_vec)
# }


obj_gamma <- function(gamma, target_pi_vec, beta1, latent_ICC, epsilon_scale = 1) {
  pi_vec <- get_pi_vec_from_gamma(gamma, beta1 = beta1,
                                  latent_ICC = latent_ICC,
                                  epsilon_scale = epsilon_scale)
  sum((pi_vec - target_pi_vec)^2)
}



epsilon_scale <- 1
n_cat <- 3

num_scena <- nrow(scenarios_for_super_pop)

calibrated_gamma <- vector("list", num_scena)

for (scena in 1:num_scena) {
  this_row <- scenarios_for_super_pop[scena, ]
  
  # starting values
  gamma_start <- c(this_row$gamma0, this_row$gamma1)
  
  target_pi_vec <- this_row %>%
    dplyr::select(contains("target_pi")) %>%
    as.numeric()
  
  beta1      <- this_row$beta_X_list
  latent_ICC <- this_row$latent_ICC
  
  fit <- optim(
    par  = gamma_start,
    fn   = obj_gamma,
    target_pi_vec = target_pi_vec,
    beta1 = beta1,
    latent_ICC = latent_ICC,
    epsilon_scale = epsilon_scale,
    method = "BFGS"
  )
  
  gamma_hat <- fit$par
  calibrated_gamma[[scena]] <- list(
    scena = this_row$scena,
    gamma0 = gamma_hat[1],
    gamma1 = gamma_hat[2],
    obj = fit$value,
    converged = (fit$convergence == 0),
    pi_vec = get_pi_vec_from_gamma(gamma_hat, beta1, latent_ICC, epsilon_scale),
    target_pi_vec = target_pi_vec
  )
}

gamma_res <- purrr::map_dfr(calibrated_gamma, ~{
  tibble(
    scena = .x$scena,
    gamma0_hat = .x$gamma0,
    gamma1_hat = .x$gamma1,
    obj = .x$obj, # sum((pi_vec - target_pi_vec)^2)
    converged = .x$converged,
    pi0_hat = .x$pi_vec[1],
    pi1_hat = .x$pi_vec[2],
    pi2_hat = .x$pi_vec[3],
    pi0_target = .x$target_pi_vec[1],
    pi1_target = .x$target_pi_vec[2],
    pi2_target = .x$target_pi_vec[3]
  )
})
gamma_res

gamma_res %>%
  mutate(
    err0 = pi0_hat - pi0_target,
    err1 = pi1_hat - pi1_target,
    err2 = pi2_hat - pi2_target
  )
gamma_res



scenarios_for_super_pop_updated <- scenarios_for_super_pop %>%
  dplyr::select(-contains("target")) %>%
  left_join(
    gamma_res %>% 
      dplyr::select(scena, gamma0_hat, gamma1_hat),
    by = "scena"
  ) %>%
  mutate(
    gamma0 = gamma0_hat,
    gamma1 = gamma1_hat
  ) %>%
  dplyr::select(-gamma0_hat, -gamma1_hat)

scenarios_for_super_pop_updated %>% write_csv("scenarios_for_superpop_updated_beta1.5and2.csv")
