### update: 02/03/26: add u-shape for ordinal_dist_type
## target distribution for 4-categorical ordinal outcome
# 1. even: (.25,.25,.25,.25)
# 2. rare: (.6,.2,.1,.1)
# 3. common: (.1,.1,.2,.6)
# 4. bell: (.1,.4,.4,.1)
# 5. U-shape: (.4,.1,.1,.4)


library(dplyr)
library(readr)
library(geepack)
library(ordinal)
library(bridgedist)
source("~/Composite_Endpoints/proj2/functions_clean0903.R")

epsilon_scale = 1; epsilon_sd=pi/sqrt(3)
n_cat <- 4

## find the target gamma 
gamma_mat_init = tibble(ordinal_dist_type = c("even","rare","common","bell",
                                              "U-shape","bell2"),
                        gamma0 = c(-.8,   1,   -1, -1,    0, -1), 
                        gamma1 = c( 0,  1.3,  -.8,  2,   .1 , 2),
                        gamma2 = c( .1, 1.4,    0,  3,   .5,  3),
                        target_pi0 = c(0.25,0.6,0.1,0.1,0.4, 0.1),
                        target_pi1 = c(0.25,0.2,0.1,0.4,0.1, 0.6),
                        target_pi2 = c(0.25,0.1,0.2,0.4,0.1, 0.2),
                        target_pi3 = c(0.25,0.1,0.6,0.1,0.4, 0.1))

# use a big cluster size for super population? or just use scenario specific?
scenarios_for_super_pop = tibble(expand.grid(number_clusters = 5000, #c(10, 20, 30, 50, 100), # remove 4
                                             cluster_sizes = 100, 
                                             beta_X_list   = c(log(1.5)#,log(2)
                                             ),
                                             latent_ICC = c( 0.01, 0.02, 0.05),   
                                             ordinal_dist_type = c("even","rare","common","bell",#"high01low2",
                                                                   "U-shape","bell2")))  %>%
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


obj_gamma <- function(gamma, target_pi_vec, beta1, latent_ICC, epsilon_scale = 1) {
  pi_vec <- get_pi_vec_from_gamma(gamma, beta1 = beta1,
                                  latent_ICC = latent_ICC,
                                  epsilon_scale = epsilon_scale)
  sum((pi_vec - target_pi_vec)^2)
}



epsilon_scale <- 1
n_cat <- 4

num_scena <- nrow(scenarios_for_super_pop)

calibrated_gamma <- vector("list", num_scena)

for (scena in 1:num_scena) {
  this_row <- scenarios_for_super_pop[scena, ]
  
  # starting values
  gamma_start <- this_row %>% dplyr::select(contains("gamma")) %>% as.numeric
  
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
  gamma_list <- as.list(gamma_hat)
  names(gamma_list) <- paste0("gamma", seq_along(gamma_hat) - 1L)  # gamma0, gamma1, ...
  
  calibrated_gamma[[scena]] <- c(list(
    scena = this_row$scena,
    obj = fit$value,
    converged = (fit$convergence == 0),
    pi_vec = get_pi_vec_from_gamma(gamma_hat, beta1, latent_ICC, epsilon_scale),
    target_pi_vec = target_pi_vec
  ), gamma_list)
}

gamma_res <- purrr::map_dfr(calibrated_gamma, function(x) {
  
  # --- pull gamma0, gamma1, ... dynamically ---
  gamma_names <- grep("^gamma\\d+$", names(x), value = TRUE)
  gamma_names <- gamma_names[order(as.integer(sub("^gamma", "", gamma_names)))]
  
  gamma_vals <- unlist(x[gamma_names], use.names = FALSE)
  gamma_df <- as.list(gamma_vals)
  names(gamma_df) <- paste0(gamma_names, "_hat")   # gamma0_hat, gamma1_hat, ...
  
  # --- pi_hat and pi_target dynamically ---
  pi_hat    <- x$pi_vec
  pi_target <- x$target_pi_vec
  K <- length(pi_hat)
  
  pi_hat_df <- as.list(pi_hat);    names(pi_hat_df) <- paste0("pi", 0:(K-1), "_hat")
  pi_tar_df <- as.list(pi_target); names(pi_tar_df) <- paste0("pi", 0:(K-1), "_target")
  
  tibble::as_tibble(c(
    list(
      scena     = x$scena,
      obj       = x$obj,
      converged = x$converged
    ),
    gamma_df,
    pi_hat_df,
    pi_tar_df
  ))
})
gamma_res

all(as.data.frame(gamma_res%>%dplyr::select(matches("^pi[0-9]_hat$")) ) - 
  as.data.frame(gamma_res%>%dplyr::select(matches("^pi[0-9]_target$")) )  < 1e-5)

scenarios_for_super_pop_updated <- scenarios_for_super_pop %>%
  dplyr::select(-contains("target"),-contains("gamma")) %>%
  left_join(
    gamma_res %>% 
      dplyr::select(scena, matches("^gamma[0-9]_hat$")),
    by = "scena"
  ) %>%
  rename_with(~ sub("_hat$", "", .x), matches("^gamma[0-9]_hat$")) 

scenarios_for_super_pop_updated %>% write_csv("scenarios_for_superpop_updated_beta1.5_4cat.csv")
