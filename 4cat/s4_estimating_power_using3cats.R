library(dplyr)
library(readr)
library(geepack)
library(ordinal)
library(bridgedist)
source("../functions/functions_clean0903.R")
source("../functions/get_pi_from_coef.R")
source("../functions/collapse_cats.R")
scenarios_addrhomat = read_csv("emp_power_scenarios_4cat.csv")
num_scena <- nrow(scenarios_addrhomat)

epsilon_scale = 1; epsilon_sd=pi/sqrt(3)
n_cat <- 4

power_results <- data.frame(
  WH_DE_latent_icc = NA,
  WH_DE_rank_icc = NA,
  WH = NA ,
  GEE_bin1 = NA,
  GEE_bin2 = NA,
  GEE_exch = NA,
  GEE_indep = NA,
  marginal_lor_by_bridge=NA,
  theta_R_bin = NA,
  theta_R_bin2 = NA
)

for (scena in 1:num_scena){
  scenarios_addrhomat[scena,] %>% as.data.frame
  gamma_vals <- scenarios_addrhomat[scena,] %>%
    dplyr::select(dplyr::matches("^gamma[0-9]$")) %>% as.numeric
  cut_points <- gamma_vals
  
  latent_ICC     <- scenarios_addrhomat$latent_ICC[scena]
  beta1          <- scenarios_addrhomat$beta_X_list[scena] # used negative because WH assumes smaller category is better condition.
  cluster_size   <- scenarios_addrhomat$cluster_sizes[scena]
  number_cluster <- scenarios_addrhomat$number_clusters[scena]
  m              <- scenarios_addrhomat$m[scena]
  N              <- scenarios_addrhomat$N[scena]
  ordinal_dist_type <- scenarios_addrhomat$ordinal_dist_type[scena]

  sd_bridge = get_b_sd(ICC=latent_ICC)
  phi0 = get_phi_bridge_from_ICC(ICC=latent_ICC) #get_phi_bridge(var_bridge=sd_bridge^2)
  phi0
  
  marginal_lor_by_bridge = phi0*beta1
  marginal_lor_by_bridge # 0.401
  
  p_C_vec <- numeric(n_cat)
  p_T_vec <- numeric(n_cat)
  
  num_cuts <- length(cut_points) 
  for (k in 0:num_cuts) {
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
  p_C_vec;p_T_vec
  pi_vec <- (p_C_vec + p_T_vec) / 2
  pi_vec
  
  ### collapse 4 cat into 3 cat
  # 1. even: (.25,.25,.25,.25)
  # 2. rare: (.6,.2,.1,.1)
  # 3. common: (.1,.1,.2,.6)
  # 4. bell: (.1,.4,.4,.1)
  # 5. bell2: (.1,.6,.2,.1)
  # 6. U-shape: (.4,.1,.1,.4)
  if (ordinal_dist_type %in% c("even", "bell","bell2")) {
    collapse_idx <- c(3, 4)
  } else {
    collapse_idx <- order(pi_vec)[1:2]
    # check if consecutive
    if (diff(sort(collapse_idx)) != 1) {
      warning(
        sprintf(
          "collapse_idx not consecutive: (%s)",
          paste(collapse_idx, collapse = ", ")
        )
      )
    }
  }
  
  pi_vec_collapse = collapse_cats(pi_vec, collapse_idx)
  p_C_vec_collapse = collapse_cats(p_C_vec, collapse_idx)   # merge 3&4
  p_T_vec_collapse = collapse_cats(p_T_vec, collapse_idx) 
  
  p0_C = p_C_vec_collapse[1] ; p1_C = p_C_vec_collapse[2]
  p0_T = p_T_vec_collapse[1] ; p1_T = p_T_vec_collapse[2]
  
  n_cat = 3
  K <- n_cat 
  
  rho_mat = as.numeric(scenarios_addrhomat[scena,c("rho11","rho12","rho22")])
  
  rho_matrix = matrix(c(rho_mat[1],rho_mat[2],
                        rho_mat[2],rho_mat[3]),ncol=2) 
  
  rho_mat11 = rho_mat[1]
  rho_bin2 = as.numeric(scenarios_addrhomat[scena,c("bin_rho01")])
  
  
  rankICC = as.numeric(scenarios_addrhomat[scena,c("mean_rank_icc")])
  
  #### estimating powers
  # WH_DE_latent_icc
  power_WH_DE_latent_icc = get_power_eq10(
    log_OR       = marginal_lor_by_bridge,   # use ordgee lor ? or marginal log-OR 
    prop_categories = pi_vec_collapse,        # average category proportions across arms
    ICC          = latent_ICC,         # latent ICC 
    cluster_size = cluster_size,
    m            = m,
    alpha        = 0.05
  )
  power_WH_DE_latent_icc
  
  ### WH_DE_rank_icc
  power_WH_DE_rank_icc = get_power_eq10(
    log_OR       = marginal_lor_by_bridge,   # use ordgee lor ? or marginal log-OR
    prop_categories = pi_vec,        # average category proportions across arms
    ICC          = rankICC,         #
    cluster_size = cluster_size,
    m            = m,
    alpha        = 0.05
  )
  power_WH_DE_rank_icc
  
  # WH
  power_WH = wh_power(theta_R=marginal_lor_by_bridge, n=N,pr=pi_vec_collapse,A=1)
  power_WH
  
  ex_V = exch_V_theta(#p0_T_oppo,p1_T_oppo,p0_C_oppo,p1_C_oppo,
    p0_T,p1_T,p0_C,p1_C,
    rho_mat=rho_matrix, 
    cluster_size=cluster_size,number_cluster = number_cluster)
  power_GEE_exch = power_gee_est( V_theta=ex_V ,theta_star = marginal_lor_by_bridge)
  power_GEE_exch
  
  # GEE indep
  indep_V = exch_V_theta(#p0_T_oppo,p1_T_oppo,p0_C_oppo,p1_C_oppo,
    p0_T,p1_T,p0_C,p1_C,
    rho_mat=matrix(0,2,2),
    cluster_size,number_cluster = number_cluster)
  power_GEE_indep = power_gee_est(V_theta=indep_V ,theta_star = marginal_lor_by_bridge)
  power_GEE_indep
  
  marginal_lor_by_bridge
  theta_R_bin <- qlogis(p0_C) - qlogis(p0_T)
  theta_R_bin
  
  power_GEE_bin = power_binary(P0 = p0_C, P1 = p0_T, theta_R=marginal_lor_by_bridge, 
                               rho = rho_mat11, alpha = 0.05,
                               cluster_size = cluster_size,
                               number_cluster=number_cluster,allocate_prop_C=.5)
  
  
  theta_R_bin2 <- qlogis(p0_C+p1_C) - qlogis(p0_T+p1_T) # same as marginal_lor_by_bridge
  theta_R_bin2
  power_GEE_bin2 = power_binary(P0 = (p0_C+p1_C), P1 = (p0_T+p1_T), 
                                theta_R=marginal_lor_by_bridge, 
                                rho = rho_bin2, alpha = 0.05,
                                cluster_size = cluster_size,
                                number_cluster=number_cluster,allocate_prop_C=.5)
  
 
  power_results[scena,] = c(power_WH_DE_latent_icc,
                            power_WH_DE_rank_icc,
                            power_WH,
                            power_GEE_bin,
                            power_GEE_bin2,
                            power_GEE_exch,
                            power_GEE_indep,
                            marginal_lor_by_bridge,
                            theta_R_bin,
                            theta_R_bin2
  )
}

power_results
power_results1 = bind_cols(scenarios_addrhomat,power_results)
power_results1 %>% write_csv("est_power_4cat.csv")


