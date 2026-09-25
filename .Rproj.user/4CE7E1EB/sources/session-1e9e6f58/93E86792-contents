library(dplyr)
library(readr)
library(geepack)
library(ordinal)
library(bridgedist)
source("./functions/functions_clean0903.R")

scenarios_addrhomat = read_csv("emp_power_scenarios_ushape_ricc_022526.csv")


num_scena <- nrow(scenarios_addrhomat)

epsilon_scale = 1; epsilon_sd=pi/sqrt(3)
n_cat <- 3

power_results <- data.frame(
  WH_DE_latent_icc = NA,
  WH_DE_rank_icc = NA,
  WH_DE_rho11_icc = NA,
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
  gamma0 <- scenarios_addrhomat$gamma0[scena]
  gamma1 <- scenarios_addrhomat$gamma1[scena]
  cut_points <- c(gamma0, gamma1)
  
  latent_ICC     <- scenarios_addrhomat$latent_ICC[scena]
  beta1          <- scenarios_addrhomat$beta_X_list[scena] # used negative because WH assumes smaller category is better condition.
  cluster_size   <- scenarios_addrhomat$cluster_sizes[scena]
  number_cluster <- scenarios_addrhomat$number_clusters[scena]
  m              <- scenarios_addrhomat$m[scena]
  N              <- scenarios_addrhomat$N[scena]
  sd_bridge = get_b_sd(ICC=latent_ICC)
  phi0 = get_phi_bridge_from_ICC(ICC=latent_ICC) #get_phi_bridge(var_bridge=sd_bridge^2)
  phi0
  
  marginal_lor_by_bridge = phi0*beta1
  marginal_lor_by_bridge # 0.401
  
  ### Derive Marginal category probabilities (per arm) via integration
  # Control
  p0_C <- ordinal_logit_marginal_prob_bridge(0, x=0, beta1=beta1, cutpoints=cut_points, epsilon_scale=epsilon_scale, b_phi=phi0)
  p1_C <- ordinal_logit_marginal_prob_bridge(1, x=0, beta1=beta1, cutpoints=cut_points, epsilon_scale=epsilon_scale, b_phi=phi0)
  p2_C <- ordinal_logit_marginal_prob_bridge(2, x=0, beta1=beta1, cutpoints=cut_points, epsilon_scale=epsilon_scale, b_phi=phi0)
  # Treatment
  p0_T <- ordinal_logit_marginal_prob_bridge(0, x=1, beta1=beta1, cutpoints=cut_points, epsilon_scale=epsilon_scale, b_phi = phi0)
  p1_T <- ordinal_logit_marginal_prob_bridge(1, x=1, beta1=beta1, cutpoints=cut_points, epsilon_scale=epsilon_scale, b_phi=phi0)
  p2_T <- ordinal_logit_marginal_prob_bridge(2, x=1, beta1=beta1, cutpoints=cut_points, epsilon_scale=epsilon_scale, b_phi=phi0)
  
  p0_C;p1_C; p2_C
  p0_T;p1_T;p2_T
  pi_vec <- c(p0_C+p0_T, p1_C+p1_T, p2_C+p2_T)/2
  pi_vec
  
  P0 <- c(p0_C,p1_C, p2_C)
  P1 <- c(p0_T,p1_T, p2_T)
  
  # cumulative probs
  C0_C <- P0[1]
  C1_C <- P0[1] + P0[2]
  
  C0_T <- P1[1]
  C1_T <- P1[1] + P1[2]
  
  # Fisher info for each cumulative split
  I1 <- 1/(C0_C*(1-C0_C)) + 1/(C0_T*(1-C0_T))
  I2 <- 1/(C1_C*(1-C1_C)) + 1/(C1_T*(1-C1_T))
  
  c(I1 = I1, I2 = I2, ratio = I1/I2)
  
  # log( p0_C*(1-p0_T)/(p0_T*(1-p0_C)) ) # 0.401
  
  rho_mat = as.numeric(scenarios_addrhomat[scena,c("rho11","rho22","rho12")])
  
  rho_matrix = matrix(c(rho_mat[1],rho_mat[3],
                        rho_mat[3],rho_mat[2]),ncol=2) 
  rho_mat11 = rho_mat[1]
  rho_bin2 = as.numeric(scenarios_addrhomat[scena,c("bin_rho01")])
  
  #mean_rho_mat = mean(rho_matrix)
  #mean_diag = mean(diag(rho_matrix))
  #mean_bin_rho = mean(as.numeric(scenarios_addrhomat[scena,c("bin_rho0","bin_rho01")]))
  
  rankICC = as.numeric(scenarios_addrhomat[scena,c("mean_rank_icc")])
  #contICC = as.numeric(scenarios_addrhomat[scena,c("mean_lme_icc")])
  #### estimating powers
  # WH_DE_latent_icc
  power_WH_DE_latent_icc = get_power_eq10(
    log_OR       = marginal_lor_by_bridge,   # use ordgee lor ? or marginal log-OR 
    prop_categories = pi_vec,        # average category proportions across arms
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
  
  ### WH_DE_rho11_icc
  power_WH_DE_rho11_icc = get_power_eq10(
    log_OR       = marginal_lor_by_bridge,   # use ordgee lor ? or marginal log-OR
    prop_categories = pi_vec,        # average category proportions across arms
    ICC          = rho_mat11,         #
    cluster_size = cluster_size,
    m            = m,
    alpha        = 0.05
  )
  power_WH_DE_rho11_icc
  
  # WH
  power_WH = wh_power(theta_R=marginal_lor_by_bridge, n=N,pr=pi_vec,A=1)
  power_WH
  
  # 
  # ### WH_DE_lme_icc 
  # power_WH_DE_lme_icc = get_power_eq10(
  #   log_OR       = marginal_lor_by_bridge,   # use ordgee lor ? or marginal log-OR 
  #   prop_categories = pi_vec,        # average category proportions across arms
  #   ICC          = contICC,         #
  #   cluster_size = cluster_size,
  #   m            = m,
  #   alpha        = 0.05
  # )
  # power_WH_DE_lme_icc
  
  ### GEE binary
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
  
  # GEE exchangeable
  # p0_T_oppo = p2_T
  # p1_T_oppo = p1_T
  # p2_T_oppo = p0_T
  # p0_C_oppo = p2_C
  # p1_C_oppo = p1_C
  # p2_C_oppo = p0_C
  
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
  
  power_results[scena,] = c(power_WH_DE_latent_icc,
                            power_WH_DE_rank_icc,
                            power_WH_DE_rho11_icc,
                            power_WH,
                            # power_WH_DE_lme_icc,
                            #power_WH_DE_clmm_icc,
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
power_results1 %>% write_csv("est_power_u_rank_022526.csv")
