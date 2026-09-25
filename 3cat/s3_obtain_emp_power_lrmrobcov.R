library(dplyr)
library(rms)
library(readr)
library(ordinal)
library(bridgedist)
source("~/Composite_Endpoints/proj2/functions_clean0903.R")

nsim = 1e4

result_path = "~/Composite_Endpoints/proj2_setup_011926/results/lrm_emp_power/"

scenarios_addrhomat = read_csv("~/Composite_Endpoints/proj2_setup_011926/results/emp_power_scenarios.csv")

index = Sys.getenv("SLURM_ARRAY_TASK_ID")
scena = as.numeric(index)   
scena

power_results <- data.frame(
  scena, 
  Empirical_lrm = NA,
  Empirical_lrm_small_K_1 = NA,
  Empirical_lrm_small_K_2 = NA
)

epsilon_scale = 1; epsilon_sd=pi/sqrt(3)
n_cat <- 3

gamma0 <- scenarios_addrhomat$gamma0[scena]
gamma1 <- scenarios_addrhomat$gamma1[scena]
cut_points <- c(gamma0, gamma1)

latent_ICC     <- scenarios_addrhomat$latent_ICC[scena]
beta1          <- scenarios_addrhomat$beta_X_list[scena]  # larger category => better condition
cluster_size   <- scenarios_addrhomat$cluster_sizes[scena]
number_cluster <- scenarios_addrhomat$number_clusters[scena]
m              <- scenarios_addrhomat$m[scena]
N              <- scenarios_addrhomat$N[scena]
small_cluster  <- scenarios_addrhomat$small_cluster[scena]

sd_bridge = get_b_sd(ICC=latent_ICC)
phi0 = get_phi_bridge_from_ICC(ICC=latent_ICC) #get_phi_bridge(var_bridge=sd_bridge^2)
phi0

stopifnot(all.equal(sd_bridge^2/(sd_bridge^2 + epsilon_sd^2), latent_ICC, tol = 1e-6))

### this simulation is just to get empirical power.
sim_outputs = data.frame(scena = scena, sim = 1:nsim, 
                         ## p values
                         p_val_lrm = NA,
                         p_val_lrm_small_K_1 = NA,
                         p_val_lrm_small_K_2 = NA)

set.seed(scena + 11532)

# IDs and treatment (cluster-randomized: half control, half treatment)
cluster_id    <- rep(1:number_cluster, each = cluster_size)
treat_cluster <- rep(0:1, each = number_cluster/2)
X1            <- rep(treat_cluster, each = cluster_size)


NA_ct_lrm = 0
#ex_notok_ct = 0
#NA_ct_indep = 0 

for (si in 1:nsim){
  if (si %% 100 == 0) message("  sim ", si)
  
  # random intercept per cluster
  b_cluster <- rep(rbridge(number_cluster, phi=phi0), each = cluster_size)
  # logistic errors on latent scale
  epsilon <- rlogis(N, location = 0, scale = epsilon_scale)
  
  # latent + thresholds -> ordinal Y in {0,1,2}
  eta <- beta1 * X1 + b_cluster + epsilon
  Y   <- cut(eta, breaks = c(-Inf, cut_points, Inf), labels = FALSE) - 1
  
  ## assemble the data 
  # NOTE: ordgee expects Y to be ordered factor
  sim_data <- data.frame(
    cluster = factor(cluster_id),
    treat   = factor(X1),
    Y       = factor(Y, ordered = TRUE)#, # this Y follows Whitehead setup
  )
  sim_data
  
  dd <- datadist(sim_data)
  options(datadist='dd')
  gee.mod <- lrm(as.formula(Y~treat),data=sim_data,y=TRUE,x=TRUE)
  gee <- robcov(gee.mod, cluster = sim_data$cluster)
  # gee
  test_mu = gee$coefficients[3]
  test_var = gee$var[3,3]
  p_val_lrm = 2*pnorm(abs(test_mu/sqrt(test_var)), lower.tail = F)
  
  t_df1 = number_cluster - 1 
  small_corr1 = number_cluster/t_df1
  t_df2 = number_cluster - 2
  small_corr2 = number_cluster/t_df2
  
  p_val_lrm_small_K_1 = 2*pt(abs(test_mu/sqrt(test_var*small_corr1)), df = t_df1, lower.tail = F)
  p_val_lrm_small_K_2 = 2*pt(abs(test_mu/sqrt(test_var*small_corr2)), df = t_df2, lower.tail = F)
  
  if(is.na(p_val_lrm))  NA_ct_lrm= NA_ct_lrm + 1
  
  sim_outputs[si,c("p_val_lrm")] = p_val_lrm
  sim_outputs[si,c("p_val_lrm_small_K_1")] = p_val_lrm_small_K_1
  sim_outputs[si,c("p_val_lrm_small_K_2")] = p_val_lrm_small_K_2
}

sim_outputs

# 1. Empirical power (use gee because clmm may have NaN)
power_emp_lrm = mean(sim_outputs$p_val_lrm < 0.05,na.rm=T)
power_emp_lrm

power_emp_lrm_K_1 = mean(sim_outputs$p_val_lrm_small_K_1 < 0.05,na.rm=T)
power_emp_lrm_K_1

power_emp_lrm_K_2 = mean(sim_outputs$p_val_lrm_small_K_2 < 0.05,na.rm=T)
power_emp_lrm_K_2

power_results[1,"Empirical_lrm"] = power_emp_lrm
power_results[1,"Empirical_lrm_small_K_1"] = power_emp_lrm_K_1
power_results[1,"Empirical_lrm_small_K_2"] = power_emp_lrm_K_2

power_results$NA_ct_lrm = NA_ct_lrm

power_results

power_results %>% write_csv(paste0(result_path,"emp_power_",scena,".csv"))
