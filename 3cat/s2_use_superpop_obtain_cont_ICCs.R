### obtain continuous ICC and rank ICC
library(aod)
library(lme4)
#library(rankICC)
library(dplyr)
library(readr)
library(geepack)
library(ordinal)
library(bridgedist)
source("~/Composite_Endpoints/proj2/functions_clean0903.R")

nsim = 1e2

epsilon_scale = 1; epsilon_sd=pi/sqrt(3)
n_cat <- 3

scenarios_for_super_pop_mini = read_csv("scenarios_for_superpop_updated_beta1.5and2.csv") %>%
  mutate(number_clusters= rep(500,30),
         N = number_clusters*cluster_sizes,
         m = N/2) 
scenarios_for_super_pop_mini
scenarios_for_super_pop = scenarios_for_super_pop_mini


index = Sys.getenv("SLURM_ARRAY_TASK_ID")
scena = as.numeric(index)   
scena # 1-15

scenarios_for_super_pop[scena,] %>% as.data.frame
gamma0 <- scenarios_for_super_pop$gamma0[scena]
gamma1 <- scenarios_for_super_pop$gamma1[scena]
cut_points <- c(gamma0, gamma1)

latent_ICC     <- scenarios_for_super_pop$latent_ICC[scena]
beta1          <- scenarios_for_super_pop$beta_X_list[scena] # used negative because WH assumes smaller category is better condition.
cluster_size   <- scenarios_for_super_pop$cluster_sizes[scena]
number_cluster <- scenarios_for_super_pop$number_clusters[scena]
m              <- scenarios_for_super_pop$m[scena]
N              <- scenarios_for_super_pop$N[scena]

phi0 = get_phi_bridge_from_ICC(ICC=latent_ICC) #get_phi_bridge(var_bridge=sd_bridge^2)
phi0


icc_result = data.frame(scena,mean_lme_icc=NA,mean_anova_icc=NA,
                        mean_lme_running_time=NA,mean_anova_running_time=NA)

icc_sim = matrix(NA,nrow = nsim,ncol=4)
colnames(icc_sim)= c("lme_icc","anova_icc","lme_running_time","anova_running_time") 

# IDs and treatment (cluster-randomized: half control, half treatment)
cluster_id    <- rep(1:number_cluster, each = cluster_size)
treat_cluster <- rep(0:1, each = number_cluster/2)
X1            <- rep(treat_cluster, each = cluster_size)

message("Scenario ", scena)
set.seed(scena + 11532)


for (si in 1:nsim){
  # random intercept per cluster
  b_cluster <- rep(rbridge(number_cluster, phi=phi0), each = cluster_size)
  # logistic errors on latent scale
  epsilon <- rlogis(N, location = 0, scale = epsilon_scale)
  
  # latent + thresholds -> ordinal Y in {0,1,2}
  eta <- beta1 * X1 + b_cluster + epsilon # with this setup, treatment => larger eta, higher probability of Y=2
  # this corresponds to PO model: logit(P(Y<=k)) = gamma_k - beta1*treat
  Y   <- cut(eta, breaks = c(-Inf, cut_points, Inf), labels = FALSE) - 1
  
  ## assemble the data 
  sim_data_cont = data.frame(
    cluster = factor(cluster_id),
    Y_cont  = Y
  ) %>% group_by(cluster) %>%
    mutate(id = as.factor(dplyr::row_number())) %>%
    ungroup
  
  ### lme
  start_t1 = Sys.time()
  fit <- lmer(Y_cont ~ 1 + (1 | cluster), data = sim_data_cont)
  vc <- as.data.frame(VarCorr(fit))
  sigma_b2 <- vc[vc$grp == "cluster", "vcov"]
  sigma_w2 <- vc[vc$grp == "Residual", "vcov"]
  
  lme_ICC <- sigma_b2 / (sigma_b2 + sigma_w2)
  end_t1 = Sys.time()
  lme_running_time = difftime(end_t1,start_t1,units='secs')
  
  ### anova
  start_t2 = Sys.time()
  aov_fit <- aov(Y_cont ~ cluster, data = sim_data_cont)
  ms <- summary(aov_fit)[[1]][, "Mean Sq"]
  MSB <- ms[1]
  MSW <- ms[2]
  mbar <- mean(table(sim_data_cont$cluster))
  
  anova_ICC <- (MSB - MSW) / (MSB + (mbar - 1) * MSW)
  
  end_t2 = Sys.time()
  anova_running_time = difftime(end_t2,start_t2,units='secs')
  
  icc_sim[si,] = as.numeric(c( lme_ICC,anova_ICC,lme_running_time ,anova_running_time) ) 
}

icc_result[1,-1] = colMeans(icc_sim)

result_path = "~/Composite_Endpoints/proj2_setup_011926/results/superpop_cont_ICC/"

icc_result %>% 
  write_csv(paste0(result_path,"superpop_contICC_",scena,".csv"))



