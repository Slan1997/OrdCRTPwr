### obtain rank ICC 
# (didn't include in the main analysis of 
# the manuscript since the performance of rank ICC was not 
# as good for ordinal outcomes with small number of categories.)

# this file was run on ACCRE cluster

## output: results/3cat/rankICCs_from_superpop.csv

library(rankICC)
library(dplyr)
library(readr)
library(geepack)
library(ordinal)
library(bridgedist)

invisible( lapply( list.files("../functions",
                              pattern = "\\.R$",
                              full.names = TRUE), source ) )

nsim = 1e2

epsilon_scale = 1; epsilon_sd=pi/sqrt(3)
n_cat <- 3

### use smaller number of clusters due to high computational cost
scenarios_for_super_pop_mini = read_csv("../../results/3cat/superpop/scenarios_for_superpop_3cat.csv") %>%
  mutate(number_clusters= rep(500,36),
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

rankicc_result = data.frame(scena,mean_rank_icc=NA,mean_running_time=NA)

rankicc_sim = matrix(NA,nrow = nsim,ncol=2)
colnames(rankicc_sim)= c("rank_icc","running_time") 

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
  # NOTE: ordgee expects Y to be ordered factor
  sim_data_cont = data.frame(
    cluster = factor(cluster_id),
    Y_cont  = Y
  ) %>% group_by(cluster) %>%
    mutate(id = as.factor(dplyr::row_number())) %>%
    ungroup
  start_t = Sys.time()
  #Assigning equal weights to clusters
  rank_icc = rankICC(sim_data_cont$Y_cont, sim_data_cont$cluster, weights = "clusters")[1]
  end_t = Sys.time()
  running_time = difftime(end_t,start_t,units='secs')
  
  rankicc_sim[si,] = as.numeric(c(rank_icc, running_time) ) 
  
  # #Assigning equal weights to observations 
  # rankICC::rankICC(dat$x, dat$cluster, weights = "obs")
  # #Iterative weighting based on the effective sample size
  # rankICC::rankICC(dat$x, dat$cluster, weights = "ess")
  # #Iterative weighting based on the combination of equal weights for clusters and equal weights for observations
  # rankICC::rankICC(dat$x, dat$cluster, weights = "combination")
}

rankicc_result[1,-1] = colMeans(rankicc_sim)

# result_path = "~/Composite_Endpoints/proj2_setup_011926/results/superpop_rankICC/"
# rankicc_result %>% 
#   write_csv(paste0(result_path,"superpop_rankICC_",scena,".csv"))


# ### combine results from each scenarios
# result_path1 = "~/Composite_Endpoints/proj2_setup_011926/results/superpop_rankICC/"
# ids = gsub("superpop_rankICC_([0-9]+).csv","\\1",list.files(result_path1)) %>% as.numeric
# (1:num_scena)[! 1:num_scena  %in% ids]
# 
# rankICC_true = NULL
# for (scena in 1:num_scena){
#   rankICC_true= bind_rows(rankICC_true,read_csv(paste0(result_path1,"superpop_rankICC_",scena,".csv")))
# }
# rankICC_true
# # > mean(as.numeric(rankICC_true$mean_running_time))
# # [1] 318
# # final1 = full_join(rankICC_true,rho_mat_true, by="scena")
# # final1 
# 
# final_out = left_join(scenarios_for_super_pop, rankICC_true #final1
#                       , by="scena")
# final_out %>% 
#   write_csv("~/Composite_Endpoints/proj2_setup_011926/results/rankICCs_from_superpop.csv")


