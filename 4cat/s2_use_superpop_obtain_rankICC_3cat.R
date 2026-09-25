library(rankICC)
library(dplyr)
library(readr)
library(geepack)
library(ordinal)
library(bridgedist)
source("~/Composite_Endpoints/proj2/functions_clean0903.R")
source("~/Composite_Endpoints/proj2/collapse_cats.R")

nsim = 1e2

epsilon_scale = 1; epsilon_sd=pi/sqrt(3)
n_cat <- 4

scenarios_for_super_pop_mini = read_csv("~/Composite_Endpoints/proj2_setup_011926/scenarios_for_superpop_updated_beta1.5_4cat.csv") %>%
  mutate(number_clusters= rep(500,
                              18),
         N = number_clusters*cluster_sizes,
         m = N/2) 
scenarios_for_super_pop_mini
scenarios_for_super_pop = scenarios_for_super_pop_mini

index = Sys.getenv("SLURM_ARRAY_TASK_ID")
scena = as.numeric(index)   
scena # 1-18

scenarios_for_super_pop[scena,] %>% as.data.frame
# extract all gamma columns for this scenario
gamma_vals <- scenarios_for_super_pop[scena,] %>%
  dplyr::select(dplyr::matches("^gamma[0-9]$")) %>%
  as.numeric
############### update everything to more than 3 categories.
# gamma0 <- scenarios_for_super_pop$gamma0[scena]
# gamma1 <- scenarios_for_super_pop$gamma1[scena]
cut_points <- gamma_vals

latent_ICC     <- scenarios_for_super_pop$latent_ICC[scena]
beta1          <- scenarios_for_super_pop$beta_X_list[scena] # used negative because WH assumes smaller category is better condition.
cluster_size   <- scenarios_for_super_pop$cluster_sizes[scena]
number_cluster <- scenarios_for_super_pop$number_clusters[scena]
m              <- scenarios_for_super_pop$m[scena]
N              <- scenarios_for_super_pop$N[scena]
ordinal_dist_type <- scenarios_for_super_pop$ordinal_dist_type[scena]

phi0 = get_phi_bridge_from_ICC(ICC=latent_ICC) #get_phi_bridge(var_bridge=sd_bridge^2)
phi0

### Derive Marginal category probabilities (per arm) via integration
num_cuts <- length(cut_points) # number of cutpoints
K <- n_cat                  # number of categories (0..num_cuts)

p_C_vec <- numeric(n_cat)
p_T_vec <- numeric(n_cat)

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

rankicc_result = data.frame(scena,mean_rank_icc=NA,mean_running_time=NA)

rankicc_sim = matrix(NA,nrow = nsim,ncol=2)
colnames(rankicc_sim)= c("rank_icc","running_time") 

# IDs and treatment (cluster-randomized: half control, half treatment)
cluster_id    <- rep(1:number_cluster, each = cluster_size)
treat_cluster <- rep(0:1, each = number_cluster/2)
X1            <- rep(treat_cluster, each = cluster_size)

collapse_ordinal_Y <- function(Y, collapse_idx) {
  
  # convert to 0-based values
  vals <- sort(collapse_idx - 1)
  
  a <- vals[1]
  b <- vals[2]
  
  Y_new <- Y
  
  # collapse upper one into lower one
  Y_new[Y == b] <- a
  
  # shift everything above b down by 1
  Y_new[Y > b] <- Y_new[Y > b] - 1
  
  Y_new
}

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
  
  Y_new <- collapse_ordinal_Y(Y, collapse_idx)
  table(Y, Y_new)
  
  ## assemble the data 
  # NOTE: ordgee expects Y to be ordered factor
  sim_data_cont <- data.frame(
    cluster = factor(cluster_id),
    Y_cont  = Y_new
   ) %>% group_by(cluster) %>%
    mutate(id = as.factor(dplyr::row_number())) %>%
    ungroup
  
  start_t = Sys.time()
  #Assigning equal weights to clusters
  rank_icc = rankICC(sim_data_cont$Y_cont, sim_data_cont$cluster, weights = "clusters")[1]
  end_t = Sys.time()
  running_time = difftime(end_t,start_t,units='secs')
  
  rankicc_sim[si,] = as.numeric(c(rank_icc, running_time) ) 
}

rankicc_result[1,-1] = colMeans(rankicc_sim)

result_path = "~/Composite_Endpoints/proj2_setup_011926/4cat/results/superpop_rankICC/" 

rankicc_result %>% 
  write_csv(paste0(result_path,"superpop_rankICC_",scena,".csv"))



