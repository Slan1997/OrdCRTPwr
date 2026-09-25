### collapse 4 cat into 3 cat
# 1. even: (.25,.25,.25,.25)
# 2. rare: (.6,.2,.1,.1)
# 3. common: (.1,.1,.2,.6)
# 4. bell: (.1,.4,.4,.1)
# 5. U-shape: (.4,.1,.1,.4)

### confirmed. the pi in Corr(Yis,Yit) has to be marginal!
## correct sign of beta1 => ordinal dist type should follow the expectations.
### the result of this file will later be combined into "~/Composite_Endpoints/proj2/results/rho_mat_from_superpop0929.csv"
library(dplyr)
library(readr)
library(geepack)
library(ordinal)
library(bridgedist)
source("~/Composite_Endpoints/proj2/functions_clean0903.R")

nsim = 1e2

epsilon_scale = 1; epsilon_sd=pi/sqrt(3)
n_cat <- 4

scenarios_for_super_pop = read_csv("~/Composite_Endpoints/proj2_setup_011926/scenarios_for_superpop_updated_beta1.5_4cat.csv")
# scenarios_for_super_pop
# tail(scenarios_for_super_pop)
index = Sys.getenv("SLURM_ARRAY_TASK_ID")
scena = as.numeric(index)   
scena # 1-15

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

add_info <- data.frame(scena = scena)
# Generate column names
pi_names    <- paste0("pi_", 0:(K - 1))
pi_C_names  <- paste0("pi_", 0:(K - 1), "_C")
pi_T_names  <- paste0("pi_", 0:(K - 1), "_T")
# Add marginal (averaged) probabilities
add_info[1, pi_names] <- pi_vec
# Add arm-specific probabilities
add_info[1, pi_C_names] <- p_C_vec
add_info[1, pi_T_names] <- p_T_vec


####### this step we need to transfer this 4-categorical outcome into 3-categorical
# and obtain the rho11, rho 22, rho12.
### do we also wanna try GEE estimator with 4 categories??

### the 3-category outcome version
rho_mat_result = data.frame(scena,rho11=NA,rho22=NA,rho12=NA)

### the full version
rho_names <- as.vector(
  outer(
    1:(n_cat - 1), 1:(n_cat - 1),
    function(i, j) paste0("rho", i, j)
  )
)
rho_mat_result_ncat <- data.frame(
  scena = scena,
  matrix(NA, nrow = 1, ncol = length(rho_names))
)

rho_mat_sim = matrix(NA,nrow = nsim,ncol=3)
colnames(rho_mat_sim) = c("rho11","rho22","rho12")

rho_mat_sim_ncat = matrix(NA,nrow = nsim,ncol=length(rho_names))
colnames(rho_mat_sim) = rho_names

sim_pi_vec = matrix(NA,nrow = nsim,ncol=9)
colnames(sim_pi_vec)=c("pi_0","pi_1","pi_2",
                       "pi_0_C","pi_1_C","pi_2_C",
                       "pi_0_T","pi_1_T","pi_2_T")


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
  sim_data <- data.frame(
    cluster = factor(cluster_id),
    treat   = factor(X1),
    Y       = factor(Y, ordered = TRUE), # this Y follows Whitehead setup
    Y_oppo = factor(case_when(Y==0~2,Y==1~1,Y==2~0), ordered = TRUE), # this follows clmm, gee setup, will give positive estimate of treatment effect size
    b_cluster) %>% mutate(Y0 = ifelse(Y==0,1,0), # Y2
                          Y01 = ifelse(Y<=1,1,0)) # Y2 & Y1
  #sim_data
  sim_pi_vec[si,c("pi_0","pi_1","pi_2")] = prop.table(table(sim_data$Y))
  sim_pi_vec[si,c("pi_0_C","pi_1_C","pi_2_C")] = prop.table(table(sim_data$Y[sim_data$treat==0]))
  sim_pi_vec[si,c("pi_0_T","pi_1_T","pi_2_T")] = prop.table(table(sim_data$Y[sim_data$treat==1]))
  
  p_i_C_and_T_orig = tibble(treat=as.factor(c(0,1)),pi_0=c(p0_C,p0_T),pi_1=c(p1_C,p1_T))
  sim_data_bin_orig = sim_data %>% transmute(cluster,treat,Y,
                                             bin_Y0 = ifelse(Y==0,1,0),
                                             bin_Y1 = ifelse(Y==1,1,0)) %>%
    full_join(p_i_C_and_T_orig , by=c("treat"))
  e_mat = sim_data_bin_orig %>% mutate(e1 =  (bin_Y0-pi_0)/sqrt(pi_0*(1-pi_0)),
                                       e2 =  (bin_Y1-pi_1)/sqrt(pi_1*(1-pi_1)))
  rho_mat_full = e_mat %>% group_by(cluster) %>% 
    summarise(
      Ti = n(),
      rho11_i = (sum(outer(e1, e1, "*")) - sum(e1^2))/ (Ti * (Ti - 1) ),
      rho22_i = (sum(outer(e2, e2, "*")) - sum(e2^2))/ (Ti * (Ti - 1) ),
      rho12_i = (sum(outer(e1, e2, "*")) - sum(e1*e2))/ (Ti * (Ti - 1) ) 
    )
  rho_mat_orig = rho_mat_full %>% dplyr::select(rho11_i ,rho22_i, rho12_i ) %>% 
    colMeans 
  
  #### mean alpha for two binary outcomes
  binY_p_i_C_and_T_orig = tibble(treat=as.factor(c(0,1)),pi_0=c(p0_C,p0_T),pi_01=c((p0_C+p1_C),(p0_T+p1_T)))
  binY_sim_data_bin_orig = sim_data %>% transmute(cluster,treat,Y,Y0,Y01) %>%
    full_join(binY_p_i_C_and_T_orig , by=c("treat"))
  e_both = binY_sim_data_bin_orig %>% mutate(e0 =  (Y0-pi_0)/sqrt(pi_0*(1-pi_0)),
                                             e01 =  (Y01-pi_01)/sqrt(pi_01*(1-pi_01)))
  rho_both_full = e_both %>% group_by(cluster) %>% 
    summarise(
      Ti = n(),
      rho0_i = (sum(outer(e0, e0, "*")) - sum(e0^2))/ (Ti * (Ti - 1) ),
      rho01_i = (sum(outer(e01, e01, "*")) - sum(e01^2))/ (Ti * (Ti - 1) )
    )
  rho_both = rho_both_full %>% dplyr::select(rho0_i,rho01_i) %>% colMeans 
  
  rho_mat_sim[si,] = as.numeric(c(rho_mat_orig, rho_both) )
}

result_path = "~/Composite_Endpoints/proj2_setup_011926/results/superpop_rhomat/"

rho_mat_result[1,-1] = colMeans(rho_mat_sim)


add_info[1,paste0(rep("mean.",9),colnames(sim_pi_vec))] = colMeans(sim_pi_vec)


rho_mat_result %>%left_join(add_info,by="scena") %>% 
  write_csv(paste0(result_path,"superpop_rho_mat_",scena,".csv"))
