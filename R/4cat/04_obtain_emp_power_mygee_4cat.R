# This file was run on ACCRE cluster
# Each SLURM_ARRAY_TASK_ID corresponds to one scenario.
## scenario-specific outputs were combined in the end into the final output: 
# ../../results/4cat/emp_power_all_4cat_052426.csv

library(dplyr)
library(rms)
library(readr)
library(ordinal)
library(bridgedist)

invisible( lapply( list.files("../functions",
                              pattern = "\\.R$",
                              full.names = TRUE), source ) )

scenarios_addrhomat = read_csv("../../results/4cat/emp_power_scenarios_4cat.csv")

print_info = F

index = Sys.getenv("SLURM_ARRAY_TASK_ID")
scena = as.numeric(index)   
scena  #1-72
scenarios_addrhomat[scena,] %>% as.data.frame

epsilon_scale = 1; epsilon_sd=pi/sqrt(3)
n_cat <- 4

target_success <- 1e4 #1e4
max_tries <- 2e4 #2e4    # safety cap
type_1_err <- 0.05
### for gee
tol <- 1e-4
tol_rho <- 1e-4
tol_abs   <- 1e-4   
tol_abs_rho <- 1e-4
tol_abs_score <- 1e-3
tot_bounce <- 1e-5
max_iter <- 1e3 #5e2
max_inc <- .5 

gamma_vals <- scenarios_addrhomat[scena,] %>%
  dplyr::select(dplyr::matches("^gamma[0-9]$")) %>%
  as.numeric
cut_points <- gamma_vals

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
# phi0*beta1 
stopifnot(all.equal(sd_bridge^2/(sd_bridge^2 + epsilon_sd^2), latent_ICC, tol = 1e-6))

# IDs and treatment (cluster-randomized: half control, half treatment)
cluster_id    <- rep(1:number_cluster, each = cluster_size)
treat_cluster <- rep(0:1, each = number_cluster/2)
X1            <- rep(treat_cluster, each = cluster_size)


#### add beta and var_beta output
est_beta_treat = numeric(0)
var_beta_treat = numeric(0)
sandwich_var_beta_treat = numeric(0)
est_rho_vec_mat = NULL
#beta_treat1_vec = NULL
pvals <- numeric(0)
n_tried <- 0L
n_bad   <- 0L
converge_ind <- NULL
count_iter <- NULL
set.seed(scena + 11532)

while (length(pvals) < target_success && n_tried < max_tries) {
  n_tried <- n_tried + 1L
  if (print_info) cat("n_tried:", n_tried, "\n")
  # ---- simulate one dataset (your existing code) ----
  b_cluster <- rep(rbridge(number_cluster, phi = phi0), each = cluster_size)
  epsilon   <- rlogis(N, location = 0, scale = epsilon_scale)
  eta <- beta1 * X1 + b_cluster + epsilon # with this setup, treatment => larger eta, higher probability of Y=2
  # this corresponds to PO model: logit(P(Y<=k)) = gamma_k - beta1*treat
  Y   <- cut(eta, breaks = c(-Inf, cut_points, Inf), labels = FALSE) - 1L
  if (print_info) cat(table(Y))
  
  sim_data <- data.frame(
    cluster = factor(cluster_id),
    treat   = factor(X1),
    Y       = factor(Y, ordered = TRUE)
  )
  
  dd <- datadist(sim_data)
  options(datadist='dd')
  gee.mod <- lrm(as.formula(Y~treat),data=sim_data,y=TRUE,x=TRUE)
  
  # ---- fit my gee ----
  beta_init = - gee.mod$coefficients  #- gee_fit0$beta
  names(beta_init) = c( paste0("y<=",(1:(n_cat-1))-1) ,  "treat=1")
  #rho_vec_init = rho_mat_update(beta=beta_init,sim_data=sim_data)
  rho_mat_init = rho_mat_update(beta=beta_init,sim_data=sim_data)
  rho_vec_init = rho_mat_init[upper.tri(rho_mat_init,diag = T)]
  
  #beta_init[3] = - max(abs(c(0.01,beta_init[3])))
  converged     <- FALSE
  bad_rho_range <- FALSE
  beta <- beta_init
  #rho_vec <- rho_vec_init
  #rho_mat <- matrix(c(rho_vec[1],rho_vec[3],rho_vec[3],rho_vec[2]),ncol=2) 
  rho_vec <- rho_vec_init
  rho_mat <- rho_mat_init
  
  beta_trace <- beta   # store the initial point
  rho_vec_trace <- rho_vec   # store the initial point
  
  score_all <- NULL
  y_by_cluster <- lapply(
    split(sim_data$Y, sim_data$cluster),
    function(yf) {
      y <- as.integer(as.character(yf))
      K <- n_cat - 1  # number of indicators
      # Build K x cluster_size matrix of indicators
      Ymat <- sapply(0:(K - 1), function(k) as.integer(y == k))
      # Stack by subject 
      c(t(Ymat))
    }
  )
  
  
  irls_ok <- TRUE
  
  tryCatch({
    
    ## ---- IRLS loop ----
    for (k in 1:max_iter) {
      step_beta <- 1
      if (print_info) cat("  k:", k, "\n")
      
      ## 1) build pi_hat from current beta
      pi_hat_info <- get_pi_from_coef(coef = beta, n_cat = n_cat)
      pi_C_hat <- pi_hat_info$pi_C_hat
      pi_T_hat <- pi_hat_info$pi_T_hat
      
      ## IMPORTANT: make treat_cluster numeric 0/1 (your treat is factor 0/1? ensure)
      treat_cluster_num <- sim_data %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(treat = as.integer(as.character(first(treat)))) %>%
        dplyr::arrange(cluster) %>%
        dplyr::pull(treat)
      
      ## 2) Fisher/IRLS increment (can throw solve() singular)
      beta_update <- beta_increment_ord(
        pi_T_hat = pi_T_hat,
        pi_C_hat = pi_C_hat,
        rho_matrix = rho_mat,
        cluster_size = cluster_size,
        number_cluster = number_cluster,
        y_by_cluster = y_by_cluster,
        treat_cluster = treat_cluster_num
      )
      
      increment <- as.numeric(beta_update$increment)
      score_now <- as.numeric(beta_update$part2)
      
      if (any(!is.finite(increment)) || any(!is.finite(score_now))) {
        warning(sprintf("Dataset %d: non-finite increment/score; skipping.", n_tried))
        bad_rho_range <- TRUE
        break
      }
      
      score_all <- rbind(score_all, score_now)
      
      if (max(abs(increment)) >= max_inc) {
        step_beta <- 1 / (max(abs(increment)) %/% max_inc) / 10
      }
      
      beta_next <- beta + step_beta * increment
      
      ## bouncing handling unchanged...
      if (k > 2) {
        beta_prev2 <- beta_trace[nrow(beta_trace) - 2, ]
        beta_prev1 <- beta_trace[nrow(beta_trace) - 1, ]
        if (max(abs(beta_next - beta_prev1)) < tot_bounce &&
            max(abs(beta - beta_prev2)) < tot_bounce) {
          cat("  Bouncing between two points at dataset", n_tried, "\n")
          step_beta <- 0.1
          beta_next <- beta + step_beta * increment
        }
      }
      
      ## 4) update rho (can also throw)
      rho_mat_next <- rho_mat_update(beta = beta_next, sim_data = sim_data)
      rho_vec_next <- rho_mat_next[upper.tri(rho_mat_next, diag = TRUE)]
      
      # range checks unchanged...
      if (any(is.nan(rho_vec_next))) {
        warning(sprintf("Dataset %d: rho has NaN; skipping.", n_tried))
        bad_rho_range <- TRUE
        break
      }
      if (min(rho_vec_next) < -1 || max(rho_vec_next) > 1) {
        warning(sprintf("Dataset %d: rho out of [-1,1]; skipping.", n_tried))
        bad_rho_range <- TRUE
        break
      }
      
      abs_change_beta <- max(abs(beta_next - beta))
      abs_change_rho_vec <- max(abs(rho_vec_next - rho_vec))
      rel_change_beta <- max(abs((beta_next - beta) / (beta + .Machine$double.eps)))
      rel_change_rho_vec <- max(abs((rho_vec_next - rho_vec) / (rho_vec + .Machine$double.eps)))
      
      if ((rel_change_beta < tol || abs_change_beta < tol_abs) &&
          (rel_change_rho_vec < tol_rho || abs_change_rho_vec < tol_abs_rho)) {
        converged <- TRUE
        beta_trace <- rbind(beta_trace, beta_next)
        rho_vec_trace <- rbind(rho_vec_trace, rho_vec_next)
        break
      } else if (max(abs(score_now)) < tol_abs_score) {
        converged <- TRUE
        beta_trace <- rbind(beta_trace, beta)
        rho_vec_trace <- rbind(rho_vec_trace, rho_vec)
        break
      }
      
      beta <- beta_next
      rho_vec <- rho_vec_next
      rho_mat <- rho_mat_next
      beta_trace <- rbind(beta_trace, beta)
      rho_vec_trace <- rbind(rho_vec_trace, rho_vec)
    }
    
  }, error = function(e) {
    irls_ok <<- FALSE
    warning(sprintf("Dataset %d: IRLS error: %s. Skipping.", n_tried, conditionMessage(e)))
  })
  
  if (print_info) {
    cat("  tail(beta_trace):\n");    print(tail(beta_trace))
    cat("  tail(rho_vec_trace):\n"); print(tail(rho_vec_trace))
  }
  
  if (!irls_ok) {
    n_bad <- n_bad + 1L
    next
  }
  
  if (bad_rho_range) {
    ## this dataset has invalid rho; treat as bad and skip everything else
    n_bad <- n_bad + 1L
    next  # go to next dataset in the while() loop
  }
  
  ## not bad_rho_range: check convergence
  count_iter    <- c(count_iter, k)
  converge_ind  <- c(converge_ind, converged)
  
  if (!converged && (k == max_iter)){
    warning("Hit max_iter without convergence")
    print(score_all)
    print(beta_trace)
    print(rho_vec_trace)
    n_bad <- n_bad + 1L
    next  # skip sandwich/p-values for this dataset
  }
  
  ## ---- only here: dataset converged and rho was in range ----
  beta_new    <- tail(beta_trace, 1)
  rho_vec_new <- tail(rho_vec_trace, 1)
  
  n_thresh <- n_cat - 1
  rho_mat_new <- matrix(0,  n_thresh,  n_thresh)
  rho_mat_new[upper.tri(rho_mat_new, diag = TRUE)] <- rho_vec_new
  rho_mat_new[lower.tri(rho_mat_new)] <- t(rho_mat_new)[lower.tri(rho_mat_new)]
  
  #### get pi_hat
  pi_hat_info_new = get_pi_from_coef(coef=beta_new, n_cat=n_cat)
  pi_C_hat = pi_hat_info_new$pi_C_hat  #pi_use[,1]
  pi_T_hat = pi_hat_info_new$pi_T_hat  #pi_use[,2]
  pi_C_hat ; pi_T_hat
  
  # gammas_est = sort(beta_new[1:2])
  # gamma1_est = gammas_est[1]
  # gamma2_est = gammas_est[2]
  # beta_treat_est = beta_new[3]
  # treat01 = c(0,1)
  # cum1_est = plogis(gamma1_est + beta_treat_est * treat01)
  # cum2_est = plogis(gamma2_est + beta_treat_est * treat01)
  # pi_1_est = cum1_est
  # pi_2_est = cum2_est - cum1_est
  # pi_hat_CT = cbind(pi_1_est,pi_2_est)  
  # pi_C_hat = matrix(as.numeric(pi_hat_CT[1,]),ncol=1)
  # pi_T_hat = matrix(as.numeric(pi_hat_CT[2,]),ncol=1)
  # pi_C_hat ; pi_T_hat
  
  ## sandwich_var returns both variance and sandwich variance
  V_beta_list = sandwich_var_ord(pi_T_hat=pi_T_hat,
                                 pi_C_hat=pi_C_hat,
                                 rho_matrix=rho_mat_new,
                                 cluster_size =cluster_size ,
                                 number_cluster=number_cluster,
                                 y_by_cluster=y_by_cluster,
                                 treat_cluster=treat_cluster)
  V_beta = V_beta_list$var_beta
  sandwich_V_beta = V_beta_list$sandwich_var_beta
  
  if(small_cluster == "yes"){
    t_df = number_cluster- n_cat
    small_corr = number_cluster/t_df
    p <- pt(abs(beta_new/sqrt(diag(sandwich_V_beta)*small_corr) ), df = t_df, lower.tail = F)[n_cat]*2
  }else{
    p <- pnorm(abs(beta_new/sqrt(diag(sandwich_V_beta) )), lower.tail = F)[n_cat]*2
  }
  
  
  if (is.na(p) | is.na(beta_new[n_cat]) |  is.na(V_beta[n_cat,n_cat]) |  is.na(sandwich_V_beta[n_cat,n_cat]) | 
      V_beta[n_cat,n_cat] < 0 |  any(diag(sandwich_V_beta) < 0)
  ){
    warning("p, beta_treat, V_beta_treat, or sandwich_V_beta_treat: NA; or V_beta_treat, or sandwich_V_beta_treat < 0")
    print(c(beta_new[n_cat],V_beta[n_cat,n_cat],sandwich_V_beta[n_cat,n_cat]))
    print(abs(beta_new/sqrt(diag(sandwich_V_beta) )))
    print(sqrt(diag(sandwich_V_beta) ))
    #break
    n_bad <- n_bad + 1L
    next  # go to next dataset in the while() loop
  }
  pvals <- c(pvals, p)
  est_beta_treat = c(est_beta_treat,beta_new[n_cat])
  var_beta_treat = c(var_beta_treat,V_beta[n_cat,n_cat])
  sandwich_var_beta_treat = c(sandwich_var_beta_treat,sandwich_V_beta[n_cat,n_cat])
  est_rho_vec_mat = rbind(est_rho_vec_mat,rho_vec_new)
  #if (is.na(V_beta[3,3])) break
  
  if (length(pvals) %% 100 == 0)
    message(sprintf("  success=%d, tried=%d, bad=%d", length(pvals), n_tried, n_bad))
}  # end while

head(count_iter)
all(converge_ind==T)
head(pvals) 

# --- results ---
power_emp <- mean(pvals < type_1_err)
# power_emp

fail_rate <- n_bad / n_tried

# var_beta_treat[which(is.na(sqrt(var_beta_treat)))] # e.g., -7.2e-05
# var_beta_treat[which(is.na(sqrt(var_beta_treat)))] = 0
test_stat_model_based = est_beta_treat/sqrt(var_beta_treat)
test_stat_sandwich = est_beta_treat/sqrt(sandwich_var_beta_treat)

p_val_noadj = pnorm(abs(test_stat_model_based), lower.tail = F)*2
mean(p_val_noadj < type_1_err)

p_val_noadj_sandwich = pnorm(abs(test_stat_sandwich), lower.tail = F)*2
mean(p_val_noadj_sandwich < type_1_err)

#### t, df= K-2
# t_df_K_2 = number_cluster - 2 #n_cat
# small_corr_K_2 = number_cluster/t_df_K_2
# # pt(abs(beta_new/sqrt(diag(sandwich_V_beta)*small_corr_K_2) ), df = t_df_K_2, lower.tail = F)[3]*2
# 
# p_val_adj_K_2 = pt(abs(test_stat_model_based/sqrt(small_corr_K_2)), df = t_df_K_2, lower.tail = F)*2
# mean(p_val_adj_K_2 < type_1_err)
# 
# p_val_adj_sandwich_K_2 = pt(abs(test_stat_sandwich/sqrt(small_corr_K_2)), df = t_df_K_2, lower.tail = F)*2
# mean(p_val_adj_sandwich_K_2 < type_1_err)

#### t, df= K-n_cat
t_df_K_ncat = number_cluster - n_cat
small_corr_K_ncat = number_cluster/t_df_K_ncat

p_val_adj_K_ncat = pt(abs(test_stat_model_based/sqrt(small_corr_K_ncat)), df = t_df_K_ncat, lower.tail = F)*2
mean(p_val_adj_K_ncat < type_1_err)

p_val_adj_sandwich_K_ncat = pt(abs(test_stat_sandwich/sqrt(small_corr_K_ncat)), df = t_df_K_ncat, lower.tail = F)*2
mean(p_val_adj_sandwich_K_ncat < type_1_err)

mean(est_beta_treat)
mean(var_beta_treat)
mean(sandwich_var_beta_treat)

mean_est_rho_vec = colMeans(est_rho_vec_mat)

rho_vec_to_named_means <- function(rho_vec, prefix = "mean_est_rho") {
  rho_vec <- as.numeric(rho_vec)
  
  m <- length(rho_vec)
  K <- as.integer((sqrt(1 + 8*m) - 1) / 2)
  stopifnot(K * (K + 1) / 2 == m)
  
  idx <- which(upper.tri(matrix(0, K, K), diag = TRUE), arr.ind = TRUE)
  nm  <- paste0(prefix, "_", idx[, 1], "_", idx[, 2])
  
  as.list(setNames(rho_vec, nm))
}
rho_cols <- as.data.frame(rho_vec_to_named_means(mean_est_rho_vec))

power_results <- cbind(data.frame(
  scena      = scena,
  small_cluster = small_cluster,
  Empirical  = power_emp, # adjusted (K-2) for small number of clusters, but non-adjusted for large number of clusters.
  
  Empirical_noadj_model_based = mean(p_val_noadj < type_1_err),
  Empirical_noadj_sandwich = mean(p_val_noadj_sandwich < type_1_err),
  # Empirical_adj_model_based_K_2 = mean(p_val_adj_K_2 < type_1_err),
  # Empirical_adj_sandwich_K_2 = mean(p_val_adj_sandwich_K_2 < type_1_err),
  Empirical_adj_model_based_K_ncat = mean(p_val_adj_K_ncat < type_1_err),
  Empirical_adj_sandwich_K_ncat = mean(p_val_adj_sandwich_K_ncat < type_1_err),
  mean_est_beta_treat = mean(est_beta_treat),
  mean_var_beta_treat = mean(var_beta_treat),
  mean_sandwich_var_beta_treat = mean(sandwich_var_beta_treat),
  n_success   = length(pvals),
  n_tried  =  n_tried,
  n_bad       = n_bad,
  fail_rate   = fail_rate,
  count_iter = mean(count_iter)
), rho_cols
)

result_path = "~/Composite_Endpoints/proj2_setup_011926/4cat/results/mygee_emp_power_051926/"
power_results %>% write_csv(paste0(result_path,"emp_power_update_func",scena,".csv"))


#########################################

###### combine scenario-specific outputs into one output file: 
result_path1 = "~/Composite_Endpoints/proj2_setup_011926/4cat/results/mygee_emp_power_051926/"

ids = gsub("emp_power_update_func([0-9]+).csv","\\1",list.files(result_path1)) %>% as.numeric

scenarios = read_csv("../../results/4cat/emp_power_scenarios_4cat.csv")

num_scena <- nrow(scenarios)

(1:num_scena)[! 1:num_scena  %in% ids]

emp_power = NULL
for (scena in 1:num_scena){
  emp_power= bind_rows(emp_power,read_csv(paste0(result_path1,"emp_power_update_func",scena,".csv")))
}
emp_power

left_join(scenarios,emp_power,by="scena") %>% 
  write_csv("../../results/4cat/emp_power_all_4cat_052426.csv")




