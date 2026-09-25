## most recent one is on local ~/Desktop/Dissertation Part2-3/proj2_codes_local/s4_comb_emp_est.R
library(dplyr)
library(readr)

power_results1 = read_csv("est_power_u_rank_022526.csv")

emp_power = read_csv("emp_power_all_addushape020626.csv") %>%
  dplyr::select(-contains("small_cluster."))
emp_power

final_results1 = full_join(emp_power %>% dplyr::select(scena,Empirical:count_iter),
                           power_results1,
                           by=c("scena" #,
                                #"cluster_sizes",'number_clusters',"beta_X_list","latent_ICC",
                                #"ordinal_dist_type","gamma0","gamma1","N","m"
                           )) %>%
  relocate(Empirical:Empirical_adj_sandwich_K_ncat,.after="GEE_indep") 

final_results1 %>% 
  write_csv("power_comparison_z_u_ricc_022526.csv")
final_results1 %>%filter(is.na(Empirical))

