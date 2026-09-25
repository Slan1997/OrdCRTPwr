result_path1 = "~/Desktop/Dissertation Part2-3/proj2_codes_local022526/4cat/superpop_rankICC/"
ids = gsub("superpop_rankICC_([0-9]+).csv","\\1",list.files(result_path1)) %>% as.numeric


#### rankICC
library(dplyr)
library(readr)
result_path = "~/Composite_Endpoints/proj2_setup_011926/results/superpop_cont_ICC/"

scenarios_for_super_pop_mini = read_csv("scenarios_for_superpop_updated_beta1.5_4cat.csv") 
scenarios_for_super_pop_mini = scenarios_for_super_pop_mini%>%
  mutate(number_clusters= rep(500,nrow( scenarios_for_super_pop_mini)),
         N = number_clusters*cluster_sizes,
         m=N/2) 
scenarios_for_super_pop_mini
scenarios_for_super_pop = scenarios_for_super_pop_mini
num_scena <- nrow(scenarios_for_super_pop)
(1:num_scena)[! 1:num_scena  %in% ids]
rankICC_true = NULL
for (scena in 1:num_scena){
  rankICC_true= bind_rows(rankICC_true,read_csv(paste0(result_path1,"superpop_rankICC_",scena,".csv")))
}
rankICC_true
mean(as.numeric(rankICC_true$mean_running_time))
# [1] 350.8994
# final1 = full_join(rankICC_true,rho_mat_true, by="scena")
# final1 


final_out = left_join(scenarios_for_super_pop, rankICC_true #final1
                      , by="scena")
final_out %>% 
  write_csv("rankICCs_from_superpop_4cat.csv")
## this file will be called in the s3_create_emp_scenarios.R
#read_csv("~/Composite_Endpoints/proj2_setup_011926/results/otherICCs_superpop.csv")