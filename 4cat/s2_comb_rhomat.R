result_path = "~/Composite_Endpoints/proj2_setup_011926/4cat/results/superpop_rhomat/"

scenarios_for_super_pop = read_csv("~/Composite_Endpoints/proj2_setup_011926/scenarios_for_superpop_updated_beta1.5_4cat.csv")
num_scena <- nrow(scenarios_for_super_pop)

rho_mat_true = NULL
for (scena in 1:num_scena){
  rho_mat_true= bind_rows(rho_mat_true,read_csv(paste0(result_path,"superpop_rho_mat_addbin",scena,".csv")))
}
rho_mat_true

left_join(scenarios_for_super_pop,rho_mat_true, by="scena")

final_out = left_join(scenarios_for_super_pop,rho_mat_true, by="scena")
final_out %>% 
  write_csv("~/Composite_Endpoints/proj2_setup_011926/4cat/results/rhomat_from_superpop.csv")
