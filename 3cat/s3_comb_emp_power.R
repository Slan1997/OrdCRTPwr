result_path1 = "~/Composite_Endpoints/proj2_setup_011926/results/mygee_emp_power_K_2/"

ids = gsub("emp_power_([0-9]+).csv","\\1",list.files(result_path1)) %>% as.numeric

#scenarios = read_csv("~/Composite_Endpoints/proj2_setup_011926/results/emp_power_scenarios.csv")
scenarios = read_csv("~/Composite_Endpoints/proj2_setup_011926/results/emp_power_scenarios_add_ushape020426.csv")

num_scena <- nrow(scenarios)

#sapply(1:10,\(x) all(scenarios[,x]==scenarios_est[,x]))

(1:num_scena)[! 1:num_scena  %in% ids]

emp_power = NULL
for (scena in 1:num_scena){
  emp_power= bind_rows(emp_power,read_csv(paste0(result_path1,"emp_power_",scena,".csv")))
}
emp_power

emp_power %>% filter(if_any(everything(), is.na))

emp_power %>% filter(fail_rate!=0)
emp_power %>% filter(is.na(Empirical))

left_join(scenarios,emp_power,by="scena") %>% filter(Empirical>.98) %>% as.data.frame()

left_join(scenarios,emp_power,by="scena") %>% 
  write_csv("~/Composite_Endpoints/proj2_setup_011926/results/emp_power_all_addushape020626.csv")


###### add lrmrobcov

result_path1 = "~/Composite_Endpoints/proj2_setup_011926/results/lrm_emp_power/"

ids = gsub("emp_power_([0-9]+).csv","\\1",list.files(result_path1)) %>% as.numeric

scenarios = read_csv("~/Composite_Endpoints/proj2_setup_011926/results/emp_power_scenarios.csv")

num_scena <- nrow(scenarios)

#sapply(1:10,\(x) all(scenarios[,x]==scenarios_est[,x]))

(1:num_scena)[! 1:num_scena  %in% ids]

emp_power = NULL
for (scena in 1:num_scena){
  emp_power= bind_rows(emp_power,read_csv(paste0(result_path1,"emp_power_",scena,".csv")))
}
emp_power
 
emp_power %>% filter(NA_ct_lrm>0) #no NA

emp_power_previous = read_csv("~/Composite_Endpoints/proj2_setup_011926/results/emp_power_K_2.csv")


emp_power_addlrm = left_join(emp_power_previous,emp_power %>% 
                               dplyr::select(-NA_ct_lrm),
                             by="scena") %>%
  relocate(Empirical_lrm:Empirical_lrm_small_K_2,.after=Empirical)
  
emp_power_addlrm %>% write_csv("~/Composite_Endpoints/proj2_setup_011926/results/emp_power_K_2_addlrm.csv")
  
  