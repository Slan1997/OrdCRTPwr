result_path1 = "./mygee_emp_power/"

ids = gsub("emp_power_([0-9]+).csv","\\1",list.files(result_path1)) %>% as.numeric

scenarios = read_csv("emp_power_scenarios_4cat.csv")

num_scena <- nrow(scenarios)

(1:num_scena)[! 1:num_scena  %in% ids]

emp_power = NULL
for (scena in 1:num_scena){
  emp_power= bind_rows(emp_power,read_csv(paste0(result_path1,"emp_power_",scena,".csv")))
}
emp_power

left_join(scenarios,emp_power,by="scena") %>% 
  write_csv("emp_power_all_4cat.csv")
