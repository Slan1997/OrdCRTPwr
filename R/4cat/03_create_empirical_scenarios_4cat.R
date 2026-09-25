### generate the scenarios for simulating empirical power
library(dplyr)
library(rms)
library(readr)
library(ordinal)
library(bridgedist)

invisible( lapply( list.files("../functions",
                              pattern = "\\.R$",
                              full.names = TRUE), source ) )

scenarios_otherICC_true = read_csv("../../results/4cat/rankICCs_from_superpop_4cat.csv")

scenarios_rhomat_true = read_csv("../../results/4cat/rhomat_from_superpop_4cat.csv") %>%
  dplyr::select(scena:m,contains("gamma"),contains("rho"))  %>%
  bind_cols( scenarios_otherICC_true %>% dplyr::select(ends_with("icc"),-latent_ICC))

scenarios_rhomat_true
scenarios_addrhomat = tibble(expand.grid(scena = 1:nrow(scenarios_rhomat_true),
                                         cluster_sizes = c(50),
                                         number_clusters = c(10, 20, 50, 100)
))  %>%
  full_join(scenarios_rhomat_true %>%
              dplyr::select(-number_clusters, -cluster_sizes),by="scena") %>%
  # filter(!(number_clusters>40 & (beta_X_list==log(2) | cluster_sizes>50))) %>%
  # filter(!(beta_X_list==log(2) & cluster_sizes>50)) %>%
  arrange(ordinal_dist_type) %>%
  mutate(N = number_clusters*cluster_sizes,
         m = N/2,
         scena = row_number()) %>%
  relocate(scena, .before = 1) %>%
  relocate(rho11:mean_rank_icc,.after=latent_ICC) %>%
  mutate(small_cluster = ifelse(number_clusters<40,"yes","no"))

scenarios_addrhomat %>% write_csv("../../results/4cat/emp_power_scenarios_4cat.csv")
