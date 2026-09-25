library(pacman)
p_load(tidyr, dplyr,readr,stringr,magrittr)
setwd("~/Desktop/Dissertation Part2-3/proj2_codes_local022526/4cat/")
out_full <- read_csv("power_comparison_4cat.csv") %>%
  select(
    number_clusters, ordinal_dist_type, latent_ICC, 
    Empirical_adj_model_based_K_ncat, Empirical_noadj_model_based,
    GEE_exch, GEE_bin1, GEE_bin2, WH_DE_latent_icc, #WH_DE_rank_icc, WH
  ) %>%
  mutate(
    Empirical = if_else(
      number_clusters < 40,
      Empirical_adj_model_based_K_ncat,
      Empirical_noadj_model_based
    )
  ) %>%
  transmute(
    `#Clusters`        = number_clusters,
    `Ordinal Dist.`    = factor(ordinal_dist_type,
                                levels = c("bell2","bell", "common", "rare", "U-shape", "even"),
                                labels =c("Bell1", "Bell2","Common", "Rare", "U-shape", "Even"),ordered = T),
    `Latent ICC`       = latent_ICC,
    `Empirical-4cat` = Empirical,
    `GEE-Ordinal-3cat`      = GEE_exch,
    `GEE-Binary1` = GEE_bin1,
    `GEE-Binary2` = GEE_bin2,
    `WH-DE-Latent-3cat`     = WH_DE_latent_icc #,
    #`WH-DE-Rank-3cat`       = WH_DE_rank_icc,
    # `WH-3cat`              = WH
  ) %>% 
  #mutate(across(Empirical:WH),~.x*100) %>% 
  arrange(`#Clusters`,`Ordinal Dist.`,`Latent ICC`)
out_full %>% filter(`#Clusters` %in% c(20,50)) %>% write_csv("table_4cat041526.csv")


est_cols <- c("GEE-Ordinal-3cat","GEE-Binary1","GEE-Binary2","WH-DE-Latent-3cat") #,"WH-DE-Rank-3cat","WH-3cat")

df_fmt <- out_full %>% filter(`#Clusters` %in% c(20,50)) %>%
  mutate(across(c("Empirical-4cat", all_of(est_cols)), ~ round(.x * 100, 1))) %>%  # if your CSV is 0-1
  rowwise() %>%
  mutate(
    min_diff = min(abs(c_across(all_of(est_cols)) - `Empirical-4cat`)),
    across(all_of(est_cols),
           ~ ifelse(abs(.x - `Empirical-4cat`) == min_diff,
                    paste0("\\textbf{", sprintf("%.1f", .x), "}"),
                    sprintf("%.1f", .x)))
  ) %>%
  ungroup()
df_fmt %>% write_csv("table_4cat_fmt041526.csv")


