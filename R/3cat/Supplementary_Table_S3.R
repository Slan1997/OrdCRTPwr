library(pacman)
p_load(tidyr, dplyr,readr,stringr,magrittr)

out_full <- read_csv("../../results/3cat/power_comparison_z_u_ricc_052426.csv") %>%
  filter(cluster_sizes == 50 & beta_X_list < log(2)) %>%
  select(
    number_clusters, ordinal_dist_type, latent_ICC, mean_rank_icc,
    Empirical_adj_model_based_K_ncat, Empirical_noadj_model_based,
    GEE_exch, GEE_bin1, GEE_bin2,
    WH_DE_latent_icc, WH_DE_rank_icc, WH
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
                                levels = c("bell", "common", "rare", "U-shape","high01low2", "even"),
                                labels =c("Bell", "Common", "Rare", "U-shape","Rare-top", "Even"),ordered = T),
    `Latent ICC`       = latent_ICC,
    `Rank ICC`         = mean_rank_icc,
    Empirical,
    `GEE-Ordinal`      = GEE_exch,
    `GEE-Binary1`      = GEE_bin1,
    `GEE-Binary2`      = GEE_bin2,
    `WH-DE-Latent`     = WH_DE_latent_icc,
    `WH-DE-Rank`       = WH_DE_rank_icc,
    WH               = WH
  ) %>% 
  #mutate(across(Empirical:WH, ~ .x * 100))%>% 
  arrange(`#Clusters`,`Ordinal Dist.`,`Latent ICC`)
#out_full %>% write_csv("table_3cat052426.csv")


est_cols <- c("GEE-Ordinal",#"WH-DE-Rank",
              "GEE-Binary1","GEE-Binary2","WH-DE-Latent"
              #,"WH"
)

df_fmt <- out_full %>%
  mutate(across(c(Empirical, all_of(est_cols)), ~ round(.x * 100, 1))) %>%  # if your CSV is 0-1
  rowwise() %>%
  mutate(
    min_diff = min(abs(c_across(all_of(est_cols)) - Empirical)),
    across(all_of(est_cols),
           ~ ifelse(abs(.x - Empirical) == min_diff,
                    paste0("\\textbf{", sprintf("%.1f", .x), "}"),
                    sprintf("%.1f", .x)))
  ) %>%
  dplyr::select(`#Clusters`,`Ordinal Dist.`,`Latent ICC`,Empirical,all_of(est_cols)) %>%
  filter(`#Clusters` %in% c(20,50)) %>%
  ungroup()
df_fmt %>% write_csv("../../manuscript_tables_and_figures/Supplementary/Supplementary_Table_S3.csv")
