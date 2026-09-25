### only keep one Empirical: adj for small number of cluster, no sandwich, df=K-ncat
library(pacman)
p_load(tidyr, dplyr,readr,stringr,magrittr,ggplot2)

# ordinal_dist_type target_pi0 target_pi1 target_pi2
# <chr>                  <dbl>      <dbl>      <dbl>
# 1 even                    0.33       0.33       0.34
# 2 rare                    0.8        0.1        0.1 
# 3 common                  0.1        0.1        0.8 
# 4 bell                    0.1        0.8        0.1 
# 5 high01low2              0.45       0.45       0.1 
# 6 U-shape                 0.45       0.1        0.45

####### main
out_full_main = read_csv("power_comparison_z_u_ricc_022526.csv") %>%
  dplyr::select(-starts_with("n_"),
                -contains("GEE_indep"),
                -contains("sandwich"),
                -contains("K_2")
  ) %>%
  filter(cluster_sizes == 50 & beta_X_list < log(2))

plot_power_by_clusterN <- function(out_full,
                                   fixed_icc = 0.01,
                                   fixed_cluster_size = 50,
                                   Y_LL = 0,
                                   Y_UL = 1,
                                   legend.position = "right") {
  
  out <- out_full %>%
    filter(
      latent_ICC == fixed_icc,
      cluster_sizes == fixed_cluster_size,
      number_clusters %in% c(10, 20, 50, 100)
    ) %>%
    mutate(
      Empirical = ifelse(
        number_clusters < 40,
        Empirical_adj_model_based_K_ncat,
        Empirical_noadj_model_based
      )
    ) %>%
    dplyr::select(
      ordinal_dist_type,
      cluster_sizes,
      latent_ICC,
      number_clusters,
      Empirical,
      GEE_exch,
      GEE_bin1,
      GEE_bin2,
      WH_DE_latent_icc,
      WH_DE_rank_icc,
      WH_DE_rho11_icc,
      WH
    )
  
  out_long <- out %>%
    pivot_longer(
      cols = c(
        Empirical,
        GEE_exch,
        GEE_bin1,
        GEE_bin2,
        WH_DE_latent_icc,
        WH_DE_rank_icc,
        WH_DE_rho11_icc,
        WH
      ),
      names_to = "method",
      values_to = "power"
    )
  
  out_long_plot <- out_long %>%
    mutate(
      ordinal_dist_type = factor(
        ordinal_dist_type,
        levels = c("bell", "common", "rare", "U-shape", "high01low2", "even"),
        labels = c("Bell", "Common", "Rare", "U-shape", "Rare-top", "Even")
      ),
      method = factor(
        method,
        levels = c(
          "Empirical",
          "GEE_exch",
          "GEE_bin1",
          "GEE_bin2",
          "WH_DE_latent_icc",
          "WH_DE_rank_icc",
          "WH_DE_rho11_icc",
          "WH"
        ),
        labels = c(
          "Empirical",
          "GEE-Ordinal",
          "GEE-Binary1",
          "GEE-Binary2",
          "WH-DE-Latent",
          "WH-DE-Rank",
          "WH-DE-rho11",
          "WH"
        )
      ),
      class = case_when(
        method == "Empirical" ~ "Empirical",
        str_detect(method, "^GEE") ~ "GEE",
        str_detect(method, "^WH") ~ "WH",
        TRUE ~ "WH"
      ),
      number_clusters = factor(
        number_clusters,
        levels = c(10, 20, 50, 100)
      )
    )
  
  cols <- c(
    Empirical       = "black",
    `GEE-Ordinal`   = "#ff7f0e",
    `GEE-Binary1`   = "#9467bd",
    `GEE-Binary2`   = "#c51b8a",
    `WH-DE-Latent`  = "#78bfc8",
    #`WH-DE-Rank`    = "#6fae6f",
    `WH-DE-rho11`   = "#0050B5"#,
   #WH              = "#A0A0A0"
  )
  shps <- c(
    Empirical       = 16,
    `GEE-Ordinal`   = 17,
    `GEE-Binary1`   = 17,
    `GEE-Binary2`   = 17,
    `WH-DE-Latent`  = 15,
    #`WH-DE-Rank`    = 15,
    `WH-DE-rho11`   = 15#,
   # WH              = 15
  )
  
  lts <- c(
    Empirical       = "solid",
    `GEE-Ordinal`   = "solid",
    `GEE-Binary1`   = "longdash",
    `GEE-Binary2`   = "dotdash",
    `WH-DE-Latent`  = "longdash",
    #`WH-DE-Rank`    = "longdash",
    `WH-DE-rho11`   = "dotdash"#,
    #WH              = "longdash"
  )
  
  lw_vals <- c(Empirical = 1.1, GEE = 0.9, WH = 0.9)
  a_vals  <- c(Empirical = 1.00, GEE = 0.95, WH = 0.95)
  
  ggplot(
    out_long_plot,
    aes(
      x = number_clusters,
      y = power * 100,
      group = method,
      color = method,
      shape = method,
      linetype = method,
      linewidth = class,
      alpha = class
    )
  ) +
    geom_line() +
    geom_point(size = 2.4, stroke = 0.2) +
    facet_wrap(~ ordinal_dist_type, ncol = 1) +
    scale_y_continuous(
      breaks = seq(Y_LL, Y_UL, 0.1) * 100,
      limits = c(Y_LL, Y_UL) * 100
    ) +
    scale_color_manual(values = cols, name = "Method") +
    scale_shape_manual(values = shps, name = "Method") +
    scale_linetype_manual(values = lts, guide = "none") +
    scale_linewidth_manual(values = lw_vals, guide = "none") +
    scale_alpha_manual(values = a_vals, guide = "none") +
    guides(
      color = guide_legend(
        override.aes = list(alpha = 1, linewidth = 1.0)
      )
    ) +
    labs(
      title = paste0(
        "Power by Number of Clusters (latent ICC = ",
        fixed_icc,
        ", cluster size = ",
        fixed_cluster_size,
        ")"
      ),
      x = "Number of clusters",
      y = "Power (%)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = legend.position,
      panel.grid.minor = element_blank()
    )
}

p_nclusters <- plot_power_by_clusterN(
  out_full = out_full_main,
  fixed_icc = 0.01,
  fixed_cluster_size = 50,
  legend.position = "right"
)

p_nclusters

p_nclusters2 <- plot_power_by_clusterN(
  out_full = out_full_main,
  fixed_icc = 0.02,
  fixed_cluster_size = 50,
  legend.position = "right"
)

p_nclusters2


p_nclusters3 <- plot_power_by_clusterN(
  out_full = out_full_main,
  fixed_icc = 0.05,
  fixed_cluster_size = 50,
  legend.position = "right"
)

p_nclusters3
