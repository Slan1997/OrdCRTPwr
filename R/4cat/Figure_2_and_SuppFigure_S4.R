library(pacman)
p_load(tidyr, dplyr,readr,stringr,magrittr,ggplot2)

####### main
out_full_main4 = read_csv("../../results/4cat/power_comparison_4cat_052426.csv") %>%
  dplyr::select(-starts_with("n_"),
                -contains("GEE_indep"),
                -contains("sandwich"),
                -contains("K_2")
  ) 

# out_full = out_full_main4
# n_clusters = 20
# Y_LL=0; Y_UL=1; legend.position = "right"

plot_power_by_nclusters4 <- function(out_full, n_clusters,
                                     Y_LL=0, Y_UL=1,
                                     legend.position = "right") {
  if (n_clusters < 40){
    out <- out_full %>%
      filter(number_clusters == n_clusters) %>%
      dplyr::select(-Empirical,-Empirical_noadj_model_based)%>%
      rename(Empirical = Empirical_adj_model_based_K_ncat)
  }else{
    out <- out_full %>%
      filter(number_clusters == n_clusters) %>%
      dplyr::select(-Empirical,-Empirical_adj_model_based_K_ncat) %>%
      rename(Empirical = Empirical_noadj_model_based)
  }
  
  out_long <- out %>%
    pivot_longer(
      cols = c("WH_DE_latent_icc", "GEE_exch", "GEE_bin1", "GEE_bin2", "Empirical"), #WH_DE_latent_icc:Empirical,
      names_to  = "method",
      values_to = "power"
    )
  
  out_long_plot <- out_long %>%
    dplyr::select(ordinal_dist_type, cluster_sizes, latent_ICC, method, power) %>%
    mutate(
      ## desired panel order left-to-right:
      ## bell, common, rare, high01low2, even
      ordinal_dist_type = factor(
        ordinal_dist_type,
        levels = c("bell2","bell", "common", "rare", "U-shape", "even"),
        labels =c("Bell1", "Bell2","Common", "Rare", "U-shape", "Even")
      ),
      cluster_sizes = factor(cluster_sizes,levels=c(50),labels=c("cluster size: 50")),
      latent_ICC    = factor(latent_ICC,levels=c(0.01,0.02,0.05),
                             labels=c(".01", ".02", ".05")),  # numeric for x-axis
      method        = factor(method, levels=c("Empirical",
                                              "GEE_exch", "GEE_bin1","GEE_bin2",
                                              "WH_DE_latent_icc"#,"WH_DE_rank_icc", "WH" 
      ), 
      labels = c("Empirical-4cat","GEE-Ordinal-3cat", "GEE-Binary1","GEE-Binary2" , "WH-DE-Latent-3cat"#,"WH-DE-Rank", "WH"
      ))
    )
  ## =========================================================
  ## 1) Manual aesthetics
  ## =========================================================
  
  # Colors (as you wrote)
  cols <- c(
    `Empirical-4cat`        = "black",
    `GEE-Ordinal-3cat`    = "#ff7f0e",
    `GEE-Binary1`    = "#9467bd",
    `GEE-Binary2`    = "#c51b8a",
    `WH-DE-Latent-3cat`   = "#78bfc8" #,
    #`WH-DE-Rank`     = "#6fae6f",
    #WH               = "#0050B5"
  )
  
  # Shapes defined PER METHOD, but following FAMILY shapes
  # Empirical ●, all GEE ▲, all WH ■
  shps <- c(
    `Empirical-4cat`        = 16,  # ●
    `GEE-Ordinal-3cat`    = 17,  # ▲
    `GEE-Binary1`    = 17,  # ▲
    `GEE-Binary2`    = 17,  # ▲
    `WH-DE-Latent-3cat`   = 15 #,  # ■
    #`WH-DE-Rank`     = 15,  # ■
    # WH               = 15   # ■
  )
  
  # Linetypes (as you wrote)
  lts <- c(
    `Empirical-4cat`      = "solid",
    `GEE-Ordinal-3cat`  = "solid",
    `GEE-Binary1`  = "longdash",
    `GEE-Binary2`  = "dotdash",
    `WH-DE-Latent-3cat` = "longdash"#,
    #`WH-DE-Rank`   = "longdash",
    #WH             = "longdash"
  )
  
  # Visual hierarchy by class
  lw_vals <- c(Empirical = 1.1, GEE = 0.9, WH = 0.9)
  a_vals  <- c(Empirical = 1.00, GEE = 0.95, WH = 0.95)
  
  ## =========================================================
  ## 2) Create class only (fam no longer needed for shape)
  ## =========================================================
  
  out_long_plot2 <- out_long_plot %>%
    mutate(
      class = case_when(
        method == "Empirical-4cat" ~ "Empirical",
        str_detect(method, "^GEE") ~ "GEE",
        str_detect(method, "^WH")  ~ "WH",
        TRUE ~ "WH"
      )
    )
  
  ## =========================================================
  ## 3) Plot (ONE legend)
  ## =========================================================
  
  ggplot(out_long_plot2,
         aes(x = latent_ICC, y = power * 100,
             group = method,
             color = method,
             shape = method,        # <-- key change (collapse legend)
             linetype = method,
             linewidth = class,
             alpha = class)) +
    
    geom_line() +
    geom_point(size = 2.4, stroke = 0.2) +
    
    # facet_grid(
    #   rows = vars(cluster_sizes),
    #   cols = vars(ordinal_dist_type)
    # ) +
    facet_grid(
      rows = vars(cluster_sizes),
      cols = vars(ordinal_dist_type),
      scales = "free_y"
    ) + 
    scale_y_continuous(
      breaks = seq(Y_LL, Y_UL, 0.05) * 100#,
      #limits = c(Y_LL, Y_UL) * 100
    ) +
    
    scale_color_manual(values = cols, name = "Method") +
    scale_shape_manual(values = shps, name = "Method") +
    
    # keep linetypes but hide their legend so we only have ONE
    scale_linetype_manual(values = lts, guide = "none") +
    
    scale_linewidth_manual(values = lw_vals, guide = "none") +
    scale_alpha_manual(values = a_vals, guide = "none") +
    
    # make the legend keys readable (full opacity/thickness)
    guides(
      color = guide_legend(
        override.aes = list(alpha = 1, linewidth = 1.0)
      )
    ) +
    
    labs(
      # title = paste(n_clusters, "clusters"),
      x = "Latent ICC",
      y = "Power (%)"
    ) +
    
    theme_bw(base_size = 12) +
    theme(
      legend.position = legend.position,
      panel.grid.minor = element_blank()
    )
  
}


#p_n10 
p_n20_4 <- plot_power_by_nclusters4(out_full_main4, #Y_LL=0.1, 
                                    n_clusters = 20,
                                    legend.position = "right")

p_n50_4 <- plot_power_by_nclusters4(out_full_main4, n_clusters = 50,
                                    #Y_LL=.3, Y_UL=1,
                                    legend.position = "right")
### Figure 2
p_n20_4
### Supplementary Figure S4
p_n50_4
