library(pacman)
p_load(tidyr, dplyr,readr,stringr,magrittr,ggplot2,
       colorspace # for plotting the heatmap
       )
setwd("~/Desktop/Dissertation Part2-3/proj2_codes_local022526")
rho_dt = read_csv("rhomat_from_superpop.csv")
rho_dt

icc_dt = read_csv("rankICCs_from_superpop.csv")
icc_dt

rho_dt_main = rho_dt %>% 
  left_join(icc_dt %>% select(beta_X_list:ordinal_dist_type,mean_rank_icc),
            by = c("beta_X_list","latent_ICC","ordinal_dist_type")) %>%
  mutate(ordinal_dist_type = factor(
    ordinal_dist_type,
    levels = c("bell", "common", "rare", "U-shape","high01low2", "even"),
    labels =c("Bell", "Common", "Rare", "U-shape","High01low2", "Even"),ordered = T
  )) %>%
  select(beta_X_list,latent_ICC,ordinal_dist_type,gamma0,gamma1,mean_rank_icc,rho11:rho12,bin_rho0, bin_rho01 ) %>%
  arrange(beta_X_list,latent_ICC,ordinal_dist_type)
rho_dt_main

rho_dt_main_fmt = rho_dt_main %>% filter(beta_X_list==log(1.5)) %>% 
  transmute( #beta_X= factor(beta_X_list,levels=c(log(1.5),log(2)),
            #                   labels=c("log(1.5)","log(2)")) ,
            `Latent ICC` = latent_ICC ,
            `Ordinal Dist.` = ordinal_dist_type,
            gamma0,gamma1,
            # `Rank ICC` = mean_rank_icc,
            rho11,rho12,rho22,bin_rho01) %>%
  mutate_at(vars(gamma0:bin_rho01),~ sprintf("%.4f", .x))

rho_dt_main_fmt %>% write_csv("from_superpop_3cat_fmt041526.csv")


rho_heatmap = rho_dt_main %>% select(-contains('gamma'),-mean_rank_icc,-contains("bin_rho"))
rho_heatmap 

p_load(forcats)
rho_heatmap <- rho_heatmap %>%
  mutate(
    ordinal_dist_type = fct_recode(
      ordinal_dist_type,
      "Rare-top" = "High01low2"
    )
  )


# rho_heatmap: beta_X_list, latent_ICC, ordinal_dist_type, rho11, rho22, rho12

# assume rho_heatmap exists with columns:
# beta_X_list, latent_ICC, ordinal_dist_type, rho11, rho12, rho22


plot_rho36_by_or <- function(rho_heatmap, OR_target = 1.5,
                             ord_levels = c("Bell","Common","Rare","U-shape","Rare-top","Even"),
                             icc_levels = c(0.01, 0.02, 0.05),
                             label_digits = 4,
                             text_thresh = 0.012,
                             cmax_val = 50) {
  
  # Filter to target OR (beta = log(OR))
  df0 <- rho_heatmap %>%
    mutate(OR = exp(beta_X_list)) %>%
    filter(abs(OR - OR_target) < 1e-8)
  
  stopifnot(nrow(df0) > 0)
  
  # Expand each scenario row into 2x2 cells
  df_panel <- df0 %>%
    mutate(
      ordinal_dist_type = factor(ordinal_dist_type, levels = ord_levels),
      latent_ICC = factor(latent_ICC, levels = icc_levels),
      icc_lab = paste0("rho[latent]==", as.character(latent_ICC))
    ) %>%
    select(ordinal_dist_type, latent_ICC, icc_lab, rho11, rho12, rho22) %>%
    tidyr::crossing(r = 1:2, c = 1:2) %>%
    mutate(
      rho = case_when(
        r == 1 & c == 1 ~ rho11,
        r == 1 & c == 2 ~ rho12,
        r == 2 & c == 1 ~ rho12,
        r == 2 & c == 2 ~ rho22,
        TRUE ~ NA_real_
      ),
      r = factor(r, levels = c(1,2)),
      c = factor(c, levels = c(1,2))
    )
  
  # Global symmetric limits (within this OR block)
  L <- max(abs(df_panel$rho), na.rm = TRUE)
  # text_thresh = 4*L
  ggplot(df_panel, aes(x = c, y = r, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(
      aes(label = sprintf(paste0("%.", label_digits, "f"), rho),
          color = abs(rho) > text_thresh),
      size = 3,
      fontface = "bold",
      na.rm = TRUE
    ) +
    scale_color_manual(values = c("gray25", "white"), guide = "none") +
    scale_fill_continuous_diverging(
      palette = "Purple-Green",
      limits  = c(-L, L),
      name    = expression(rho),
      rev     = FALSE,
      cmax    = cmax_val
    ) +
    scale_y_discrete(limits = rev(levels(df_panel$r))) +  # r=1 on top
    coord_equal() +
    facet_grid(
      icc_lab~ordinal_dist_type,
      labeller = label_parsed,
      switch = "y"
    ) +
    labs(
      title = bquote(
        "Conditional treatment effect: " ~ theta[c] == log(.(OR_target))
      ),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(size = 12),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      strip.placement = "outside",
      strip.text.x = element_text(size = 11, face = "bold"),
      strip.text.y.left = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "grey95")
    )
}

p_or15_hori <- plot_rho36_by_or(rho_heatmap, OR_target = 1.5)
p_or20_hori <- plot_rho36_by_or(rho_heatmap, OR_target = 2)
p_or15_hori
p_or20_hori

p_or15_vert
p_or20_vert




