library(pacman)
pacman::p_load(readxl,tidyr, dplyr,readr,stringr,magrittr,purrr,tibble,ggplot2,colorspace)

# read in the file
chtn = read_excel("~/Desktop/Dissertation Part2-3/proj2_new_app_VUMC/chtn_location.xlsx")
source("~/Desktop/Dissertation Part2-3/proj2_codes_local022526/functions/collapse_cats.R")

chtn1 = chtn %>% filter(newOB==1 &  visityear %in% (2021:2024))
chtn1$location %>% table
# 0         Bellevue     Birth Center        Brentwood         Columbia 
# 1               16                1               17                2 
# Cool Springs         Franklin   Hendersonville          Melrose       Northcrest 
# 43                7               32               16                1 
# OHO    Pleasant View           Smyrna Thompson Station 
# 196                2               16               15 
chtn1 %>% filter(location=="0") 

# combine the location based on the geographic closeness:
# o	OHO: 196
# o	Cool Springs: 43
# o	Hendersonville: 32
# o	Bellevue+Pleasant View+Northcrest?: 16+2 = 18
# o	Brentwood+Melrose: 17+16 = 33
# o	Smyrna: 16
# o	Franklin+ Thompson Station+Columbia?: 7+15=22
### exclude 0 and birth center for now

chtn2 = chtn1 %>%
  filter(!location %in% c("0", "Birth Center")) %>%
  mutate(
    cluster = case_when(
      location == "OHO" ~ "OHO",
      location == "Cool Springs" ~ "Cool Springs",
      location == "Hendersonville" ~ "Hendersonville",
      
      location %in% c("Bellevue", "Pleasant View", "Northcrest") ~ "Bellevue/Pleasant View/Northcrest",
      
      location %in% c("Brentwood", "Melrose") ~ "Brentwood/Melrose",
      
      location == "Smyrna" ~ "Smyrna",
      
      location %in% c("Franklin", "Thompson Station", "Columbia") ~ "Franklin/Thompson Station/Columbia",
      
      TRUE ~ location
    ),
    outcome = factor(
      chapoutcome_prim_ch,
      levels = c(
        "None",
        "Placental abruption",
        "Preeclampsia with severe features",
        "Medically indicated preterm birth at <35 weeks",
        "Stillbirth or neonatal death <28 days"
      ),
      ordered = TRUE
    ),
    bin_outcome = chapoutcome_prim_ch_any
  ) %>%
  dplyr::select(outcome,bin_outcome,cluster)


chtn2 %>% write_csv("~/Desktop/Dissertation Part2-3/proj2_new_app_VUMC/chtn2.csv")
table(chtn2$cluster)

table(chtn2$outcome)
prop.table(table(chtn2$outcome))

str(chtn2)
# Medically indicated preterm birth at <35 weeks                                           None 
# 45                                            238 
# Placental abruption              Preeclampsia with severe features 
# 6                                             68 
# Stillbirth or neonatal death <28 days 
# 6 

# estimate of order of severity from 0 to 4 based on degree of permanent outcome. 
# 
# 0. None
# 1. Abruption
# 2. Preeclampsia with severe features
# 3. Medically indicated preterm birth <35 weeks
# 4. Stillbirth, maternal death, or neonatal death <28 days

categories_full = c(
  "None",
  "Placental abruption",
  "Preeclampsia with severe features",
  "Medically indicated preterm birth at <35 weeks",
  "Stillbirth or neonatal death <28 days"
)

pi_all = prop.table(table(chtn2$outcome))
pi_all

source("~/Desktop/Dissertation Part2-3/proj2_codes_local022526/app_stratify/estimate_rho_matrix.R")
source("~/Desktop/Dissertation Part2-3/proj2_codes_local022526/app_stratify/make_pi_df.R")
source("~/Desktop/Dissertation Part2-3/proj2_codes_local022526/functions/design_effect_clustered.R")
source("~/Desktop/Dissertation Part2-3/proj2_codes_local022526/app_stratify/solve_number_cluster.R")


## ---------------------------------------------------------
## Define scenarios (assumes categories_full is ordered as the 5-cat outcome)
## ---------------------------------------------------------
categories_full = c(
  "None",
  "Placental abruption",
  "Preeclampsia with severe features",
  "Medically indicated preterm birth at <35 weeks",
  "Stillbirth or neonatal death <28 days"
)
scenarios <- tibble(
  scenario_id = c(
    "Binary1",
    "Binary2",
    "Binary3",
    "Ordinal-3cat1",
    "Ordinal-3cat2",
    "Ordinal-3cat3",
    "Ordinal-5cat"
  ),
  outcome_type = c(
    "Binary",
    "Binary",
    "Binary",
    "Ordinal-3",
    "Ordinal-3",
    "Ordinal-3",
    "Ordinal-5"
  ),
  
  # collapse_idx refers to indices in categories_full
  # - NULL means keep all original categories (no collapse)
  collapse_idx = list(
    2:5,               # Binary1: merge all events into 1
    list(1:2, 3:5),    # Binary2: More severe event vs none + less severe
    list(1:3, 4:5),    # Binary3: Severe/fatal endpoint vs all others
    list(2:3, 4:5),    # Ordinal-3cat1: none, serious maternal complication, more severe perinatal outcome
    2:4,               # Ordinal-3cat2
    list(1:2,4:5),     # Ordinal-3cat3
    NULL               # Ordinal-5cat: no collapse
  ),
  
  # Human-readable mapping (just for documentation)
  level_map = c(
    "0=None; 1=Any adverse outcome",
    "0=None or placental abruption; 1=More severe events",
    "0=None, placental abruption or preeclampsia; 1=Severe/fatal events",
    "0=None; 1=Placental abruption or preeclampsia with severe features; 2=Medically indicated preterm birth <35 weeks or stillbirth/neonatal death <28 days",
    "0=None; 1=Placental abruption, preeclampsia with severe features, or medically indicated preterm birth <35 weeks; 2=Stillbirth or neonatal death <28 days",
    "0=None or placental abruption; 1=Preeclampsia with severe features; 2=Medically indicated preterm birth <35 weeks or stillbirth/neonatal death <28 days",
    "0=None; 1=Placental abruption; 2=Preeclampsia with severe features; 3=Medically indicated preterm birth <35 weeks; 4=Stillbirth or neonatal death <28 days"
  )
) %>%
  mutate(
    # K after collapsing
    K = map_int(collapse_idx, ~{
      if (is.null(.x)) length(categories_full) else length(collapse_cats(rep(1, length(categories_full)), .x))
    })
  )

scenarios


scenarios_out <- scenarios %>%
  mutate(
    # collapse pi to match scenario
    pi_vec = map(collapse_idx, ~{
      if (is.null(.x)) as.numeric(pi_all) else as.numeric(collapse_cats(p = as.numeric(pi_all), idx = .x))
    }),
    
    # build pi_df (control only)
    pi_df = map(pi_vec, make_pi_df),
    
    # collapse data Y
    data_collapsed = map(collapse_idx, ~{
      if (is.null(.x)) {
        chtn2 %>%
          transmute(
            cluster,
            treat   = factor(0, levels = c(0,1)),
            Y = factor(as.integer(factor(outcome,
                                         levels = categories_full, ordered = TRUE)) - 1L,
                       ordered = TRUE)
          )
      } else {
        chtn2 %>%
          transmute(
            cluster,
            treat   = factor(0, levels = c(0,1)),
            Y = collapse_factor(outcome,
                                idx = .x, levels_ref = categories_full)
          )
      }
    }),
    
    # estimate rho for each scenario
    rho_est_weighted = map2(data_collapsed, pi_df, ~estimate_rho_matrix(data = .x, pi_df = .y,
                                                                        weights = "cluster_size")),
    
    # convenience: show category counts for sanity checks
    Y_counts = map(data_collapsed, ~as.data.frame(table(.x$Y)))
  )

# View scenario summary (without dumping big matrices)
scenarios_out %>%
  select(scenario_id, outcome_type, K, collapse_idx, level_map)

scenarios_out$rho_est_weighted

################################## First, just use empirical rho matrix (most fair comparison)
################################################## application 1: calculate number of clusters for each of settings

## combined from different:
# scenarios of collapsing
# treatment effect size: log(1.25), log(1.5), log(1.75)
# use clusters size: 200
# target power: 80%

scenarios_out1 <- scenarios_out %>%
  mutate(
    n_cat  = map_int(pi_vec, length),
    pi_thr = map(pi_vec, ~ .x[1:(length(.x) - 1)])
  )
scenarios_out1$pi_thr

# use clusters size 200
cluster_size <- 200
treat_lor <- log(c(1.5, 1.75,2))

grid_results <- crossing(
  scenarios_out1 %>% select(scenario_id, n_cat, pi_vec, pi_thr, rho_est_weighted),
  theta = treat_lor
) %>%
  mutate(
    N_out = pmap(
      list(n_cat, pi_vec, pi_thr, rho_est_weighted, theta),
      function(n_cat, pi_vec, pi_thr, rho_mat, theta) {
        
        if (n_cat == 2) {
          # Binary: P0 = Pr(Y=0) (baseline category)
          P0 <- as.numeric(pi_vec[1])
          
          # rho_mat is 1x1; extract scalar
          rho_scalar <- as.numeric(rho_mat[1, 1])
          
          solve_clusters_binary(
            P0 = P0, # control prevalence
            theta_R = theta,
            rho = rho_scalar,
            cluster_size = cluster_size
          )
          
        } else {
          # Ordinal: use your existing solve_number_cluster
          # solve_number_cluster expects pi_C_hat length = (n_cat-1) (your pi_thr)
          solve_number_cluster(
            pi_C_hat = as.numeric(pi_thr),
            theta = theta,
            cluster_size = cluster_size,
            rho_matrix = rho_mat
          )
        }
      }
    ),
    OR = exp(theta),
    N_continuous = map_dbl(N_out, ~ .x$N_required_continuous),
    N_integer    = map_int(N_out, ~ .x$N_required_integer),
    achieved_power = map_dbl(N_out, ~ .x$achieved_power), #,
    
    # unified extraction:
    pi_C_hat = map(N_out, ~ .x$pi_C_hat),
    pi_T_hat = map(N_out, ~ .x$pi_T_hat)
  ) %>%
  select(scenario_id, n_cat, theta,OR, N_continuous, N_integer, achieved_power ,pi_C_hat,pi_T_hat)

grid_results

grid_results %>% select(scenario_id:achieved_power) %>% as.data.frame()


source("~/Desktop/Dissertation Part2-3/proj2_codes_local022526/functions/sandwich_var.R")
assumed_rhos <- c(0.005, 0.01, 0.02 ) #, 0.05)   # up to 3 values
included_methods <- c("Binary1", "Ordinal-3cat1" ) #, "Ordinal-3cat2")


if (length(assumed_rhos) > 3) {
  stop("Please provide no more than 3 assumed rho values.")
}

# =========================================================
# HELPER: build assumed rho matrix from a single scalar rho
# =========================================================
build_rho_assumed <- function(n_cat, rho_entry) {
  if (n_cat == 2) {
    matrix(rho_entry, nrow = 1, ncol = 1)
  } else if (n_cat == 3) {
    matrix(
      c(rho_entry, -rho_entry,
        -rho_entry, rho_entry),
      nrow = 2,
      byrow = TRUE
    )
  } else {
    stop("This design-stage grid is only intended for binary (n_cat = 2) and ordinal-3cat (n_cat = 3).")
  }
}

# =========================================================
# DESIGN-BASED RESULTS FOR USER-SPECIFIED ASSUMED RHOS
# =========================================================

scenarios_design <- scenarios_out1 %>%
  filter(outcome_type %in% c("Binary", "Ordinal-3")) %>%
  select(scenario_id, outcome_type, n_cat, pi_vec, pi_thr)

grid_design <- crossing(
  scenarios_design,
  theta = treat_lor,
  rho_entry = assumed_rhos
) %>%
  mutate(
    OR = exp(theta),
    rho_matrix = pmap(list(n_cat, rho_entry), build_rho_assumed),
    out = pmap(
      list(n_cat, pi_vec, pi_thr, rho_matrix, theta),
      function(n_cat, pi_vec, pi_thr, rho_mat, theta) {
        if (n_cat == 2) {
          P0 <- as.numeric(pi_vec[1])
          rho_scalar <- as.numeric(rho_mat[1, 1])
          
          solve_clusters_binary(
            P0 = P0,
            theta_R = theta,
            rho = rho_scalar,
            cluster_size = cluster_size
          )
        } else {
          solve_number_cluster(
            pi_C_hat = as.numeric(pi_thr),
            theta = theta,
            cluster_size = cluster_size,
            rho_matrix = rho_mat
          )
        }
      }
    ),
    N_continuous   = map_dbl(out, "N_required_continuous"),
    N_integer      = map_int(out, "N_required_integer"),
    achieved_power = map_dbl(out, "achieved_power"),
    pi_C_hat       = map(out, "pi_C_hat"),
    pi_T_hat       = map(out, "pi_T_hat")
  ) %>%
  select(
    scenario_id, outcome_type, n_cat,
    OR, theta, rho_entry,
    N_continuous, N_integer, achieved_power,
    pi_C_hat, pi_T_hat
  )

# =========================================================
# PLOT DATA PREP
# =========================================================

or_levels <- c(1.5, 1.75, 2)
or_labels <- c("log(1.5)", "log(1.75)", "log(2)")

assumed_rho_labels <- paste0("\u03C1 = ", format(assumed_rhos, trim = TRUE, scientific = FALSE))

df_assumed <- grid_design %>%
  transmute(
    scenario_id = as.character(scenario_id),
    OR = factor(OR, levels = or_levels, labels = or_labels),
    N = N_continuous,
    rho_value = rho_entry,
    rho_label = paste0("\u03C1 = ", format(rho_entry, trim = TRUE, scientific = FALSE)),
    rho_source = "Assumed"
  )

df_true <- grid_results %>%
  filter(scenario_id %in% included_methods) %>%
  transmute(
    scenario_id = as.character(scenario_id),
    OR = factor(OR, levels = or_levels, labels = or_labels),
    N = N_continuous,
    rho_value = NA_real_,
    rho_label = "\u03C1 (data)",
    rho_source = "Data-driven"
  )


plot_df <- bind_rows(df_assumed, df_true) %>%
  mutate(
    scenario_id = if_else(is.na(scenario_id) | scenario_id == "NA",
                          "Ordinal-5cat",
                          scenario_id),
    scenario_id = factor(scenario_id, levels = included_methods),
    rho_source = factor(rho_source, levels = c("Assumed", "Data-driven")),
    rho_label = factor(
      rho_label,
      levels = c(assumed_rho_labels, "\u03C1 (data)")
    ),
    color_group = if_else(rho_source == "Data-driven", "\u03C1 (data)", as.character(rho_label))
  ) %>%
  filter(!is.na(scenario_id))

# =========================================================
# COLORS
# data-driven fixed to black
# assumed rho values get up to 3 distinct colors
# =========================================================

base_assumed_cols <- c("#0072B2", "#CC79A7", "#E69F00")  # blue, purple, orange
assumed_color_map <- setNames(base_assumed_cols[seq_along(assumed_rho_labels)], assumed_rho_labels)

color_map <- c(
  assumed_color_map,
  "\u03C1 (data)" = "#000000"
)

# =========================================================
# PLOT FUNCTION
# =========================================================

make_cluster_plot <- function(df,
                              included_methods = levels(df$scenario_id),
                              method_labels = NULL) {
  
  d <- df %>%
    filter(as.character(scenario_id) %in% included_methods) %>%
    mutate(
      scenario_id = factor(as.character(scenario_id), levels = included_methods)
    )
  
  if (is.null(method_labels)) {
    method_labels <- setNames(included_methods, included_methods)
  }
  
  ggplot(
    d,
    aes(
      x = scenario_id,
      y = N,
      color = color_group,
      shape = rho_source
    )
  ) +
    geom_point(
      #position = position_dodge(width = 0.45),
      size = 2.4,
      alpha = 0.9
    ) +
    facet_wrap(~ OR, nrow = 1) +
    scale_color_manual(
      values = color_map,
      breaks = c(assumed_rho_labels, "\u03C1 (data)"),
      name = expression(rho)
    ) +
    scale_shape_manual(
      values = c("Assumed" = 16, "Data-driven" = 17),
      name = expression(rho~source)
    ) +
    scale_x_discrete(labels = method_labels, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, NA),
      breaks = function(x) seq(
        0,
        ceiling(max(x, na.rm = TRUE) / 5) * 5,
        by = 5
      ),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
      x = "Outcome definition",
      y = "Calculated number of clusters (continuous)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      #axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.box = "vertical"
    )
}

# =========================================================
# EXAMPLE CALL
# =========================================================

make_cluster_plot(
  plot_df ,#%>% filter(rho_value==0.01 | is.na(rho_value)),
  included_methods = c("Binary1", "Ordinal-3cat1" ), #, "Ordinal-3cat2"),
  method_labels = c(
    "Binary1" = "Binary",
    "Ordinal-3cat1" = "Ordinal-3cat"#,
    #"Ordinal-3cat2" = "Ordinal-3cat2"
  )
)
