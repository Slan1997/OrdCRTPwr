# ============================================================
# CHTN application
# ============================================================
#
# This script reproduces the sample-size calculations and
# application figure reported in the manuscript.
#
# The original patient-level pilot data are not publicly
# available. Therefore, this script begins with the derived
# quantities reported in the manuscript: the three-category
# outcome distribution and pilot-data-based within-cluster
# correlation estimates.
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

source("../functions/solve_number_cluster.R")
source("../functions/make_x_mat_helpers.R")

# ------------------------------------------------------------
# 1. Design inputs
# ------------------------------------------------------------

# Three-category ordinal outcome:
#   0 = No event
#   1 = Placental abruption or preeclampsia with severe features
#   2 = Medically indicated preterm birth, stillbirth,
#       maternal death, or neonatal death within 28 days

pi_ordinal <- c(
  0.6556,
  0.2039,
  0.1405
)

# Binary collapse:
#   0 = No event
#   1 = Any adverse event

pi_binary <- c(
  0.6556,
  0.3444
)

# Pilot-data-based within-cluster correlation estimates.
# These are derived quantities from the non-public pilot data
# and are reported in the manuscript.

rho_ordinal_data <- matrix(
  c(
    0.016641160, -0.005457226,
    -0.005457226, -0.004476678
  ),
  nrow = 2,
  byrow = TRUE
)

rho_binary_data <- 0.01651473


# Planned design parameters

cluster_size <- 200

treatment_effects <- log(c(1.5, 1.75, 2))

assumed_rhos <- c(0.005, 0.01, 0.02)


# ------------------------------------------------------------
# 2. Helper function for assumed correlation structures
# ------------------------------------------------------------

# For the binary outcome, rho is scalar.
#
# For the three-category ordinal outcome, the working
# correlation matrix has diagonal entries rho and
# off-diagonal entries -rho.

build_rho_assumed <- function(outcome, rho) {
  
  if (outcome == "Binary") {
    
    matrix(rho, nrow = 1, ncol = 1)
    
  } else if (outcome == "Ordinal-3cat") {
    
    matrix(
      c(
        rho, -rho,
        -rho,  rho
      ),
      nrow = 2,
      byrow = TRUE
    )
    
  } else {
    
    stop("Unknown outcome definition.")
    
  }
}


# ------------------------------------------------------------
# 3. Sample-size calculation under assumed correlations
# ------------------------------------------------------------

design_grid <- crossing(
  outcome = c("Binary", "Ordinal-3cat"),
  theta = treatment_effects,
  rho = assumed_rhos
)

results_assumed <- design_grid %>%
  mutate(
    
    result = pmap(
      list(outcome, theta, rho),
      function(outcome, theta, rho) {
        
        if (outcome == "Binary") {
          
          solve_clusters_binary(
            P0 = pi_binary[1],
            theta_R = theta,
            rho = rho,
            cluster_size = cluster_size
          )
          
        } else {
          
          rho_matrix <- build_rho_assumed(
            outcome = "Ordinal-3cat",
            rho = rho
          )
          
          solve_number_cluster(
            # First K-1 category probabilities
            pi_C_hat = pi_ordinal[1:2],
            theta = theta,
            cluster_size = cluster_size,
            rho_matrix = rho_matrix
          )
        }
      }
    ),
    
    N_continuous = map_dbl(
      result,
      ~ .x$N_required_continuous
    ),
    
    N_integer = map_int(
      result,
      ~ .x$N_required_integer
    ),
    
    OR = exp(theta)
  )


# ------------------------------------------------------------
# 4. Sample-size calculation using pilot-data-based correlations
# ------------------------------------------------------------

results_data <- crossing(
  outcome = c("Binary", "Ordinal-3cat"),
  theta = treatment_effects
) %>%
  mutate(
    
    result = map2(
      outcome,
      theta,
      function(outcome, theta) {
        
        if (outcome == "Binary") {
          
          solve_clusters_binary(
            P0 = pi_binary[1],
            theta_R = theta,
            rho = rho_binary_data,
            cluster_size = cluster_size
          )
          
        } else {
          
          solve_number_cluster(
            pi_C_hat = pi_ordinal[1:2],
            theta = theta,
            cluster_size = cluster_size,
            rho_matrix = rho_ordinal_data
          )
        }
      }
    ),
    
    N_continuous = map_dbl(
      result,
      ~ .x$N_required_continuous
    ),
    
    N_integer = map_int(
      result,
      ~ .x$N_required_integer
    ),
    
    OR = exp(theta)
  )


# ------------------------------------------------------------
# 5. Display calculated sample sizes
# ------------------------------------------------------------

results_assumed %>%
  select(
    outcome,
    OR,
    rho,
    N_continuous,
    N_integer
  ) %>%
  print(n = Inf)

results_data %>%
  select(
    outcome,
    OR,
    N_continuous,
    N_integer
  ) %>%
  print(n = Inf)


# ------------------------------------------------------------
# 6. Prepare Figure 3
# ------------------------------------------------------------

or_levels <- c(1.5, 1.75, 2)

or_labels <- c(
  "log(1.5)",
  "log(1.75)",
  "log(2)"
)

results_assumed_plot <- results_assumed %>%
  transmute(
    outcome,
    OR = factor(
      OR,
      levels = or_levels,
      labels = or_labels
    ),
    N = N_continuous,
    rho_label = paste0(
      "\u03C1 = ",
      format(
        rho,
        trim = TRUE,
        scientific = FALSE
      )
    ),
    rho_source = "Assumed"
  )

results_data_plot <- results_data %>%
  transmute(
    outcome,
    OR = factor(
      OR,
      levels = or_levels,
      labels = or_labels
    ),
    N = N_continuous,
    rho_label = "\u03C1 (data)",
    rho_source = "Data-driven"
  )


# ------------------------------------------------------------
# 7. Reproduce application figure
# ------------------------------------------------------------

# ------------------------------------------------------------
# Colors
# ------------------------------------------------------------

assumed_rho_labels <- paste0(
  "\u03C1 = ",
  format(assumed_rhos, trim = TRUE, scientific = FALSE)
)

base_assumed_cols <- c(
  "#0072B2",  # rho = 0.005
  "#CC79A7",  # rho = 0.010
  "#E69F00"   # rho = 0.020
)

assumed_color_map <- setNames(
  base_assumed_cols[seq_along(assumed_rho_labels)],
  assumed_rho_labels
)

color_map <- c(
  assumed_color_map,
  "\u03C1 (data)" = "#000000"
)

plot_df <- bind_rows(
  results_assumed_plot,
  results_data_plot
) %>%
  mutate(
    outcome = factor(
      outcome,
      levels = c("Binary", "Ordinal-3cat")
    ),
    rho_source = factor(
      rho_source,
      levels = c("Assumed", "Data-driven")
    ),
    rho_label = factor(
      rho_label,
      levels = c(
        assumed_rho_labels,
        "\u03C1 (data)"
      )
    ),
    color_group = as.character(rho_label)
  )

### Figure 3
ggplot(
  plot_df,
  aes(
    x = outcome,
    y = N,
    color = color_group,
    shape = rho_source
  )
) +
  geom_point(
    size = 2.4,
    alpha = 0.9
  ) +
  facet_wrap(
    ~ OR,
    nrow = 1
  ) +
  scale_color_manual(
    values = color_map,
    breaks = c(
      assumed_rho_labels,
      "\u03C1 (data)"
    ),
    name = expression(rho)
  ) +
  scale_shape_manual(
    values = c(
      "Assumed" = 16,
      "Data-driven" = 17
    ),
    name = expression(rho~source)
  ) +
  scale_x_discrete(
    labels = c(
      "Binary" = "Binary",
      "Ordinal-3cat" = "Ordinal-3cat"
    )
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    breaks = function(x) {
      seq(
        0,
        ceiling(max(x, na.rm = TRUE) / 5) * 5,
        by = 5
      )
    },
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = "Outcome definition",
    y = "Calculated number of clusters (continuous)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.box = "vertical"
  )
