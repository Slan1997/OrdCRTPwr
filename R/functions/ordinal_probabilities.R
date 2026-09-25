## This file contains:
# inv_logit()
# logit()
# ordinal_logit_cond_prob()
# ordinal_logit_marginal_prob_bridge()
# get_marg_lor()

# Logistic CDF (equals plogis(x, location = 0, scale))
inv_logit <- function(x, scale=1){
  1 / (1 + exp(-x/scale))
}

# handy logit
logit <- function(p) log(p/(1 - p))

# Conditional probability for ordinal y ∈ {0,1,...,m} given random effect b
# where m = length(cutpoints)
ordinal_logit_cond_prob <- function(y, x, beta1, cutpoints, epsilon_scale, b) {
  m <- length(cutpoints)          # number of cutpoints
  eta_lin <- beta1 * x + b
  
  if (y < 0 || y > m) stop("Invalid category y")
  
  # cumulative probs P(Y <= j) for j = 0..m-1 correspond to cutpoints[1..m]
  #cum <- vapply(cutpoints, function(cp) inv_logit(cp - eta_lin, epsilon_scale), numeric(1))
  
  if (y == 0) {
    return(inv_logit(cutpoints[1] - eta_lin, epsilon_scale))
  } else if (y == m) {
    return(1 - inv_logit(cutpoints[m] - eta_lin, epsilon_scale))
  } else {
    cum_y = inv_logit(cutpoints[y+1] - eta_lin, epsilon_scale)
    cum_y_prev = inv_logit(cutpoints[y] - eta_lin, epsilon_scale)
    return(cum_y - cum_y_prev)
  }
}

# Marginal probability: integrate over b ~ bridge(phi)
ordinal_logit_marginal_prob_bridge <- function(y, x, beta1, cutpoints, epsilon_scale, b_phi){
  integrand <- function(b) {
    ordinal_logit_cond_prob(y, x, beta1, cutpoints, epsilon_scale, b) * dbridge(b, phi=b_phi)
  }
  integrate(integrand, lower = -Inf, upper = Inf, rel.tol = 1e-8)$value
}

get_marg_lor = function(p0_C, p1_C, p0_T, p1_T){
  cumC_0 <- p0_C;      cumT_0 <- p0_T
  cumC_1 <- p0_C+p1_C; cumT_1 <- p0_T+p1_T
  #cumC_2 <- 1-cumC_1; cumT_2 <- 1- cumT_1
  lor0   <- logit(cumT_0) - logit(cumC_0)
  lor1   <- logit(cumT_1) - logit(cumC_1)
  #lor2   <- logit(cumT_2) - logit(cumC_2)
  mean(c(lor0, lor1), na.rm = TRUE)
}
