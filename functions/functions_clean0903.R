################ update to bridge
get_var_bridge=function(phi){
  pi^2 * (phi^(-2) -1) /3
}

# solving for phi
# var_bridge = pi^2 * (phi^(-2) -1) /3
get_phi_bridge = function(var_bridge){
  (var_bridge*3/ pi^2+1 )^(-1/2)
}

get_phi_bridge_from_ICC = function(ICC) sqrt(1-ICC)
# get_phi_bridge_from_ICC(ICC=.2)
# sd_bridge = get_b_sd(ICC=.2)
# sd_bridge^2
# # 0.822
# phi0 = get_phi_bridge(var_bridge=sd_bridge^2)
# phi0
# # 0.894
# get_var_bridge(phi = phi0)
# # 0.822

### power calculation based on equation 10 (ordinal)
get_power_eq10 <- function(log_OR, prop_categories, ICC, cluster_size, m, alpha = 0.05){
  x1 <- 1 - sum(prop_categories^3)
  x2 <- 1 + (cluster_size - 1) * ICC
  x3 <- 6 / (log_OR^2)
  z_alpha_2 <- qnorm(1 - alpha/2)
  pnorm( sqrt(m * x1 / x2 / x3) - z_alpha_2 ) # log_OR should not be smaller than 0?
}

### power calculation based on equation 10, use t-dist (ordinal)
get_power_eq10_small_clus <- function(log_OR, prop_categories, ICC, cluster_size, number_cluster, m, ncat=3,alpha = 0.05){
  x1 <- 1 - sum(prop_categories^3)
  x2 <- 1 + (cluster_size - 1) * ICC
  x3 <- 6 / (log_OR^2)
  t_alpha_2 <- qt(1 - alpha/2,df=number_cluster-ncat) #qnorm(1 - alpha/2)
  pt( sqrt(m * x1 / x2 / x3) - t_alpha_2,df=number_cluster-ncat ) # log_OR should not be smaller than 0?
}

get_power_eq10_small_clus_nct <- function(log_OR, prop_categories, ICC,
                                          cluster_size, number_cluster, m,
                                          ncat = 3, alpha = 0.05) {
  
  x1 <- 1 - sum(prop_categories^3)
  x2 <- 1 + (cluster_size - 1) * ICC
  x3 <- 6 / (log_OR^2)
  
  nu  <- number_cluster - ncat
  tcr <- qt(1 - alpha/2, df = nu)
  
  # noncentrality parameter (positive by construction here)
  ncp <- sqrt(m * x1 / x2 / x3)
  
  # two-sided power: P(|T| > tcr) for T ~ nct(df=nu, ncp=ncp)
  1 - (pt(tcr, df = nu, ncp = ncp) - pt(-tcr, df = nu, ncp = ncp))
}

get_b_sd <- function(ICC) {
  if (ICC <= 0 || ICC >= 1) stop("ICC must be between 0 and 1")
  var_logistic <- pi^2 / 3
  b_var <- (ICC / (1 - ICC)) * var_logistic
  sqrt(b_var)
}

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

# # Marginal probability: integrate over b ~ N(0, sigma_b^2)
# ordinal_logit_marginal_prob <- function(y, x, beta1, cutpoints, epsilon_scale, sigma_b){
#   integrand <- function(b) {
#     ordinal_logit_cond_prob(y, x, beta1, cutpoints, b, epsilon_scale) * dnorm(b, 0, sigma_b)
#   }
#   integrate(integrand, lower = -Inf, upper = Inf, rel.tol = 1e-8)$value
# }

get_marg_lor = function(p0_C, p1_C, p0_T, p1_T){
  cumC_0 <- p0_C;      cumT_0 <- p0_T
  cumC_1 <- p0_C+p1_C; cumT_1 <- p0_T+p1_T
  #cumC_2 <- 1-cumC_1; cumT_2 <- 1- cumT_1
  lor0   <- logit(cumT_0) - logit(cumC_0)
  lor1   <- logit(cumT_1) - logit(cumC_1)
  #lor2   <- logit(cumT_2) - logit(cumC_2)
  mean(c(lor0, lor1), na.rm = TRUE)
}
  
wh_power = function(theta_R,n, pr,A){
  # (est_beta_X = fit$coefficients['X'] %>% as.numeric)
  # n = length(ordinalY_new)
  # pr =VC_ordinal_dist_final
  # theta_R = max(0.01,-est_beta_X)
  sum_cube = sum(pr^3)
  # A is the ratio between control and treatment (randomization: B is treatment, C is placebo)
  # A = as.numeric(table(dt_CE2$randomization)[2]/table(dt_CE2$randomization)[1])
  mu_beta1 = theta_R*sqrt( A*n^3*(1-sum_cube) /3/(n+1)^2/(A+1)^2 ) - qnorm(.975)
  pnorm(mu_beta1)
}

wh_power_small_cluster = function(theta_R,n, pr,A,alpha=.05,
                                  number_cluster,ncat=3){
  # (est_beta_X = fit$coefficients['X'] %>% as.numeric)
  # n = length(ordinalY_new)
  # pr =VC_ordinal_dist_final
  # theta_R = max(0.01,-est_beta_X)
  sum_cube = sum(pr^3)
  t_alpha_2 <- qt(1 - alpha/2,df=number_cluster-ncat) #qnorm(1 - alpha/2)
  # A is the ratio between control and treatment (randomization: B is treatment, C is placebo)
  # A = as.numeric(table(dt_CE2$randomization)[2]/table(dt_CE2$randomization)[1])
  mu_beta1 = theta_R*sqrt( A*n^3*(1-sum_cube) /3/(n+1)^2/(A+1)^2 ) - t_alpha_2 
  pt(mu_beta1,df=number_cluster-ncat )
}

wh_power_small_cluster_nct <- function(theta_R, n, pr, A, alpha = 0.05,
                                       number_cluster, ncat=3) {
  df = number_cluster-ncat
  sum_cube <- sum(pr^3)
  nu  <- df
  tcr <- qt(1 - alpha/2, df = nu)
  
  # noncentrality parameter (signed)
  ncp <- theta_R * sqrt(A * n^3 * (1 - sum_cube) /
                          (3 * (n + 1)^2 * (A + 1)^2))
  
  1 - (pt(tcr, df = nu, ncp = ncp) - pt(-tcr, df = nu, ncp = ncp))
}

wh_V = function(A,n,pr){
  sum_cube = sum(pr^3)
  A*n*(1-sum_cube)/3/(A+1)^2
}


#### Binary power formula
power_binary = function(P0, P1, theta_R, rho,alpha = 0.05,
                        cluster_size,number_cluster,allocate_prop_C=.5){
  DE = 1+ (cluster_size-1)*rho
  var_term = 1 / ( allocate_prop_C*P0*(1-P0) ) + 1 / ( (1-allocate_prop_C)*P1*(1-P1) )
  ex_V = DE/cluster_size * var_term / number_cluster
  z_alpha_2 <- qnorm(1 - alpha/2)
  final_term1_numerator = cluster_size*number_cluster*theta_R^2
  final_term1_denominator = DE*var_term
  final_term1 = sqrt(final_term1_numerator/final_term1_denominator)
  pnorm(final_term1 - z_alpha_2  )
}

power_binary_nct <- function(P0, P1, theta_R, rho,
                           alpha = 0.05,
                           cluster_size, number_cluster,
                           allocate_prop_C = 0.5,
                           df = number_cluster - 2) {
  
  DE <- 1 + (cluster_size - 1) * rho
  
  var_term <- 1/(allocate_prop_C * P0 * (1 - P0)) +
    1/((1 - allocate_prop_C) * P1 * (1 - P1))
  
  # Noncentrality parameter (same "delta" you already compute)
  ncp <- sqrt(cluster_size * number_cluster * theta_R^2 / (DE * var_term))
  
  nu  <- df
  tcr <- qt(1 - alpha/2, df = nu)
  
  # Two-sided power: P(|T| > tcr), T ~ nct(df=nu, ncp=ncp)
  1 - (pt(tcr, df = nu, ncp = ncp) - pt(-tcr, df = nu, ncp = ncp))
  # pt(ncp-tcr, df = nu)
}
#### GEE functions
exch_V_theta = function(pi_T1, pi_T2, pi_C1, pi_C2,
                        rho_mat, # 2 by 2 # if rho_mat is 2 by 2 0, this is independence.
                        # J: cluster size 
                        cluster_size,number_cluster){
  cdf_T2 = pi_T1 + pi_T2
  cdf_C2 = pi_C1 + pi_C2
  
  a_T1 = (pi_T1*(1-pi_T1))^(-1/2); a_T2 = (pi_T2*(1-pi_T2))^(-1/2)
  a_C1 = (pi_C1*(1-pi_C1))^(-1/2); a_C2 = (pi_C2*(1-pi_C2))^(-1/2)
  a_mat_T = matrix(c( a_T1,0,0, a_T2),nrow=2)
  a_mat_C = matrix(c( a_C1,0,0, a_C2),nrow=2)
  
  r_T = -pi_T1* pi_T2 *a_T1*a_T2 #-pi_T1* pi_T2 / sqrt(pi_T1*(1-pi_T1)*pi_T2*(1-pi_T2))
  r_C = -pi_C1* pi_C2 *a_C1*a_C2
  
  Q_mat_T = matrix(c( 1,r_T ,r_T , 1),nrow=2)
  Q_mat_C = matrix(c( 1,r_C ,r_C , 1),nrow=2)  
  
  d_mat_T = matrix(c(pi_T1*(1-pi_T1),-pi_T1*(1-pi_T1),
                     0, cdf_T2*(1-cdf_T2),
                     pi_T1*(1-pi_T1),  cdf_T2*(1-cdf_T2)-pi_T1*(1-pi_T1)),
                   ncol=3) # 2*3
  
  d_mat_C = matrix(c(pi_C1*(1-pi_C1),-pi_C1*(1-pi_C1),
                     0, cdf_C2*(1-cdf_C2),
                     0,  0) ,
                   ncol=3)
  inv_part_T = solve(Q_mat_T+(cluster_size-1)*rho_mat)
  inv_part_C = solve(Q_mat_C+(cluster_size-1)*rho_mat)
  
  DVD_T = cluster_size* t(d_mat_T) %*% a_mat_T %*% inv_part_T %*%a_mat_T %*% d_mat_T
  DVD_C = cluster_size* t(d_mat_C) %*% a_mat_C %*% inv_part_C %*%a_mat_C %*% d_mat_C
  DVD_T
  DVD_C 
  
  #V_mat = 2 * solve(DVD_T+DVD_C)
  V_mat = solve(number_cluster/2*DVD_T+number_cluster/2*DVD_C)
  V_theta = V_mat[3,3]
  V_theta
}

power_gee_est = function(alpha=0.05,V_theta,theta_star){
  z_alpha_2 <- qnorm(1 - alpha/2)
  pnorm(sqrt(theta_star^2/V_theta) - z_alpha_2)
}

power_gee_est_small_cluster = function(alpha=0.05,V_theta,theta_star,
                                       number_cluster,ncat=3){
  # z_alpha_2 <- qnorm(1 - alpha/2)
  # pnorm(sqrt(theta_star^2/V_theta) - z_alpha_2)
  t_alpha_2 <- qt(1 - alpha/2,df=number_cluster-ncat) #qnorm(1 - alpha/2)
  pt( sqrt(theta_star^2/V_theta)- t_alpha_2,df=number_cluster-ncat ) # log_OR should not be smaller than 0?
}

power_gee_est_small_cluster_nct <- function(alpha = 0.05,
                                        V_theta,
                                        theta_star,
                                        number_cluster,
                                        ncat=3) {
  df = number_cluster-ncat
  nu  <- df
  tcr <- qt(1 - alpha/2, df = nu)
  ncp <- theta_star / sqrt(V_theta)
  
  1 - (pt(tcr, df = nu, ncp = ncp) -
         pt(-tcr, df = nu, ncp = ncp))
}

get_ordgee_est_rho_mat = function(gee_fit,sim_data){ 
  # extract
  # real gamma = - intercept
  b <- gee_fit$beta
  alpha <- -unname(b[grepl("^Inter:", names(b))])  # length K-1
  beta_treat <- unname(b["treat1"])
  treat_vec <- model.matrix(~ treat, data = sim_data)[, "treat1"]
  
  # cumulative probs
  eta <- outer(rep(1, nrow(sim_data)), alpha) - treat_vec %o% rep(beta_treat, length(alpha))
  cum <- plogis(eta)                               # n x (K-1)
  
  # category probs (K=3 here)
  pi_hat <- cbind(cum[,1], cum[,2]-cum[,1], 1 - cum[,2])
  pi_hat[c(1,2000),]
  # check
  prop.table(table(sim_data$Y_oppo[sim_data$treat==0]))
  prop.table(table(sim_data$Y_oppo[sim_data$treat==1]))
  # # (optional) clip tiny numerics and renormalize rows
  # pi_hat <- pmax(pi_hat, 0); pi_hat <- pi_hat / rowSums(pi_hat)
  
  # one-hot for observed Y
  Y_onehot <- model.matrix(~ Y_oppo - 1, data = sim_data)  # columns match levels(sim_data$Y)
  
  # Pearson residuals per category
  resid_mat <- (Y_onehot - pi_hat) / sqrt(pi_hat * (1 - pi_hat)) #  Treatment & Control are separated
  head(resid_mat)
  resid_mat_12 = resid_mat[,1:2]
  
  rho_hat <- sim_data |>
    mutate(e1 = resid_mat_12[,1],
           e2 = resid_mat_12[,2]) |>
    group_by(cluster) |>
    summarise(
      Ti = n(),
      rho11_i = (sum(outer(e1, e1, "*")) - sum(e1^2)) / (Ti * (Ti - 1) ),
      rho22_i = (sum(outer(e2, e2, "*")) - sum(e2^2)) / (Ti * (Ti - 1)),
      rho12_i = (sum(outer(e1, e2, "*")) - sum(e1*e2)) / (Ti * (Ti - 1))
    ) |>
    summarise(
      rho11 = mean(rho11_i, na.rm = TRUE),
      rho22 = mean(rho22_i, na.rm = TRUE),
      rho12 = mean(rho12_i, na.rm = TRUE)
    )
  
  head(rho_hat)
  
  est_rho_mat = matrix(c(rho_hat$rho11,rho_hat$rho12,
                         rho_hat$rho12,rho_hat$rho22),ncol=2)
  est_rho_mat 
}
