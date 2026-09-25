## This file contains:
# get_power_eq10()
# wh_power()
# wh_V()

# Archived: get_power_eq10_small_clus()
# Archived: wh_power_small_cluster()


### power calculation based on equation 10 (ordinal)
get_power_eq10 <- function(log_OR, prop_categories, ICC, cluster_size, m, alpha = 0.05){
  x1 <- 1 - sum(prop_categories^3)
  x2 <- 1 + (cluster_size - 1) * ICC
  x3 <- 6 / (log_OR^2)
  z_alpha_2 <- qnorm(1 - alpha/2)
  pnorm( sqrt(m * x1 / x2 / x3) - z_alpha_2 ) # log_OR should not be smaller than 0?
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

wh_V = function(A,n,pr){
  sum_cube = sum(pr^3)
  A*n*(1-sum_cube)/3/(A+1)^2
}


# Archived
### power calculation based on equation 10, use t-dist (ordinal)
get_power_eq10_small_clus <- function(log_OR, prop_categories, ICC, cluster_size, number_cluster, m, ncat=3,alpha = 0.05){
  x1 <- 1 - sum(prop_categories^3)
  x2 <- 1 + (cluster_size - 1) * ICC
  x3 <- 6 / (log_OR^2)
  t_alpha_2 <- qt(1 - alpha/2,df=number_cluster-ncat) #qnorm(1 - alpha/2)
  pt( sqrt(m * x1 / x2 / x3) - t_alpha_2,df=number_cluster-ncat ) # log_OR should not be smaller than 0?
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






