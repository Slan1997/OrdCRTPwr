## This file contains:
# exch_V_theta()
# power_gee_est()
# Archived: power_gee_est_small_cluster()

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

# power_gee_est_small_cluster = function(alpha=0.05,V_theta,theta_star,
#                                        number_cluster,ncat=3){
#   # z_alpha_2 <- qnorm(1 - alpha/2)
#   # pnorm(sqrt(theta_star^2/V_theta) - z_alpha_2)
#   t_alpha_2 <- qt(1 - alpha/2,df=number_cluster-ncat) #qnorm(1 - alpha/2)
#   pt( sqrt(theta_star^2/V_theta)- t_alpha_2,df=number_cluster-ncat ) # log_OR should not be smaller than 0?
# }


