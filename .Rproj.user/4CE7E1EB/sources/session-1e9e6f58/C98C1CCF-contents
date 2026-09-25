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
