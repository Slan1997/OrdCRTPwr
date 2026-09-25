## This file contains:
# get_var_bridge()
# get_phi_bridge()
# get_phi_bridge_from_ICC()
# get_b_sd()

get_phi_bridge_from_ICC = function(ICC) sqrt(1-ICC)

get_var_bridge=function(phi){
  pi^2 * (phi^(-2) -1) /3
}

# solving for phi
# var_bridge = pi^2 * (phi^(-2) -1) /3
get_phi_bridge = function(var_bridge){
  (var_bridge*3/ pi^2+1 )^(-1/2)
}

get_b_sd <- function(ICC) {
  if (ICC <= 0 || ICC >= 1) stop("ICC must be between 0 and 1")
  var_logistic <- pi^2 / 3
  b_var <- (ICC / (1 - ICC)) * var_logistic
  sqrt(b_var)
}

# get_phi_bridge_from_ICC(ICC=.2)
# sd_bridge = get_b_sd(ICC=.2)
# sd_bridge^2
# # 0.822
# phi0 = get_phi_bridge(var_bridge=sd_bridge^2)
# phi0
# # 0.894
# get_var_bridge(phi = phi0)
# # 0.822