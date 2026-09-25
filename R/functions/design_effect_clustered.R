# unequal cluster size design effects
design_effect_clustered <- function(rho,
                                    cluster_sizes = NULL,
                                    mean_cluster_size = NULL,
                                    sd_cluster_size = NULL) {
  # rho: ICC
  # cluster_sizes: numeric vector of cluster sizes
  # mean_cluster_size, sd_cluster_size: summary inputs if cluster_sizes unknown
  
  if (!is.numeric(rho) || length(rho) != 1 || !is.finite(rho)) {
    stop("rho must be a finite numeric scalar.")
  }
  
  # Case 1: full cluster sizes provided
  if (!is.null(cluster_sizes)) {
    if (!is.numeric(cluster_sizes) || length(cluster_sizes) == 0) {
      stop("cluster_sizes must be a non-empty numeric vector.")
    }
    if (any(!is.finite(cluster_sizes)) || any(cluster_sizes <= 0)) {
      stop("All cluster sizes must be finite and > 0.")
    }
    
    c_num <- length(cluster_sizes)
    n_bar <- mean(cluster_sizes)
    sd_n  <- stats::sd(cluster_sizes)
    cv    <- if (n_bar > 0) sd_n / n_bar else NA_real_
    
    de_exact <- (n_bar * c_num) / sum(cluster_sizes / (1 + (cluster_sizes - 1) * rho))
    de_cv <- 1 + (((cv^2 + 1) * n_bar) - 1) * rho
    
    return(list(
      input_type = "cluster_sizes",
      rho = rho,
      n_clusters = c_num,
      mean_cluster_size = n_bar,
      sd_cluster_size = sd_n,
      CV = cv,
      DE_exact = de_exact,
      DE_CV = de_cv
    ))
  }
  
  # Case 2: only mean and sd provided
  if (is.null(mean_cluster_size) || is.null(sd_cluster_size)) {
    stop("Provide either cluster_sizes, or both mean_cluster_size and sd_cluster_size.")
  }
  
  if (!is.numeric(mean_cluster_size) || length(mean_cluster_size) != 1 ||
      !is.finite(mean_cluster_size) || mean_cluster_size <= 0) {
    stop("mean_cluster_size must be a finite numeric scalar > 0.")
  }
  
  if (!is.numeric(sd_cluster_size) || length(sd_cluster_size) != 1 ||
      !is.finite(sd_cluster_size) || sd_cluster_size < 0) {
    stop("sd_cluster_size must be a finite numeric scalar >= 0.")
  }
  
  cv <- sd_cluster_size / mean_cluster_size
  de_cv <- 1 + (((cv^2 + 1) * mean_cluster_size) - 1) * rho
  
  list(
    input_type = "mean_sd_only",
    rho = rho,
    mean_cluster_size = mean_cluster_size,
    sd_cluster_size = sd_cluster_size,
    CV = cv,
    DE_CV = de_cv
  )
}
