# Utilities for collapsing ordinal category probabilities and observed ordinal outcomes.
# Used in the 4+ category simulation/sensitivity analyses.
## =========================================================
## Collapse probabilities p (vector) by merging multiple index blocks
## idx can be:
##   - integer vector (single block): c(2,3,4)
##   - list of integer vectors (multiple blocks): list(2:3, 4:6)
## =========================================================
collapse_cats <- function(p, idx) {
  stopifnot(is.numeric(p), is.vector(p))
  n <- length(p)
  
  # normalize idx input
  if (is.numeric(idx)) idx <- list(idx)
  stopifnot(is.list(idx), length(idx) >= 1)
  
  idx <- lapply(idx, function(v) sort(unique(as.integer(v))))
  stopifnot(all(vapply(idx, function(v) all(v >= 1 & v <= n), logical(1))))
  
  # no overlap allowed (prevents ambiguous merges)
  all_idx <- unlist(idx, use.names = FALSE)
  stopifnot(length(all_idx) == length(unique(all_idx)))
  
  # build group labels: default group = its own index; merged blocks get group = min(block)
  grp <- seq_len(n)
  for (block in idx) grp[block] <- min(block)
  
  # stable ordering by first appearance
  ugrp <- unique(grp)
  
  # sum probabilities within each group
  out <- vapply(ugrp, function(g) sum(p[grp == g]), numeric(1))
  as.numeric(out)
}


## =========================================================
## Collapse observed outcome y (factor/character) into ordered 0:(Knew-1)
## idx same format rules as above
## levels_ref = original ordered category names
## =========================================================
collapse_factor <- function(y, idx, levels_ref) {
  # coerce to ordered factor with specified reference levels
  y <- factor(y, levels = levels_ref, ordered = TRUE)
  n <- length(levels(y))
  
  # normalize idx input
  if (is.numeric(idx)) idx <- list(idx)
  stopifnot(is.list(idx), length(idx) >= 1)
  
  idx <- lapply(idx, function(v) sort(unique(as.integer(v))))
  stopifnot(all(vapply(idx, function(v) all(v >= 1 & v <= n), logical(1))))
  
  # no overlap
  all_idx <- unlist(idx, use.names = FALSE)
  stopifnot(length(all_idx) == length(unique(all_idx)))
  
  # group map old levels -> collapsed groups
  grp <- seq_len(n)
  for (block in idx) grp[block] <- min(block)
  
  # re-rank to consecutive 1..Knew (preserving order)
  grp_new <- match(grp, unique(grp))
  
  # map y to 0-based ordered factor
  y_new <- grp_new[as.integer(y)] - 1L
  factor(y_new, levels = sort(unique(y_new)), ordered = TRUE)
}