# Common helper functions for nuclear-norm matrix completion
#
# These routines are shared by the fixed-effects matrix completion benchmark
# and the ESFNNMC estimator. They provide the lambda grid, matrix reconstruction,
# singular-value soft-thresholding, objective evaluation, and RMSE computation.

logsp <- function(start_log, end_log, num_points) {
  if (num_points == 1) return(10^end_log)
  step <- (end_log - start_log) / (num_points - 1)
  10^(start_log + seq(0, num_points - 1) * step)
}

# compute L + u 1^T + 1 v^T  (ComputeMatrix)
compute_matrix <- function(L, u, v) {
  L + outer(u, rep(1, length(v))) + outer(rep(1, length(u)), v)
}

# SVT: soft-threshold singular values and reconstruct (SVT)
SVT_reconstruct <- function(mat, thr) {
  sv <- svd(mat)
  d <- pmax(sv$d - thr, 0)
  # If all truncated to zero, return zero matrix of right dims
  if (all(d == 0)) return(matrix(0, nrow(mat), ncol(mat)))
  sv$u %*% (diag(d, nrow = length(d), ncol = length(d)) %*% t(sv$v))
}

# Compute objective value for L,u,v (Compute_objval in C++)
compute_objval <- function(M, mask, L, u, v, lambda_L) {
  train_size <- sum(mask)
  est <- compute_matrix(L, u, v)
  err <- (est - M) * mask
  sum_sq <- sum(err * err)
  sum_sing <- sum(svd(L)$d)
  (1 / train_size) * sum_sq + lambda_L * sum_sing
}

# RMSE on observed entries (Compute_RMSE)
compute_RMSE <- function(M, mask, L, u, v) {
  valid_size <- sum(mask)
  est <- compute_matrix(L, u, v)
  err <- (est - M) * mask
  sqrt( (1 / valid_size) * sum(err * err) )
}

