# Eigenvector Spatial Filtering Nuclear-Norm Matrix Completion (ESFNNMC)
#
# This file implements the proposed ESFNNMC estimator. Unit fixed effects are
# replaced by a spatially structured component A %*% alpha, where A contains
# selected Moran eigenvectors. The low-rank component is estimated by
# singular-value soft-thresholding along a warm-start lambda path.

compute_matrix_with_A <- function(L, A, alpha, v) {
  u <- as.vector(A %*% alpha)
  L + outer(u, rep(1, length(v))) + outer(rep(1, nrow(L)), v)
}

# Update alpha by weighted least squares over the observed entries.
update_alpha <- function(M, mask, L, v, A) {
  n <- nrow(M)
  q <- ncol(A)
  # Number of observed entries in each row.
  n_i <- rowSums(mask)
  # Mean residual for each spatial unit.
  mean_r <- numeric(n)
  for (i in seq_len(n)) {
    obs <- which(mask[i, ] == 1)
    if (length(obs) == 0L) {
      mean_r[i] <- 0
    } else {
      r_ij <- L[i, obs] + v[obs] - M[i, obs]
      mean_r[i] <- mean(r_ij)
    }
  }
  # Least-squares target.
  y <- - mean_r
  # Use row observation counts as weights and keep rows with at least one observation.
  keep <- which(n_i > 0)
  if (length(keep) == 0L) return(rep(0, q))
  W_sqrt <- sqrt(n_i[keep])
  Aw <- A[keep, , drop = FALSE] * W_sqrt
  yw <- y[keep] * W_sqrt
  # Solve the weighted least-squares problem using QR decomposition.
  qr_Aw <- qr(Aw)
  alpha <- qr.coef(qr_Aw, yw)
  if (any(is.na(alpha))) {
    # Fallback to a small ridge correction if QR returns missing coefficients.
    ridge <- 1e-6
    alpha <- solve(t(Aw) %*% Aw + ridge * diag(ncol(Aw)), t(Aw) %*% yw)
  }
  as.numeric(alpha)
}

# Update column effects using A %*% alpha instead of unrestricted row effects.
update_v_with_A <- function(M, mask, L, A, alpha) {
  p <- ncol(M)
  v <- numeric(p)
  u_vec <- as.vector(A %*% alpha)
  for (j in seq_len(p)) {
    obs <- which(mask[, j] == 1)
    if (length(obs) == 0L) {
      v[j] <- 0
    } else {
      b <- L[obs, j] + u_vec[obs] - M[obs, j]
      v[j] <- - mean(b)
    }
  }
  v
}

# Soft-impute update for the low-rank component with the spatial component fixed.
soft_impute_update_L_with_A <- function(M, mask, L, A, alpha, v, lambda_L) {
  H <- compute_matrix_with_A(L, A, alpha, v)
  Pomega <- (M - H) * mask
  proj <- Pomega + L
  thr <- lambda_L * sum(mask) / 2
  SVT_reconstruct(proj, thr)
}

# Initialize alpha and column effects, and compute the maximum lambda value.
initialize_alpha_v <- function(M, mask, A, to_estimate_alpha = TRUE, to_estimate_v = TRUE,
                               niter = 1000, rel_tol = 1e-5) {
  n <- nrow(M); p <- ncol(M); q <- ncol(A)
  alpha <- numeric(q)
  v <- numeric(p)
  L <- matrix(0, n, p)
  u_vec <- as.vector(A %*% alpha)
  obj_old <- compute_objval(M, mask, L, u_vec, v, 0)
  for (it in seq_len(niter)) {
    if (to_estimate_alpha) alpha <- update_alpha(M, mask, L, v, A) else alpha <- numeric(q)
    if (to_estimate_v) v <- update_v_with_A(M, mask, L, A, alpha) else v <- numeric(p)
    u_vec <- as.vector(A %*% alpha)
    obj_new <- compute_objval(M, mask, L, u_vec, v, 0)
    rel_error <- if (obj_old == 0) 0 else (obj_new - obj_old) / obj_old
    obj_old <- obj_new
    if (rel_error < rel_tol && rel_error >= 0) break
  }
  est <- compute_matrix_with_A(L, A, alpha, v)
  Pomega <- (M - est) * mask
  svals <- svd(Pomega)$d
  lambda_L_max <- 2 * max(svals) / sum(mask)
  list(alpha = alpha, v = v, lambda_L_max = lambda_L_max)
}

# Fit ESFNNMC for a single lambda value.
NNM_fit_with_A <- function(M, mask, L_init, alpha_init, v_init, A,
                           to_estimate_alpha = TRUE, to_estimate_v = TRUE,
                           lambda_L, niter = 1000, rel_tol = 1e-5, is_quiet = TRUE) {
  L <- L_init
  alpha <- alpha_init
  v <- v_init
  u_vec <- as.vector(A %*% alpha)
  obj_old <- compute_objval(M, mask, L, u_vec, v, lambda_L)
  for (it in seq_len(niter)) {
    if (to_estimate_alpha) alpha <- update_alpha(M, mask, L, v, A) else alpha <- numeric(ncol(A))
    if (to_estimate_v) v <- update_v_with_A(M, mask, L, A, alpha) else v <- numeric(ncol(M))
    L <- soft_impute_update_L_with_A(M, mask, L, A, alpha, v, lambda_L)
    u_vec <- as.vector(A %*% alpha)
    obj_new <- compute_objval(M, mask, L, u_vec, v, lambda_L)
    rel_error <- if (obj_old == 0) 0 else (obj_old - obj_new) / obj_old
    if (!is_quiet && it %% 50 == 0) message("iter ", it, " obj: ", round(obj_new, 6))
    if (obj_new < 1e-8) break
    if (rel_error < rel_tol && rel_error >= 0) break
    obj_old <- obj_new
  }
  list(L = L, alpha = alpha, v = v)
}

# Fit ESFNNMC along a warm-start path over a sequence of lambda values.
NNM_with_A_alpha_init <- function(M, mask, A, alpha_init, v_init,
                                  to_estimate_alpha = TRUE, to_estimate_v = TRUE,
                                  lambda_Ls, niter = 1000, rel_tol = 1e-5, is_quiet = TRUE) {
  num_lam <- length(lambda_Ls)
  n <- nrow(M); p <- ncol(M)
  res <- vector("list", num_lam)
  L_init <- matrix(0, n, p)
  alpha_cur <- alpha_init
  v_cur <- v_init
  for (i in seq_len(num_lam)) {
    lam <- lambda_Ls[i]
    fit <- NNM_fit_with_A(M, mask, L_init, alpha_cur, v_cur, A,
                         to_estimate_alpha, to_estimate_v, lam, niter, rel_tol, is_quiet)
    res[[i]] <- list(L = fit$L, alpha = fit$alpha, v = fit$v, lambda_L = lam)
    L_init <- fit$L
    alpha_cur <- fit$alpha
    v_cur <- fit$v
  }
  res
}

# Wrapper for the ESFNNMC warm-start path.
NNM_A <- function(M, mask, A, num_lam_L = 100, lambda_L = NULL,
                  to_estimate_alpha = TRUE, to_estimate_v = TRUE,
                  niter = 1000, rel_tol = 1e-5, is_quiet = TRUE) {
  tmp_av <- initialize_alpha_v(M, mask, A, to_estimate_alpha, to_estimate_v, niter, rel_tol)
  if (is.null(lambda_L)) {
    max_lam_L <- tmp_av$lambda_L_max
    if (num_lam_L == 1) lambda_Ls <- 0 else {
      lambda_without_zero <- logsp(log10(max_lam_L), log10(max_lam_L) - 3, num_lam_L - 1)
      lambda_Ls <- c(lambda_without_zero, 0)
    }
  } else {
    lambda_Ls <- lambda_L
    num_lam_L <- length(lambda_Ls)
  }
  tmp_res <- NNM_with_A_alpha_init(M, mask, A, tmp_av$alpha, tmp_av$v,
                                   to_estimate_alpha, to_estimate_v,
                                   lambda_Ls, niter, rel_tol, is_quiet)
  tmp_res
}

# Create cross-validation folds and store fold-specific alpha and column-effect initial values.
create_folds_with_A <- function(M, mask, A, to_estimate_alpha = TRUE, to_estimate_v = TRUE,
                                niter = 1000, rel_tol = 1e-5, cv_ratio = 0.8, num_folds = 5) {
  n <- nrow(M); p <- ncol(M)
  out <- vector("list", num_folds)
  for (k in seq_len(num_folds)) {
    ma_new <- matrix(rbinom(n * p, 1, cv_ratio), n, p)
    fold_mask <- mask * ma_new
    M_tr <- M * fold_mask
    tmp <- initialize_alpha_v(M_tr, fold_mask, A, to_estimate_alpha, to_estimate_v, niter, rel_tol)
    out[[k]] <- list(alpha = tmp$alpha,
                     v = tmp$v,
                     lambda_L_max = tmp$lambda_L_max,
                     fold_mask = fold_mask)
  }
  out
}

# Main cross-validation routine for ESFNNMC.
mcnnm_cv_R_with_A <- function(M, mask, A,
                              to_estimate_alpha = TRUE, to_estimate_v = TRUE,
                              num_lam = 20, niter = 1000, rel_tol = 1e-5,
                              cv_ratio = 0.6, num_folds = 5, is_quiet = TRUE) {
  # input checks
  if (!is.matrix(M) || !is.numeric(M)) stop("M must be numeric matrix")
  if (!is.matrix(mask) || !is.numeric(mask)) stop("mask must be numeric matrix")
  if (!all(dim(M) == dim(mask))) stop("M and mask must match dims")
  if (!all(mask %in% c(0, 1))) stop("mask must be 0/1")
  if (num_lam < 2) stop("num_lam should be >= 2")
  n <- nrow(M); p <- ncol(M)
  # create folds using A-aware initializer
  confgs <- create_folds_with_A(M, mask, A, to_estimate_alpha, to_estimate_v, niter, rel_tol, cv_ratio, num_folds)
  # find largest lambda_L_max across folds
  max_lam_L <- max(sapply(confgs, function(x) x$lambda_L_max))
  # create lambda grid (decreasing, last element 0)
  if (num_lam == 1) lambda_Ls <- 0 else {
    lambda_Ls_wo_zero <- logsp(log10(max_lam_L), log10(max_lam_L) - 3, num_lam - 1)
    lambda_Ls <- c(lambda_Ls_wo_zero, 0)
  }
  MSEmat <- matrix(NA_real_, nrow = num_lam, ncol = num_folds)
  for (k in seq_len(num_folds)) {
    if (!is_quiet) message("Fold number ", k, " started")
    h <- confgs[[k]]
    mask_training <- h$fold_mask
    M_tr <- M * mask_training
    mask_validation <- mask * (1 - mask_training)
    # train path starting from fold's alpha/v
    train_configs <- NNM_with_A_alpha_init(M_tr, mask_training, A, h$alpha, h$v,
                                          to_estimate_alpha, to_estimate_v, lambda_Ls, niter, rel_tol, is_quiet)
    for (i in seq_len(num_lam)) {
      this_config <- train_configs[[i]]
      L_use <- this_config$L
      alpha_use <- this_config$alpha
      v_use <- this_config$v
      u_use <- as.vector(A %*% alpha_use)
      rmse <- compute_RMSE(M, mask_validation, L_use, u_use, v_use)
      MSEmat[i, k] <- rmse^2
    }
  }
  Avg_MSE <- rowMeans(MSEmat)
  Avg_RMSE <- sqrt(Avg_MSE)
  minindex <- which.min(Avg_RMSE)
  minRMSE <- Avg_RMSE[minindex]
  best_lambda <- lambda_Ls[minindex]
  if (!is_quiet) {
    message("Minimum RMSE achieved on validation set: ", round(minRMSE, 6))
    message("Optimum value of lambda_L: ", best_lambda)
    message("Fitting to the full dataset using optimum lambda_L.")
  }
  # Re-fit on full data: follow approach of training decreasing lambdas >= chosen
  lambda_Ls_n <- lambda_Ls[lambda_Ls >= best_lambda]
  # run NNM_A on full data but limit to first length(lambda_Ls_n) lambdas
  final_configs <- NNM_A(M, mask, A, length(lambda_Ls_n), lambda_Ls, to_estimate_alpha, to_estimate_v, niter, rel_tol, TRUE)
  # find position of best_lambda within trained subset
  subset_lambdas <- lambda_Ls[seq_len(length(lambda_Ls_n))]
  pick_pos <- which(abs(subset_lambdas - best_lambda) < .Machine$double.eps^0.5)
  if (length(pick_pos) == 0L) pick_pos <- which.min(abs(subset_lambdas - best_lambda))
  z <- final_configs[[pick_pos]]
  L_fin <- z$L
  alpha_fin <- z$alpha
  v_fin <- z$v
  u_fin <- as.vector(A %*% alpha_fin)
  list(L = L_fin,
       alpha = alpha_fin,
       u = u_fin,                # also return explicit u = A %*% alpha for convenience
       v = v_fin,
       Avg_RMSE = Avg_RMSE,
       best_lambda = best_lambda,
       min_RMSE = minRMSE,
       lambda_L = lambda_Ls)
}


