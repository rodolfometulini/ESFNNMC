# Fixed-effects nuclear-norm matrix completion benchmark
#
# This file contains the R implementation of the fixed-effects nuclear-norm
# matrix completion estimator used as benchmark. The code follows the logic of
# the C++ implementation in Athey et al. (2021), using alternating updates for
# row effects, column effects, and the low-rank component via soft-impute.

update_u <- function(M, mask, L, v) {
  n <- nrow(M)
  u <- numeric(n)
  for (i in seq_len(n)) {
    obs <- which(mask[i, ] == 1)
    if (length(obs) == 0L) {
      u[i] <- 0
    } else {
      b <- L[i, obs] + v[obs] - M[i, obs]
      u[i] <- - mean(b)
    }
  }
  u
}

update_v <- function(M, mask, L, u) {
  p <- ncol(M)
  v <- numeric(p)
  for (j in seq_len(p)) {
    obs <- which(mask[, j] == 1)
    if (length(obs) == 0L) {
      v[j] <- 0
    } else {
      b <- L[obs, j] + u[obs] - M[obs, j]
      v[j] <- - mean(b)
    }
  }
  v
}

# soft-impute update for L (update_L in C++)
soft_impute_update_L <- function(M, mask, L, u, v, lambda_L) {
  # P_omega(M - (u1^T + 1v^T + L)) + L  then SVT with threshold = lambda * |Omega| / 2
  H <- compute_matrix(L, u, v)
  Pomega <- (M - H) * mask
  proj <- Pomega + L
  thr <- lambda_L * sum(mask) / 2
  SVT_reconstruct(proj, thr)
}

### --- Initialization and warm-start path routines --- ###

# initialize u and v with L = 0 and return lambda_L_max (initialize_uv)
initialize_uv <- function(M, mask, to_estimate_u = TRUE, to_estimate_v = TRUE,
                          niter = 1000, rel_tol = 1e-5) {
  n <- nrow(M); p <- ncol(M)
  u <- numeric(n)
  v <- numeric(p)
  L <- matrix(0, n, p)
  obj_old <- compute_objval(M, mask, L, u, v, 0)
  for (it in seq_len(niter)) {
    if (to_estimate_u) u <- update_u(M, mask, L, v) else u <- numeric(n)
    if (to_estimate_v) v <- update_v(M, mask, L, u) else v <- numeric(p)
    obj_new <- compute_objval(M, mask, L, u, v, 0)
    rel_error <- if (obj_old == 0) 0 else (obj_new - obj_old) / obj_old
    obj_old <- obj_new
    if (rel_error < rel_tol && rel_error >= 0) break
  }
  # Compute lambda_L_max as in C++: 2 * max singular value of P_omega / |Omega|
  est <- compute_matrix(L, u, v)
  Pomega <- (M - est) * mask
  svals <- svd(Pomega)$d
  lambda_L_max <- 2 * max(svals) / sum(mask)
  list(u = u, v = v, lambda_L_max = lambda_L_max)
}

# NNM_fit: fit for ONE lambda starting from given L_init, u_init, v_init (NNM_fit)
NNM_fit <- function(M, mask, L_init, u_init, v_init,
                    to_estimate_u = TRUE, to_estimate_v = TRUE,
                    lambda_L, niter = 1000, rel_tol = 1e-5, is_quiet = TRUE) {
  L <- L_init
  u <- u_init
  v <- v_init
  # initial obj uses singular values of L_init
  obj_old <- compute_objval(M, mask, L, u, v, lambda_L)
  for (it in seq_len(niter)) {
    if (to_estimate_u) u <- update_u(M, mask, L, v) else u <- numeric(nrow(M))
    if (to_estimate_v) v <- update_v(M, mask, L, u) else v <- numeric(ncol(M))
    L <- soft_impute_update_L(M, mask, L, u, v, lambda_L)
    obj_new <- compute_objval(M, mask, L, u, v, lambda_L)
    # relative improvement as in C++: (obj_old - obj_new) / obj_old
    rel_error <- if (obj_old == 0) 0 else (obj_old - obj_new) / obj_old
    if (!is_quiet && it %% 50 == 0) message("iter ", it, " obj: ", round(obj_new, 6))
    if (obj_new < 1e-8) break
    if (rel_error < rel_tol && rel_error >= 0) break
    obj_old <- obj_new
  }
  list(L = L, u = u, v = v)
}

# NNM_with_uv_init: warm-start across a decreasing sequence of lambda_L (NNM_with_uv_init)
NNM_with_uv_init <- function(M, mask, u_init, v_init, to_estimate_u = TRUE, to_estimate_v = TRUE,
                             lambda_Ls, niter = 1000, rel_tol = 1e-5, is_quiet = TRUE) {
  num_lam <- length(lambda_Ls)
  n <- nrow(M); p <- ncol(M)
  res <- vector("list", num_lam)
  L_init <- matrix(0, n, p)
  u_cur <- u_init
  v_cur <- v_init
  for (i in seq_len(num_lam)) {
    lam <- lambda_Ls[i]
    fit <- NNM_fit(M, mask, L_init, u_cur, v_cur,
                   to_estimate_u, to_estimate_v, lam, niter, rel_tol, is_quiet)
    res[[i]] <- list(L = fit$L, u = fit$u, v = fit$v, lambda_L = lam)
    # warm start for next lambda
    L_init <- fit$L
    u_cur <- fit$u
    v_cur <- fit$v
  }
  res
}

# NNM wrapper that computes initialize_uv and then returns list of fitted configs
NNM <- function(M, mask, num_lam_L = 100, lambda_L = NULL,
                to_estimate_u = TRUE, to_estimate_v = TRUE,
                niter = 1000, rel_tol = 1e-5, is_quiet = TRUE) {
  tmp_uv <- initialize_uv(M, mask, to_estimate_u, to_estimate_v, niter, rel_tol)
  if (is.null(lambda_L)) {
    max_lam_L <- tmp_uv$lambda_L_max
    if (num_lam_L == 1) lambda_Ls <- 0 else {
      lambda_without_zero <- logsp(log10(max_lam_L), log10(max_lam_L) - 3, num_lam_L - 1)
      lambda_Ls <- c(lambda_without_zero, 0)
    }
  } else {
    lambda_Ls <- lambda_L
    num_lam_L <- length(lambda_Ls)
  }
  if (to_estimate_u || to_estimate_v) {
    tmp_res <- NNM_with_uv_init(M, mask, tmp_uv$u, tmp_uv$v, to_estimate_u, to_estimate_v,
                                lambda_Ls, niter, rel_tol, is_quiet)
  } else {
    tmp_res <- NNM_with_uv_init(M, mask, numeric(nrow(M)), numeric(ncol(M)), to_estimate_u, to_estimate_v,
                                lambda_Ls, niter, rel_tol, is_quiet)
  }
  # return list where each element has L,u,v,lambda_L
  tmp_res
}

### --- Cross-validation helpers --- ###

# create_folds: produce num_folds random training masks (Bernoulli per observed cell) like C++ create_folds
create_folds <- function(M, mask, to_estimate_u = TRUE, to_estimate_v = TRUE,
                         niter = 1000, rel_tol = 1e-5, cv_ratio = 0.8, num_folds = 5) {
  n <- nrow(M); p <- ncol(M)
  out <- vector("list", num_folds)
  for (k in seq_len(num_folds)) {
    ma_new <- matrix(rbinom(n * p, 1, cv_ratio), n, p)
    fold_mask <- mask * ma_new
    M_tr <- M * fold_mask
    tmp_uv <- initialize_uv(M_tr, fold_mask, to_estimate_u, to_estimate_v, niter, rel_tol)
    out[[k]] <- list(u = tmp_uv$u,
                     v = tmp_uv$v,
                     lambda_L_max = tmp_uv$lambda_L_max,
                     fold_mask = fold_mask)
  }
  out
}

### --- Main NNM_CV translation to R --- ###
# This mirrors the logic in C++ NNM_CV
mcnnm_cv_R <- function(M, mask,
                       to_estimate_u = TRUE, to_estimate_v = TRUE,
                       num_lam = 20, niter = 1000, rel_tol = 1e-5,
                       cv_ratio = 0.6, num_folds = 5, is_quiet = TRUE) {
  # input checks
  if (!is.matrix(M) || !is.numeric(M)) stop("M must be numeric matrix")
  if (!is.matrix(mask) || !is.numeric(mask)) stop("mask must be numeric matrix")
  if (!all(dim(M) == dim(mask))) stop("M and mask must match dims")
  if (!all(mask %in% c(0, 1))) stop("mask must be 0/1")
  if (num_lam < 2) stop("num_lam should be >= 2")
  n <- nrow(M); p <- ncol(M)
  confgs <- create_folds(M, mask, to_estimate_u, to_estimate_v, niter, rel_tol, cv_ratio, num_folds)
  # find largest lambda_L_max across folds
  max_lam_L <- max(sapply(confgs, function(x) x$lambda_L_max))
  # create lambda grid (decreasing, last element 0)
  if (num_lam == 1) lambda_Ls <- 0 else {
    lambda_Ls_wo_zero <- logsp(log10(max_lam_L), log10(max_lam_L) - 3, num_lam - 1)
    lambda_Ls <- c(lambda_Ls_wo_zero, 0)
  }
  # matrix to store squared RMSE for each lambda x fold
  MSEmat <- matrix(NA_real_, nrow = num_lam, ncol = num_folds)
  for (k in seq_len(num_folds)) {
    if (!is_quiet) message("Fold number ", k, " started")
    h <- confgs[[k]]
    mask_training <- h$fold_mask
    M_tr <- M * mask_training
    mask_validation <- mask * (1 - mask_training)
    # For training we warm-start from the initial u/v stored in fold
    train_configs <- NNM_with_uv_init(M_tr, mask_training, h$u, h$v, to_estimate_u, to_estimate_v,
                                      lambda_Ls, niter, rel_tol, is_quiet)
    for (i in seq_len(num_lam)) {
      this_config <- train_configs[[i]]
      L_use <- this_config$L
      u_use <- this_config$u
      v_use <- this_config$v
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
  # re-fit on full data: follow C++ approach of using all lambdas >= chosen (warm-start path)
  lambda_Ls_n <- lambda_Ls[lambda_Ls >= best_lambda]  # keeps decreasing larger values up to chosen
  # We'll pass the full lambda_Ls vector but limit to first length(lambda_Ls_n) as C++ did:
  final_configs <- NNM(M, mask, length(lambda_Ls_n), lambda_Ls, to_estimate_u, to_estimate_v, niter, rel_tol, TRUE)
  # mapping: the chosen lambda appears at position minindex among original lambda_Ls,
  # but we only trained the first length(lambda_Ls_n) entries of lambda_Ls; ensure index maps:
  # Find position of best_lambda within the subset used for final fit:
  subset_lambdas <- lambda_Ls[seq_len(length(lambda_Ls_n))]
  pick_pos <- which(abs(subset_lambdas - best_lambda) < .Machine$double.eps^0.5)
  if (length(pick_pos) == 0L) {
    # fallback: take nearest
    pick_pos <- which.min(abs(subset_lambdas - best_lambda))
  }
  z <- final_configs[[pick_pos]]
  L_fin <- z$L
  u_fin <- z$u
  v_fin <- z$v
  list(L = L_fin,
       u = u_fin,
       v = v_fin,
       Avg_RMSE = Avg_RMSE,
       best_lambda = best_lambda,
       min_RMSE = minRMSE,
       lambda_L = lambda_Ls)
}

##### END FUNCTIONS #####

