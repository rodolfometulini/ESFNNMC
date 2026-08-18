# Eigenvector Spatial Filtering Nuclear-Norm Matrix Completion (ESFNNMC)
#
# This file implements the proposed ESFNNMC estimator. Unit fixed effects are
# replaced by a spatially structured component A %*% alpha, where A contains
# selected Moran eigenvectors. The low-rank component is estimated by
# singular-value soft-thresholding along a warm-start lambda path.

######### R function translated from C++

# --- keep existing helpers (SVT, compute_RMSE, compute_objval expects u vector) ---
logsp <- function(start_log, end_log, num_points) {
  if (num_points == 1) return(10^end_log)
  step <- (end_log - start_log) / (num_points - 1)
  10^(start_log + seq(0, num_points - 1) * step)
}

compute_matrix <- function(L, u, v) {
  L + outer(u, rep(1, length(v))) + outer(rep(1, length(u)), v)
}

SVT_reconstruct <- function(mat, thr) {
  sv <- svd(mat)
  d <- pmax(sv$d - thr, 0)
  if (all(d == 0)) return(matrix(0, nrow(mat), ncol(mat)))
  sv$u %*% (diag(d, nrow = length(d), ncol = length(d)) %*% t(sv$v))
}

compute_objval <- function(M, mask, L, u, v, lambda_L) {
  train_size <- sum(mask)
  est <- compute_matrix(L, u, v)
  err <- (est - M) * mask
  sum_sq <- sum(err * err)
  sum_sing <- sum(svd(L)$d)
  (1 / train_size) * sum_sq + lambda_L * sum_sing
}

compute_RMSE <- function(M, mask, L, u, v) {
  valid_size <- sum(mask)
  est <- compute_matrix(L, u, v)
  err <- (est - M) * mask
  sqrt((1 / valid_size) * sum(err * err))
}

# --- New helpers using A (spatial eigenvectors) and alpha --- #

# compute matrix using A %*% alpha as the unit effects
compute_matrix_with_A <- function(L, A, alpha, v) {
  u <- as.vector(A %*% alpha)
  L + outer(u, rep(1, length(v))) + outer(rep(1, nrow(L)), v)
}

# update alpha by weighted least squares: minimize sum_i n_i ((A alpha)_i + mean_r_i)^2
update_alpha <- function(M, mask, L, v, A) {
  n <- nrow(M)
  q <- ncol(A)
  # counts per row
  n_i <- rowSums(mask)
  # compute mean residual per row (like original update_u internal b -> mean)
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
  # target y = -mean_r
  y <- - mean_r
  # weights w = n_i; keep only rows with positive weight
  keep <- which(n_i > 0)
  if (length(keep) == 0L) return(rep(0, q))
  W_sqrt <- sqrt(n_i[keep])
  Aw <- A[keep, , drop = FALSE] * W_sqrt
  yw <- y[keep] * W_sqrt
  # solve least squares Aw %*% alpha = yw robustly (qr)
  qr_Aw <- qr(Aw)
  alpha <- qr.coef(qr_Aw, yw)
  if (any(is.na(alpha))) {
    # fallback with small ridge
    ridge <- 1e-6
    alpha <- solve(t(Aw) %*% Aw + ridge * diag(ncol(Aw)), t(Aw) %*% yw)
  }
  as.numeric(alpha)
}

# update v using A %*% alpha instead of u
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

# soft-impute update for L with A/alpha
soft_impute_update_L_with_A <- function(M, mask, L, A, alpha, v, lambda_L) {
  H <- compute_matrix_with_A(L, A, alpha, v)
  Pomega <- (M - H) * mask
  proj <- Pomega + L
  thr <- lambda_L * sum(mask) / 2
  SVT_reconstruct(proj, thr)
}

# initialize alpha and v (replaces initialize_uv)
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

# NNM fit for a single lambda using A/alpha
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

# Warm-start path across lambda sequence with A/alpha
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

# Wrapper NNM equivalent for A: compute initialize_alpha_v then run warm-start path
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

# create_folds for CV storing alpha init
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

# Main CV routine adapted to use A (replaces mcnnm_cv_R)
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


# Select top q eigenvectors by Moran's I
select_spatial_evec80 <- function(W, q = NULL, only_positive = TRUE, threshold = 1e-6) {
  # symmetrize W just to be safe
  W <- (W + t(W)) / 2
  n <- nrow(W)
  H <- diag(n) - matrix(1, n, n) / n
  C <- H %*% W %*% H
  
  eig <- eigen(C, symmetric = TRUE)
  vals <- eig$values
  vecs <- eig$vectors

  # filter eigenvectors by eigenvalue sign
  idx <- which(vals > threshold)
  Afull <- vecs[, idx, drop = FALSE]

  # compute Moran's I for each
  Ivals <- moran_I_evec(Afull, W)

  # keep only positive Moran I if requested
  if (only_positive) {
    keep <- which(Ivals > 0)
    Afull <- Afull[, keep, drop = FALSE]
    Ivals <- Ivals[keep]
  }

  # rank by Moran I magnitude (highest first)
  ord <- order(Ivals, decreasing = TRUE)
  Afull <- Afull[, ord, drop = FALSE]
  Ivals <- Ivals[ord]

  # choose q (default: retain enough to explain ~80% autocorr signal)
  if (is.null(q)) {
    cum <- cumsum(Ivals) / sum(Ivals)
    q <- max(1, which(cum <= 0.80))
  } else {
    q <- min(q, ncol(Afull))
  }

  A <- Afull[, seq_len(q), drop = FALSE]

  list(A = A, MoranI = Ivals[seq_len(q)], idx = idx[seq_len(q)], q = q)
}

select_spatial_evec <- function(W, q = NULL, only_positive = TRUE, threshold = 1e-6) {
  # symmetrize W just to be safe
  W <- (W + t(W)) / 2
  n <- nrow(W)
  H <- diag(n) - matrix(1, n, n) / n
  C <- H %*% W %*% H
  
  eig <- eigen(C, symmetric = TRUE)
  vals <- eig$values
  vecs <- eig$vectors

  # filter eigenvectors by eigenvalue sign
  idx <- which(vals > threshold)
  Afull <- vecs[, idx, drop = FALSE]

  # compute Moran's I for each
  Ivals <- moran_I_evec(Afull, W)

  # keep only positive Moran I if requested
  if (only_positive) {
    keep <- which(Ivals > 0)
    Afull <- Afull[, keep, drop = FALSE]
    Ivals <- Ivals[keep]
  }

  # rank by Moran I magnitude (highest first)
  ord <- order(Ivals, decreasing = TRUE)
  Afull <- Afull[, ord, drop = FALSE]
  Ivals <- Ivals[ord]

  # choose q (default: retain enough to explain ~90% autocorr signal)
  if (is.null(q)) {
    cum <- cumsum(Ivals) / sum(Ivals)
    q <- max(1, which(cum <= 0.90))
  } else {
    q <- min(q, ncol(Afull))
  }

  A <- Afull[, seq_len(q), drop = FALSE]

  list(A = A, MoranI = Ivals[seq_len(q)], idx = idx[seq_len(q)], q = q)
}

select_spatial_evec95 <- function(W, q = NULL, only_positive = TRUE, threshold = 1e-6) {
  # symmetrize W just to be safe
  W <- (W + t(W)) / 2
  n <- nrow(W)
  H <- diag(n) - matrix(1, n, n) / n
  C <- H %*% W %*% H
  
  eig <- eigen(C, symmetric = TRUE)
  vals <- eig$values
  vecs <- eig$vectors

  # filter eigenvectors by eigenvalue sign
  idx <- which(vals > threshold)
  Afull <- vecs[, idx, drop = FALSE]

  # compute Moran's I for each
  Ivals <- moran_I_evec(Afull, W)

  # keep only positive Moran I if requested
  if (only_positive) {
    keep <- which(Ivals > 0)
    Afull <- Afull[, keep, drop = FALSE]
    Ivals <- Ivals[keep]
  }

  # rank by Moran I magnitude (highest first)
  ord <- order(Ivals, decreasing = TRUE)
  Afull <- Afull[, ord, drop = FALSE]
  Ivals <- Ivals[ord]

  # choose q (default: retain enough to explain ~95% autocorr signal)
  if (is.null(q)) {
    cum <- cumsum(Ivals) / sum(Ivals)
    q <- max(1, which(cum <= 0.95))
  } else {
    q <- min(q, ncol(Afull))
  }

  A <- Afull[, seq_len(q), drop = FALSE]

  list(A = A, MoranI = Ivals[seq_len(q)], idx = idx[seq_len(q)], q = q)
}



## compute Moran's I
moran_I_evec <- function(A, W) {
  n <- nrow(A)
  # row-standardize W
  rs <- rowSums(W)
  Wrs <- W / ifelse(rs == 0, 1, rs)

  S0 <- sum(Wrs)

  Ivals <- sapply(seq_len(ncol(A)), function(j) {
    a <- A[, j]
    a <- a - mean(a)
    num <- t(a) %*% Wrs %*% a
    den <- sum(a^2)
    as.numeric((n / S0) * (num / den))
  })
  Ivals
}


### low rank matrix

generate_low_rank <- function(n = 10, rank = 3, singular_scale = c(50, 15, 5), noise_sd = 1e-3, seed = 42) {
  # n: matrix dimension (n x n)
  # rank: target low rank (<< n)
  # singular_scale: numeric vector length rank giving singular values (large -> strong structure)
  # noise_sd: standard deviation of added Gaussian noise (small -> stronger low-rank signal)
  # seed: for reproducibility

  if(length(singular_scale) != rank) {
    # if a single number provided, generate a decaying sequence
    if(length(singular_scale) == 1) {
      singular_scale <- 10^(seq(log10(singular_scale), -1, length.out = rank))
    } else {
      stop("singular_scale must be length 1 or equal to rank")
    }
  }

 # set.seed(seed)
  U <- matrix(rnorm(n * rank), n, rank)
  V <- matrix(rnorm(n * rank), n, rank)

  # orthonormalize columns to make singular values more interpretable
  Uq <- qr.Q(qr(U))
  Vq <- qr.Q(qr(V))

  S <- diag(singular_scale, nrow = rank, ncol = rank)
  M0 <- Uq %*% S %*% t(Vq)          # low-rank part
  noise <- matrix(rnorm(n * n, sd = noise_sd), n, n)
  M <- M0 + noise
  return(list(M = M, low_rank_part = M0, singular_values = svd(M)$d))
}


generate_low_rank_nonsquared <- function(m = 10, n = 8, rank = 3,
                              singular_scale = c(50, 15, 5),
                              noise_sd = 1e-3) {
  # m: number of rows
  # n: number of columns
  # rank: target low rank (<< min(m, n))
  # singular_scale: vector length `rank` giving singular values
  # noise_sd: Gaussian noise SD
  # seed: for reproducibility

  if(rank > min(m, n)) stop("rank must be <= min(m, n)")

  if(length(singular_scale) != rank) {
    # if a single number provided, generate a decaying sequence
    if(length(singular_scale) == 1) {
      singular_scale <- 10^(seq(log10(singular_scale), -1, length.out = rank))
    } else {
      stop("singular_scale must be length 1 or equal to rank")
    }
  }

  U <- matrix(rnorm(m * rank), m, rank)
  V <- matrix(rnorm(n * rank), n, rank)

  # Orthonormalize columns
  Uq <- qr.Q(qr(U))
  Vq <- qr.Q(qr(V))

  S <- diag(singular_scale, nrow = rank, ncol = rank)

  M0 <- Uq %*% S %*% t(Vq)               # low-rank matrix (m x n)
  noise <- matrix(rnorm(m * n, sd = noise_sd), m, n)
  M <- M0 + noise

  return(list(
    M = M,
    low_rank_part = M0,
    singular_values = svd(M)$d
  ))
}


#added 1/04/2026

generate_low_rank_spatiotemp <- function(
  n = 10,
  p = 8,
  rank = 3,
  W,
  rho = 0.3,
  phi = 0.7,
  singular_scale = c(50, 15, 5),
  noise_sd = 1e-3,
  offset = 0,
  enforce_positive = FALSE,
  min_value = 1,
  row_standardize = TRUE,
  check_stability = TRUE,
  time_effect = TRUE,
  regime_breaks = NULL,          # e.g. c(3, 6) gives regimes 1:3, 4:6, 7:p
  regime_means = NULL,           # e.g. c(0, 10, -5)
  regime_sd = 0                  # within-regime variation in column means
) {

  # --- checks ---
  if (!all(dim(W) == c(n, n))) {
    stop("W must be an n x n matrix")
  }

  if (rank > min(n, p)) {
    stop("rank must be <= min(n, p)")
  }

  if (abs(phi) >= 1) {
    stop("phi must be < 1")
  }

  if (length(singular_scale) != rank) {
    stop("singular_scale must have length equal to rank")
  }

  if (min_value <= 0) {
    stop("min_value must be > 0")
  }

  if (!is.null(regime_breaks)) {
    if (any(regime_breaks < 1) || any(regime_breaks >= p)) {
      stop("regime_breaks must be between 1 and p-1")
    }
    if (is.unsorted(regime_breaks)) {
      stop("regime_breaks must be sorted in increasing order")
    }
  }

  # --- row-standardize W ---
  if (row_standardize) {
    rs <- rowSums(W)
    rs[rs == 0] <- 1
    W <- W / rs
  }

  # --- SAR stability ---
  if (check_stability) {
    lambda_max <- max(Mod(eigen(W, only.values = TRUE)$values))
    if (abs(rho) >= 1 / lambda_max) {
      stop("rho violates stability")
    }
  }

  # --- spatial operator ---
  A_sp <- solve(diag(n) - rho * W)

  # --- row scores (spatial units) ---
  A <- matrix(rnorm(n * rank), n, rank)

  # --- latent temporal factors with AR(1) dynamics ---
  B <- matrix(0, p, rank)
  B[1, ] <- rnorm(rank, sd = 1 / sqrt(1 - phi^2))
  for (t in 2:p) {
    B[t, ] <- phi * B[t - 1, ] + rnorm(rank)
  }

  # --- scale factors ---
  S <- diag(singular_scale, rank, rank)

  # --- exact rank-r latent matrix ---
  M0 <- A %*% S %*% t(B)

  # --- apply spatial dependence ---
  M_latent <- A_sp %*% M0

  # --- time effects / regimes ---
  gamma <- rep(0, p)

  if (time_effect) {
    if (is.null(regime_breaks)) {
      # default: one constant level across all columns if regime_means provided,
      # otherwise random column effects
      if (is.null(regime_means)) {
        gamma <- rnorm(p, mean = 0, sd = 5)
      } else if (length(regime_means) == 1) {
        gamma <- rep(regime_means, p)
      } else if (length(regime_means) == p) {
        gamma <- regime_means
      } else {
        stop("Without regime_breaks, regime_means must have length 1 or p")
      }
    } else {
      regime_id <- cut(
        1:p,
        breaks = c(0, regime_breaks, p),
        labels = FALSE,
        right = TRUE
      )
      n_regimes <- length(unique(regime_id))

      if (is.null(regime_means)) {
        regime_means <- seq(0, by = 10, length.out = n_regimes)
      }

      if (length(regime_means) != n_regimes) {
        stop("length(regime_means) must equal number of regimes")
      }

      for (g in 1:n_regimes) {
        idx <- which(regime_id == g)
        gamma[idx] <- regime_means[g] + rnorm(length(idx), sd = regime_sd)
      }
    }
  }

  # --- add time effects ---
  time_component <- matrix(rep(gamma, each = n), nrow = n, ncol = p)

  # --- add observation noise and optional fixed offset ---
  M <- M_latent + time_component + matrix(rnorm(n * p, sd = noise_sd), n, p) + offset

  # --- optional positivity shift for stable MAPE ---
  positivity_shift <- 0
  if (enforce_positive) {
    positivity_shift <- max(0, min_value - min(M))
    M <- M + positivity_shift
  }

  return(list(
    M = M,
    latent_signal = M_latent,
    low_rank_part = M0,
    time_effects = gamma,
    time_component = time_component,
    W = W,
    rho = rho,
    phi = phi,
    spatial_operator = A_sp,
    offset = offset,
    enforce_positive = enforce_positive,
    min_value = min_value,
    positivity_shift = positivity_shift,
    singular_values = svd(M)$d
  ))
}



## 24/06/26:
# generare processo con unit fixed effects al posto di dipendenza spaziale
generate_low_rank_fixed_effects <- function(
  n = 10,
  p = 8,
  rank = 3,
  phi = 0.7,
  singular_scale = c(50, 15, 5),
  unit_fe_sd = 5,
  noise_sd = 1e-3,
  offset = 0,
  enforce_positive = FALSE,
  min_value = 1,
  time_effect = TRUE,
  regime_breaks = NULL,
  regime_means = NULL,
  regime_sd = 0
) {
  
  # --- checks ---
  if (rank > min(n, p)) {
    stop("rank must be <= min(n, p)")
  }
  
  if (abs(phi) >= 1) {
    stop("phi must be < 1")
  }
  
  if (length(singular_scale) != rank) {
    stop("singular_scale must have length equal to rank")
  }
  
  if (min_value <= 0) {
    stop("min_value must be > 0")
  }
  
  if (!is.null(regime_breaks)) {
    if (any(regime_breaks < 1) || any(regime_breaks >= p)) {
      stop("regime_breaks must be between 1 and p-1")
    }
    if (is.unsorted(regime_breaks)) {
      stop("regime_breaks must be sorted in increasing order")
    }
  }
  
  # --- row scores / unit loadings ---
  A <- matrix(rnorm(n * rank), n, rank)
  
  # --- latent temporal factors with AR(1) dynamics ---
  B <- matrix(0, p, rank)
  B[1, ] <- rnorm(rank, sd = 1 / sqrt(1 - phi^2))
  
  for (t in 2:p) {
    B[t, ] <- phi * B[t - 1, ] + rnorm(rank)
  }
  
  # --- scale factors ---
  S <- diag(singular_scale, rank, rank)
  
  # --- exact rank-r latent matrix ---
  M_latent <- A %*% S %*% t(B)
  
  # --- non-spatial unit fixed effects ---
  u <- rnorm(n, mean = 0, sd = unit_fe_sd)
  unit_component <- matrix(rep(u, times = p), nrow = n, ncol = p)
  
  # --- time effects / regimes ---
  gamma <- rep(0, p)
  
  if (time_effect) {
    if (is.null(regime_breaks)) {
      if (is.null(regime_means)) {
        gamma <- rnorm(p, mean = 0, sd = 5)
      } else if (length(regime_means) == 1) {
        gamma <- rep(regime_means, p)
      } else if (length(regime_means) == p) {
        gamma <- regime_means
      } else {
        stop("Without regime_breaks, regime_means must have length 1 or p")
      }
    } else {
      regime_id <- cut(
        1:p,
        breaks = c(0, regime_breaks, p),
        labels = FALSE,
        right = TRUE
      )
      
      n_regimes <- length(unique(regime_id))
      
      if (is.null(regime_means)) {
        regime_means <- seq(0, by = 10, length.out = n_regimes)
      }
      
      if (length(regime_means) != n_regimes) {
        stop("length(regime_means) must equal number of regimes")
      }
      
      for (g in 1:n_regimes) {
        idx <- which(regime_id == g)
        gamma[idx] <- regime_means[g] + rnorm(length(idx), sd = regime_sd)
      }
    }
  }
  
  time_component <- matrix(rep(gamma, each = n), nrow = n, ncol = p)
  
  # --- noise ---
  noise_component <- matrix(rnorm(n * p, mean = 0, sd = noise_sd), n, p)
  
  # --- final matrix ---
  M <- M_latent +
    unit_component +
    time_component +
    noise_component +
    offset
  
  # --- optional positivity shift for stable MAPE ---
  positivity_shift <- 0
  
  if (enforce_positive) {
    positivity_shift <- max(0, min_value - min(M))
    M <- M + positivity_shift
  }
  
  return(list(
    M = M,
    latent_signal = M_latent,
    low_rank_part = M_latent,
    unit_effects = u,
    unit_component = unit_component,
    time_effects = gamma,
    time_component = time_component,
    noise_component = noise_component,
    phi = phi,
    unit_fe_sd = unit_fe_sd,
    offset = offset,
    enforce_positive = enforce_positive,
    min_value = min_value,
    positivity_shift = positivity_shift,
    singular_values = svd(M)$d
  ))
}

##### END FUNCTIONS #####


#### simulare B replicazioni diverse al variare di:
# sparsità di W
# percentuale di missing in M
# struttura low-rank / high-rank
# matrici quadrate, rettangolari con molte righe e rettangolari molte colonne

### Marzo 2026
#### con matrice generate con dipendenza spaziale (in ciascuna colonna) e temporale (in ciascuna riga)
# e con rango ridotto
#### simulare B replicazioni diverse al variare di:
#al variare del numero di missing (sparsity of W: 0.3, phi = 0.5)

# UNIQUE (31/03/26)

rho_vals <- c(0, 0.4, 0.8)
phi_vals <- c(0, 0.4, 0.8)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")

log_file <- "simulation_summary_10_10.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 100
    nrow_mat <- 10
    ncol_mat <- 10
    r <- 3
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = ncol_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(50, 15, 5),
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(3, 6),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) Spatial MC without time fixed effects
        res_sp_ntf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss]) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss]) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_allmodels_rho_", rho_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}

sink()

### END UNIQUE

### 01/04/26: with rank = 5, much harder to estimate (less dominant eigenvalues, less rapid decay)

rho_vals <- c(0, 0.4, 0.8)
phi_vals <- c(0, 0.4, 0.8)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")

log_file <- "simulation_summary_10_10_rank5.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 100
    nrow_mat <- 10
    ncol_mat <- 10
    r <- 5
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = ncol_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(50, 25, 15, 10, 5), # harder to estimate
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(3, 6),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) Spatial MC without time fixed effects
        res_sp_ntf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss]) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss]) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_allmodels_rank5_rho_", rho_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}

sink()

## 24 06 24: simulazione su matrice con effetti fissi e senza dipendenza spaziale
# dove dovrebbe fare meglio FENNMC
# Lo faccio su una 10 x 10 con rank = 5 e 3 regimi (equivalente tabella 1)

unit_fe_sd_vals <- c(0, 5, 15)
phi_vals <- c(0, 0.4, 0.8)

extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")

log_file <- file.path(
  "G:/Il mio Drive/00. PRIN_PNRR_2022/02. Analisi/Spatial matrix completion/Dataset",
  "simulation_summary_10_10_rank5_nonspatial_FE.txt"
)
sink(log_file, split = TRUE)

for (unit_fe_sd in unit_fe_sd_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 100
    nrow_mat <- 10
    ncol_mat <- 10
    r <- 5
    
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)
    
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix used only to construct ESF filters ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = ncol_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data with non-spatial unit fixed effects ---
        lr <- generate_low_rank_fixed_effects(
          n = nrow_mat,
          p = ncol_mat,
          rank = r,
          phi = ar,
          unit_fe_sd = unit_fe_sd,
          noise_sd = 1e-3,
          singular_scale = c(50, 25, 15, 10, 5),
          enforce_positive = TRUE,
          regime_breaks = c(3, 6),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) ESFNNMC with time fixed effects
        res_sp_tf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) ESFNNMC without time fixed effects
        res_sp_ntf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss], na.rm = TRUE) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss], na.rm = TRUE) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss], na.rm = TRUE) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss], na.rm = TRUE) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("unit_fe_sd =", unit_fe_sd, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(
      nfilters,
      miss_seq,
      paste0("Number of spatial filters (ncol(A)) - unit_fe_sd=", unit_fe_sd, ", phi=", ar)
    )
    
    # --- dynamic filename ---
    u_str <- gsub("-", "m", gsub("\\.", "_", as.character(unit_fe_sd)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
#    pdf(paste0("sim_missing_allmodels_rank5_nonspatialFE_unitSD_", u_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (unit FE sd =", unit_fe_sd, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}


sink()


#### 05 06 2026: analisi sensitività, per tau = 0.80

rho_vals <- c(0, 0.4, 0.8)
phi_vals <- c(0, 0.4, 0.8)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")

log_file <- "simulation_summary_10_10_rank5_tau80.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 100
    nrow_mat <- 10
    ncol_mat <- 10
    r <- 5
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = ncol_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(50, 25, 15, 10, 5), # harder to estimate
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(3, 6),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec80(W)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) Spatial MC without time fixed effects
        res_sp_ntf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss]) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss]) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_allmodels_rank5_tau80_rho_", rho_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}

sink()

#### 05 06 2026: analisi sensitività, per tau = 0.95

rho_vals <- c(0, 0.4, 0.8)
phi_vals <- c(0, 0.4, 0.8)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")

log_file <- "simulation_summary_10_10_rank5_tau95.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 100
    nrow_mat <- 10
    ncol_mat <- 10
    r <- 5
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = ncol_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(50, 25, 15, 10, 5), # harder to estimate
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(3, 6),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec95(W)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) Spatial MC without time fixed effects
        res_sp_ntf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss]) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss]) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_allmodels_rank5_tau95_rho_", rho_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}

sink()


#### 02 06 2026: without the two without time FE
rho_vals <- c(0.7)
phi_vals <- c(0.1)
set.seed(150583)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")

log_file <- "simulation_summary_10_10_rank5.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 100
    nrow_mat <- 10
    ncol_mat <- 10
    r <- 5
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)

    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = ncol_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(50, 25, 15, 10, 5), # harder to estimate
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(3, 6),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
               
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
           }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)

    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
  
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_twomodels_rank5_rho_", rho_str, "_phi_", phi_str, ".pdf"), width = 8, height = 6)
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf,
      med_mcfe_tf, med_sp_tf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "darkblue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "darkblue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "darkorange3",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "darkorange3",       lty = 3, lwd = 1.2)
   
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "darkblue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "darkorange3",       lwd = 2, lty = 1)

    legend(
      "topleft",
      legend = c(
        "FENNMC",
        "ESFNNMC",
        "25%-75% band"
      ),
      col = c("darkblue", "darkorange3", "black"),
      lwd = c(2, 2, 1.2),
      lty = c(1, 1, 3),
      pch = c(19, 19, NA),
      bty = "n"
    )
    
    dev.off()
  }
}
sink()

 
# next: 10 righe x 50 colonne (rank = 5)

rho_vals <- c(0, 0.4, 0.8)
phi_vals <- c(0, 0.4, 0.8)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")

log_file <- "simulation_summary_10_50_rank5.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 500
    nrow_mat <- 10
    ncol_mat <- 50
    r <- 5
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = nrow_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(50, 25, 15, 10, 5), # harder to estimate
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(15, 30),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) Spatial MC without time fixed effects
        res_sp_ntf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss]) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss]) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_allmodels_10_50_rank5_rho_", rho_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}

sink()



# next: 10 righe x 50 colonne (rank = 3)

rho_vals <- c(0, 0.4, 0.8)
phi_vals <- c(0, 0.4, 0.8)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")

log_file <- "simulation_summary_10_50_rank3.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 500
    nrow_mat <- 10
    ncol_mat <- 50
    r <- 3
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = nrow_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(50, 15, 5), # harder to estimate
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(15, 30),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) Spatial MC without time fixed effects
        res_sp_ntf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss]) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss]) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_allmodels_10_50_rank3_rho_", rho_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}

sink()



# next: 50 righe x 10 colonne (rank = 5)

rho_vals <- c(0, 0.4, 0.8)
phi_vals <- c(0, 0.4, 0.8)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")

log_file <- "simulation_summary_50_10_rank5.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 500
    nrow_mat <- 50
    ncol_mat <- 10
    r <- 5
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = nrow_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(50, 25, 15, 10, 5),
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(3, 6),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) Spatial MC without time fixed effects
        res_sp_ntf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss]) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss]) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_allmodels_50_10_rank5_rho_", rho_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}

sink()



rho_vals <- c(0.8)
phi_vals <- c(0, 0.4, 0.8)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")

log_file <- "simulation_summary_50_10_rank10_bis.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 500
    nrow_mat <- 50
    ncol_mat <- 10
    r <- 10
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = nrow_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(40, 30, 25, 20, 15, 12, 10, 7, 5, 2),
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(3, 6),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) Spatial MC without time fixed effects
        res_sp_ntf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss]) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss]) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_allmodels_50_10_rank10_rho_", rho_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}

sink()


# 30 x 30 (rank = 10)

rho_vals <- c(0.8)
phi_vals <- c(0, 0.4, 0.8)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")

log_file <- "simulation_summary_30_30_rank10_bis.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 900
    nrow_mat <- 30
    ncol_mat <- 30
    r <- 10
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = nrow_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(40, 30, 25, 20, 15, 12, 10, 7, 5, 2),
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(10, 20),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) Spatial MC without time fixed effects
        res_sp_ntf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss]) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss]) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_allmodels_30_30_rank10_rho_", rho_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}

sink()


# 30 x 30 (rank = 5)

rho_vals <- c(0, 0.4, 0.8)
phi_vals <- c(0, 0.4, 0.8)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")

log_file <- "simulation_summary_30_30_rank5.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 900
    nrow_mat <- 30
    ncol_mat <- 30
    r <- 5
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = nrow_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(40, 25, 15, 10, 5),
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(10, 20),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) Spatial MC without time fixed effects
        res_sp_ntf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss]) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss]) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_allmodels_30_30_rank5_rho_", rho_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}

sink()


# 50 x 50 (rank = 5)

rho_vals <- c(0, 0.4, 0.8)
phi_vals <- c(0, 0.4, 0.8)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")

log_file <- "simulation_summary_50_50_rank5.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 2500
    nrow_mat <- 50
    ncol_mat <- 50
    r <- 5
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = nrow_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(40, 25, 15, 10, 5),
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(15, 30),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- mcnnm_cv_R(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) Spatial MC without time fixed effects
        res_sp_ntf <- mcnnm_cv_R_with_A(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss]) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss]) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_allmodels_50_50_rank5_rho_", rho_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}

sink()



#### 01/04/26: Alternative: choose the number of filters, not the amount of Moran's I



### su matrice generata con AR(phi = 0.5)

# with both time and individual time effect

rho_vals <- c(0, 0.2, 0.4, 0.6, 0.8)

for (rho in rho_vals) {

  start <- Sys.time()
  
  B = 200
  C = length(seq(2,20,2))
  Wspar = 0.3
  n = 100
  col = 10
  row = 10
  ar = 0.5 
  ns = sqrt(n)
  r = 3
  mape_missing = matrix(0, C, B)
  mape_missing_sp = matrix(0, C, B)
  
  set.seed(180326)
  
  for(i in 1:10) {
    for(j in 1:B) {
      
      W = matrix(rbinom(n,1,Wspar),nrow=row, ncol=col)
      W[lower.tri(W)] = t(W)[lower.tri(W)]
      diag(W) = 0
      
      lr = generate_low_rank_spatiotemp(n = row, p = col, rank = r, W = W, rho = rho, phi = ar, noise_sd = 1e-3, singular_scale = c(50, 15, 5), row_standardize = TRUE, check_stability = TRUE)
      mat = lr$M
      
      mask = matrix(rbinom(n,1,1-((i*2)/100)),row,col)
      
      sel = select_spatial_evec(W)
      A = sel$A
      
      res = mcnnm_cv_R(mat, mask, num_lam = 20,
                       to_estimate_u = TRUE, to_estimate_v = TRUE,
                       num_folds = 5, cv_ratio = 0.6,
                       niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred = res$L + matrix(rep(res$u,col),row,col) + 
             t(matrix(rep(res$v,row),col,row))
      
      res_sp = mcnnm_cv_R_with_A(mat, mask, A, num_lam = 20,
                                to_estimate_alpha = TRUE, to_estimate_v = TRUE,
                                num_folds = 5, cv_ratio = 0.6,
                                niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred_sp <- res_sp$L +
        (A %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L)) +
        matrix(1, nrow(res_sp$L), 1) %*% t(res_sp$v)
      
      mape_missing[i,j] = mean(abs((mat - pred) / mat)[mask == 0]) * 100
      mape_missing_sp[i,j] = mean(abs((mat - pred_sp) / mat)[mask == 0]) * 100
    }
  }
  
  end <- Sys.time()
  print(end - start)
  
  # --- summaries ---
  rownames(mape_missing) <- seq(2,20,2)
  colnames(mape_missing) <- 1:B
  rownames(mape_missing_sp) <- seq(2,20,2)
  colnames(mape_missing_sp) <- 1:B
  
  Median <- apply(mape_missing, 1, median, na.rm=TRUE)
  Median_sp <- apply(mape_missing_sp, 1, median, na.rm=TRUE)
  p25 <- apply(mape_missing, 1, quantile, probs = 0.25, na.rm = TRUE)
  p25_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.25, na.rm = TRUE)
  p75 <- apply(mape_missing, 1, quantile, probs = 0.75, na.rm = TRUE)
  p75_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.75, na.rm = TRUE)
  
  # --- dynamic filename ---
  rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
  # e.g. -0.5 -> "m0_5", 0.3 -> "0_3"
  
  setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")
  pdf(paste0("sim_missing_ar_rho_", rho_str, ".pdf"))
  
  x <- as.numeric(rownames(mape_missing))
  
  plot(x, Median, type="n",
       xlab="Percentage of missings", ylab="MAPE",
       main=paste("MAPE (rho =", rho, ")"),
       ylim = c(0,100))
  
  polygon(c(x, rev(x)), c(p25, rev(p75)),
          col = rgb(0, 0, 1, 0.15), border = NA)
  polygon(c(x, rev(x)), c(p25_sp, rev(p75_sp)),
          col = rgb(1, 0, 0, 0.15), border = NA)
  
  lines(x, Median, type="b", pch=19, col="blue", lwd=2)
  lines(x, Median_sp, type="b", pch=19, col="red", lwd=2)
  
  legend("topleft",
         legend = c("MCFE", "ESF-NNMC"),
         col = c("blue", "red"),
         lwd = 2,
         pch = 19,
         bty = "n")
  
  dev.off()
}


# Senza time fixed effects

rho_vals <- c(0, 0.2, 0.4, 0.6, 0.8)

for (rho in rho_vals) {

  start <- Sys.time()
  
  B = 200
  C = length(seq(2,20,2))
  Wspar = 0.3
  n = 100
  col = 10
  row = 10
  ar = 0.5 
  ns = sqrt(n)
  r = 3
  mape_missing = matrix(0, C, B)
  mape_missing_sp = matrix(0, C, B)
  
  set.seed(180326)
  
  for(i in 1:10) {
    for(j in 1:B) {
      
      W = matrix(rbinom(n,1,Wspar),nrow=row, ncol=col)
      W[lower.tri(W)] = t(W)[lower.tri(W)]
      diag(W) = 0
      
      lr = generate_low_rank_spatiotemp(n = row, p = col, rank = r, W = W, rho = rho, phi = ar, noise_sd = 1e-3, singular_scale = c(50, 15, 5), row_standardize = TRUE, check_stability = TRUE)
      mat = lr$M
      
      mask = matrix(rbinom(n,1,1-((i*2)/100)),row,col)
      
      sel = select_spatial_evec(W)
      A = sel$A
      
      res = mcnnm_cv_R(mat, mask, num_lam = 20,
                       to_estimate_u = TRUE, to_estimate_v = FALSE,
                       num_folds = 5, cv_ratio = 0.6,
                       niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred = res$L + matrix(rep(res$u,col),row,col) 
      
      res_sp = mcnnm_cv_R_with_A(mat, mask, A, num_lam = 20,
                                to_estimate_alpha = TRUE, to_estimate_v = FALSE,
                                num_folds = 5, cv_ratio = 0.6,
                                niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred_sp <- res_sp$L +
        (A %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L))
      
      mape_missing[i,j] = mean(abs((mat - pred) / mat)[mask == 0]) * 100
      mape_missing_sp[i,j] = mean(abs((mat - pred_sp) / mat)[mask == 0]) * 100
    }
  }
  
  end <- Sys.time()
  print(end - start)
  
  # --- summaries ---
  rownames(mape_missing) <- seq(2,20,2)
  colnames(mape_missing) <- 1:B
  rownames(mape_missing_sp) <- seq(2,20,2)
  colnames(mape_missing_sp) <- 1:B
  
  Median <- apply(mape_missing, 1, median, na.rm=TRUE)
  Median_sp <- apply(mape_missing_sp, 1, median, na.rm=TRUE)
  p25 <- apply(mape_missing, 1, quantile, probs = 0.25, na.rm = TRUE)
  p25_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.25, na.rm = TRUE)
  p75 <- apply(mape_missing, 1, quantile, probs = 0.75, na.rm = TRUE)
  p75_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.75, na.rm = TRUE)
  
  # --- dynamic filename ---
  rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
  # e.g. -0.5 -> "m0_5", 0.3 -> "0_3"
  
  setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")
  pdf(paste0("sim_missing_notime_ar_rho_", rho_str, ".pdf"))
  
  x <- as.numeric(rownames(mape_missing))
  
  plot(x, Median, type="n",
       xlab="Percentage of missings", ylab="MAPE",
       main=paste("MAPE (rho =", rho, ")"),
       ylim = c(0,100))
  
  polygon(c(x, rev(x)), c(p25, rev(p75)),
          col = rgb(0, 0, 1, 0.15), border = NA)
  polygon(c(x, rev(x)), c(p25_sp, rev(p75_sp)),
          col = rgb(1, 0, 0, 0.15), border = NA)
  
  lines(x, Median, type="b", pch=19, col="blue", lwd=2)
  lines(x, Median_sp, type="b", pch=19, col="red", lwd=2)
  
  legend("topleft",
         legend = c("MCFE", "ESF-NNMC"),
         col = c("blue", "red"),
         lwd = 2,
         pch = 19,
         bty = "n")
  
  dev.off()
}


######### 23 03 2026: 10 righe e 50 colonne 

rho_vals <- c(0, 0.2, 0.4, 0.6, 0.8)

for (rho in rho_vals) {

  start <- Sys.time()
  
  B = 200
  C = length(seq(2,20,2))
  Wspar = 0.3
  n = 500
  col = 50
  row = 10 
  ns = sqrt(n)
  r = 3
  mape_missing = matrix(0, C, B)
  mape_missing_sp = matrix(0, C, B)
  
  set.seed(180326)
  
  for(i in 1:10) {
    for(j in 1:B) {
      
      W = matrix(rbinom(n,1,Wspar),nrow=row, ncol=row)
      W[lower.tri(W)] = t(W)[lower.tri(W)]
      diag(W) = 0
      
      lr = generate_low_rank_spatial(n = row, p = col, rank = r, W = W, rho = rho, noise_sd = 1e-3, singular_scale = c(50, 15, 5), row_standardize = TRUE, check_stability = TRUE)
      mat = lr$M
      
      mask = matrix(rbinom(n,1,1-((i*2)/100)),row,col)
      
      sel = select_spatial_evec(W)
      A = sel$A
      
      res = mcnnm_cv_R(mat, mask, num_lam = 20,
                       to_estimate_u = TRUE, to_estimate_v = TRUE,
                       num_folds = 5, cv_ratio = 0.6,
                       niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred = res$L + matrix(rep(res$u,col),row,col) + 
             t(matrix(rep(res$v,row),col,row))
      
      res_sp = mcnnm_cv_R_with_A(mat, mask, A, num_lam = 20,
                                to_estimate_alpha = TRUE, to_estimate_v = TRUE,
                                num_folds = 5, cv_ratio = 0.6,
                                niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred_sp <- res_sp$L +
        (A %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L)) +
        matrix(1, nrow(res_sp$L), 1) %*% t(res_sp$v)
      
      mape_missing[i,j] = mean(abs((mat - pred) / mat)[mask == 0]) * 100
      mape_missing_sp[i,j] = mean(abs((mat - pred_sp) / mat)[mask == 0]) * 100
    }
  }
  
  end <- Sys.time()
  print(end - start)
  
  # --- summaries ---
  rownames(mape_missing) <- seq(2,20,2)
  colnames(mape_missing) <- 1:B
  rownames(mape_missing_sp) <- seq(2,20,2)
  colnames(mape_missing_sp) <- 1:B
  
  Median <- apply(mape_missing, 1, median, na.rm=TRUE)
  Median_sp <- apply(mape_missing_sp, 1, median, na.rm=TRUE)
  p25 <- apply(mape_missing, 1, quantile, probs = 0.25, na.rm = TRUE)
  p25_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.25, na.rm = TRUE)
  p75 <- apply(mape_missing, 1, quantile, probs = 0.75, na.rm = TRUE)
  p75_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.75, na.rm = TRUE)
  
  # --- dynamic filename ---
  rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
  # e.g. -0.5 -> "m0_5", 0.3 -> "0_3"
  
  setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")
  pdf(paste0("sim_missing_10_50_rho_", rho_str, ".pdf"))
  
  x <- as.numeric(rownames(mape_missing))
  
  plot(x, Median, type="n",
       xlab="Percentage of missings", ylab="MAPE",
       main=paste("MAPE (rho =", rho, ")"),
       ylim = c(0,20))
  
  polygon(c(x, rev(x)), c(p25, rev(p75)),
          col = rgb(0, 0, 1, 0.15), border = NA)
  polygon(c(x, rev(x)), c(p25_sp, rev(p75_sp)),
          col = rgb(1, 0, 0, 0.15), border = NA)
  
  lines(x, Median, type="b", pch=19, col="blue", lwd=2)
  lines(x, Median_sp, type="b", pch=19, col="red", lwd=2)
  
  legend("topleft",
         legend = c("MCFE", "ESF-NNMC"),
         col = c("blue", "red"),
         lwd = 2,
         pch = 19,
         bty = "n")
  
  dev.off()
}


### su matrice generata con AR(phi = 0.5)

rho_vals <- c(0, 0.2, 0.4, 0.6, 0.8)

for (rho in rho_vals) {

  start <- Sys.time()
  
  B = 200
  C = length(seq(2,20,2))
  Wspar = 0.3
  n = 500
  col = 50
  row = 10
  ar = 0.5 
  ns = sqrt(n)
  r = 3
  mape_missing = matrix(0, C, B)
  mape_missing_sp = matrix(0, C, B)
  
  set.seed(180326)
  
  for(i in 1:10) {
    for(j in 1:B) {
      
      W = matrix(rbinom(n,1,Wspar),nrow=row, ncol=row)
      W[lower.tri(W)] = t(W)[lower.tri(W)]
      diag(W) = 0
      
      lr = generate_low_rank_spatiotemp(n = row, p = col, rank = r, W = W, rho = rho, phi = ar, noise_sd = 1e-3, singular_scale = c(50, 15, 5), row_standardize = TRUE, check_stability = TRUE)
      mat = lr$M
      
      mask = matrix(rbinom(n,1,1-((i*2)/100)),row,col)
      
      sel = select_spatial_evec(W)
      A = sel$A
      
      res = mcnnm_cv_R(mat, mask, num_lam = 20,
                       to_estimate_u = TRUE, to_estimate_v = TRUE,
                       num_folds = 5, cv_ratio = 0.6,
                       niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred = res$L + matrix(rep(res$u,col),row,col) + 
             t(matrix(rep(res$v,row),col,row))
      
      res_sp = mcnnm_cv_R_with_A(mat, mask, A, num_lam = 20,
                                to_estimate_alpha = TRUE, to_estimate_v = TRUE,
                                num_folds = 5, cv_ratio = 0.6,
                                niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred_sp <- res_sp$L +
        (A %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L)) +
        matrix(1, nrow(res_sp$L), 1) %*% t(res_sp$v)
      
      mape_missing[i,j] = mean(abs((mat - pred) / mat)[mask == 0]) * 100
      mape_missing_sp[i,j] = mean(abs((mat - pred_sp) / mat)[mask == 0]) * 100
    }
  }
  
  end <- Sys.time()
  print(end - start)
  
  # --- summaries ---
  rownames(mape_missing) <- seq(2,20,2)
  colnames(mape_missing) <- 1:B
  rownames(mape_missing_sp) <- seq(2,20,2)
  colnames(mape_missing_sp) <- 1:B
  
  Median <- apply(mape_missing, 1, median, na.rm=TRUE)
  Median_sp <- apply(mape_missing_sp, 1, median, na.rm=TRUE)
  p25 <- apply(mape_missing, 1, quantile, probs = 0.25, na.rm = TRUE)
  p25_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.25, na.rm = TRUE)
  p75 <- apply(mape_missing, 1, quantile, probs = 0.75, na.rm = TRUE)
  p75_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.75, na.rm = TRUE)
  
  # --- dynamic filename ---
  rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
  # e.g. -0.5 -> "m0_5", 0.3 -> "0_3"
  
  setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")
  pdf(paste0("sim_missing_10_50_ar_rho_", rho_str, ".pdf"))
  
  x <- as.numeric(rownames(mape_missing))
  
  plot(x, Median, type="n",
       xlab="Percentage of missings", ylab="MAPE",
       main=paste("MAPE (rho =", rho, ")"),
       ylim = c(0,20))
  
  polygon(c(x, rev(x)), c(p25, rev(p75)),
          col = rgb(0, 0, 1, 0.15), border = NA)
  polygon(c(x, rev(x)), c(p25_sp, rev(p75_sp)),
          col = rgb(1, 0, 0, 0.15), border = NA)
  
  lines(x, Median, type="b", pch=19, col="blue", lwd=2)
  lines(x, Median_sp, type="b", pch=19, col="red", lwd=2)
  
  legend("topleft",
         legend = c("MCFE", "ESF-NNMC"),
         col = c("blue", "red"),
         lwd = 2,
         pch = 19,
         bty = "n")
  
  dev.off()
}



# Senza time fixed effects

rho_vals <- c(0, 0.2, 0.4, 0.6, 0.8)

for (rho in rho_vals) {

  start <- Sys.time()
  
  B = 200
  C = length(seq(2,20,2))
  Wspar = 0.3
  n = 500
  col = 50
  row = 10
  ns = sqrt(n)
  r = 3
  mape_missing = matrix(0, C, B)
  mape_missing_sp = matrix(0, C, B)
  
  set.seed(180326)
  
  for(i in 1:10) {
    for(j in 1:B) {
      
      W = matrix(rbinom(n,1,Wspar),nrow=row, ncol=row)
      W[lower.tri(W)] = t(W)[lower.tri(W)]
      diag(W) = 0
      
      lr = generate_low_rank_spatial(n = row, p = col, rank = r, W = W, rho = rho, noise_sd = 1e-3, singular_scale = c(50, 15, 5), row_standardize = TRUE, check_stability = TRUE)
      mat = lr$M
      
      mask = matrix(rbinom(n,1,1-((i*2)/100)),row,col)
      
      sel = select_spatial_evec(W)
      A = sel$A
      
      res = mcnnm_cv_R(mat, mask, num_lam = 20,
                       to_estimate_u = TRUE, to_estimate_v = FALSE,
                       num_folds = 5, cv_ratio = 0.6,
                       niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred = res$L + matrix(rep(res$u,col),row,col) 
           # + t(matrix(rep(res$v,row),col,row))
      
      res_sp = mcnnm_cv_R_with_A(mat, mask, A, num_lam = 20,
                                to_estimate_alpha = TRUE, to_estimate_v = FALSE,
                                num_folds = 5, cv_ratio = 0.6,
                                niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred_sp <- res_sp$L +
        (A %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L))
       # + matrix(1, nrow(res_sp$L), 1) %*% t(res_sp$v)
      
      mape_missing[i,j] = mean(abs((mat - pred) / mat)[mask == 0]) * 100
      mape_missing_sp[i,j] = mean(abs((mat - pred_sp) / mat)[mask == 0]) * 100
    }
  }
  
  end <- Sys.time()
  print(end - start)
  
  # --- summaries ---
  rownames(mape_missing) <- seq(2,20,2)
  colnames(mape_missing) <- 1:B
  rownames(mape_missing_sp) <- seq(2,20,2)
  colnames(mape_missing_sp) <- 1:B
  
  Median <- apply(mape_missing, 1, median, na.rm=TRUE)
  Median_sp <- apply(mape_missing_sp, 1, median, na.rm=TRUE)
  p25 <- apply(mape_missing, 1, quantile, probs = 0.25, na.rm = TRUE)
  p25_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.25, na.rm = TRUE)
  p75 <- apply(mape_missing, 1, quantile, probs = 0.75, na.rm = TRUE)
  p75_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.75, na.rm = TRUE)
  
  # --- dynamic filename ---
  rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
  # e.g. -0.5 -> "m0_5", 0.3 -> "0_3"
  
  setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")
  pdf(paste0("sim_missing_notime_10_50_rho_", rho_str, ".pdf"))
  
  x <- as.numeric(rownames(mape_missing))
  
  plot(x, Median, type="n",
       xlab="Percentage of missings", ylab="MAPE",
       main=paste("MAPE (rho =", rho, ")"),
       ylim = c(0,20))
  
  polygon(c(x, rev(x)), c(p25, rev(p75)),
          col = rgb(0, 0, 1, 0.15), border = NA)
  polygon(c(x, rev(x)), c(p25_sp, rev(p75_sp)),
          col = rgb(1, 0, 0, 0.15), border = NA)
  
  lines(x, Median, type="b", pch=19, col="blue", lwd=2)
  lines(x, Median_sp, type="b", pch=19, col="red", lwd=2)
  
  legend("topleft",
         legend = c("MCFE", "ESF-NNMC"),
         col = c("blue", "red"),
         lwd = 2,
         pch = 19,
         bty = "n")
  
  dev.off()
}

### su matrice generata con AR(phi = 0.5)

rho_vals <- c(0.8)

for (rho in rho_vals) {

  start <- Sys.time()
  
  B = 200
  C = length(seq(2,20,2))
  Wspar = 0.3
  n = 500
  col = 50
  row = 10
  ar = 0.5 
  ns = sqrt(n)
  r = 3
  mape_missing = matrix(0, C, B)
  mape_missing_sp = matrix(0, C, B)
  
  set.seed(180326)
  
  for(i in 1:10) {
    for(j in 1:B) {
      
      W = matrix(rbinom(n,1,Wspar),nrow=row, ncol=row)
      W[lower.tri(W)] = t(W)[lower.tri(W)]
      diag(W) = 0
      
      lr = generate_low_rank_spatiotemp(n = row, p = col, rank = r, W = W, rho = rho, phi = ar, noise_sd = 1e-3, singular_scale = c(50, 15, 5), row_standardize = TRUE, check_stability = TRUE)
      mat = lr$M
      
      mask = matrix(rbinom(n,1,1-((i*2)/100)),row,col)
      
      sel = select_spatial_evec(W)
      A = sel$A
      
      res = mcnnm_cv_R(mat, mask, num_lam = 20,
                       to_estimate_u = TRUE, to_estimate_v = FALSE,
                       num_folds = 5, cv_ratio = 0.6,
                       niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred = res$L + matrix(rep(res$u,col),row,col) 
      
      res_sp = mcnnm_cv_R_with_A(mat, mask, A, num_lam = 20,
                                to_estimate_alpha = TRUE, to_estimate_v = FALSE,
                                num_folds = 5, cv_ratio = 0.6,
                                niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred_sp <- res_sp$L +
        (A %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L))
      
      mape_missing[i,j] = mean(abs((mat - pred) / mat)[mask == 0]) * 100
      mape_missing_sp[i,j] = mean(abs((mat - pred_sp) / mat)[mask == 0]) * 100
    }
  }
  
  end <- Sys.time()
  print(end - start)
  
  # --- summaries ---
  rownames(mape_missing) <- seq(2,20,2)
  colnames(mape_missing) <- 1:B
  rownames(mape_missing_sp) <- seq(2,20,2)
  colnames(mape_missing_sp) <- 1:B
  
  Median <- apply(mape_missing, 1, median, na.rm=TRUE)
  Median_sp <- apply(mape_missing_sp, 1, median, na.rm=TRUE)
  p25 <- apply(mape_missing, 1, quantile, probs = 0.25, na.rm = TRUE)
  p25_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.25, na.rm = TRUE)
  p75 <- apply(mape_missing, 1, quantile, probs = 0.75, na.rm = TRUE)
  p75_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.75, na.rm = TRUE)
  
  # --- dynamic filename ---
  rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
  # e.g. -0.5 -> "m0_5", 0.3 -> "0_3"
  
  setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")
  pdf(paste0("sim_missing_notime_10_50_ar_rho_", rho_str, ".pdf"))
  
  x <- as.numeric(rownames(mape_missing))
  
  plot(x, Median, type="n",
       xlab="Percentage of missings", ylab="MAPE",
       main=paste("MAPE (rho =", rho, ")"),
       ylim = c(0,20))
  
  polygon(c(x, rev(x)), c(p25, rev(p75)),
          col = rgb(0, 0, 1, 0.15), border = NA)
  polygon(c(x, rev(x)), c(p25_sp, rev(p75_sp)),
          col = rgb(1, 0, 0, 0.15), border = NA)
  
  lines(x, Median, type="b", pch=19, col="blue", lwd=2)
  lines(x, Median_sp, type="b", pch=19, col="red", lwd=2)
  
  legend("topleft",
         legend = c("MCFE", "ESF-NNMC"),
         col = c("blue", "red"),
         lwd = 2,
         pch = 19,
         bty = "n")
  
  dev.off()
}




######### 26 03 2026: 50 righe e 10 colonne 

rho_vals <- c(0.4, 0.6, 0.8)

for (rho in rho_vals) {

  start <- Sys.time()
  B = 200
  C = length(seq(2,20,2))
  Wspar = 0.3
  n = 500
  col = 10
  row = 50 
  ns = sqrt(n)
  r = 3
  mape_missing = matrix(0, C, B)
  mape_missing_sp = matrix(0, C, B)
  
  set.seed(180326)
  
  for(i in 1:10) {
    for(j in 1:B) {
      
      W = matrix(rbinom(n,1,Wspar),nrow=row, ncol=row)
      W[lower.tri(W)] = t(W)[lower.tri(W)]
      diag(W) = 0
      
      lr = generate_low_rank_spatial(n = row, p = col, rank = r, W = W, rho = rho, noise_sd = 1e-3, singular_scale = c(50, 15, 5), row_standardize = TRUE, check_stability = TRUE)
      mat = lr$M
      
      mask = matrix(rbinom(n,1,1-((i*2)/100)),row,col)
      
      sel = select_spatial_evec(W)
      A = sel$A
      
      res = mcnnm_cv_R(mat, mask, num_lam = 20,
                       to_estimate_u = TRUE, to_estimate_v = TRUE,
                       num_folds = 5, cv_ratio = 0.6,
                       niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred = res$L + matrix(rep(res$u,col),row,col) + 
             t(matrix(rep(res$v,row),col,row))
      
      res_sp = mcnnm_cv_R_with_A(mat, mask, A, num_lam = 20,
                                to_estimate_alpha = TRUE, to_estimate_v = TRUE,
                                num_folds = 5, cv_ratio = 0.6,
                                niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred_sp <- res_sp$L +
        (A %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L)) +
        matrix(1, nrow(res_sp$L), 1) %*% t(res_sp$v)
      
      mape_missing[i,j] = mean(abs((mat - pred) / mat)[mask == 0]) * 100
      mape_missing_sp[i,j] = mean(abs((mat - pred_sp) / mat)[mask == 0]) * 100
    }
  }
  
  end <- Sys.time()
  print(end - start)
  
  # --- summaries ---
  rownames(mape_missing) <- seq(2,20,2)
  colnames(mape_missing) <- 1:B
  rownames(mape_missing_sp) <- seq(2,20,2)
  colnames(mape_missing_sp) <- 1:B
  
  Median <- apply(mape_missing, 1, median, na.rm=TRUE)
  Median_sp <- apply(mape_missing_sp, 1, median, na.rm=TRUE)
  p25 <- apply(mape_missing, 1, quantile, probs = 0.25, na.rm = TRUE)
  p25_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.25, na.rm = TRUE)
  p75 <- apply(mape_missing, 1, quantile, probs = 0.75, na.rm = TRUE)
  p75_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.75, na.rm = TRUE)
  
  # --- dynamic filename ---
  rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
  # e.g. -0.5 -> "m0_5", 0.3 -> "0_3"
  
  setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")
  pdf(paste0("sim_missing_50_10_rho_", rho_str, ".pdf"))
  
  x <- as.numeric(rownames(mape_missing))
  
  plot(x, Median, type="n",
       xlab="Percentage of missings", ylab="MAPE",
       main=paste("MAPE (rho =", rho, ")"),
       ylim = c(0,20))
  
  polygon(c(x, rev(x)), c(p25, rev(p75)),
          col = rgb(0, 0, 1, 0.15), border = NA)
  polygon(c(x, rev(x)), c(p25_sp, rev(p75_sp)),
          col = rgb(1, 0, 0, 0.15), border = NA)
  
  lines(x, Median, type="b", pch=19, col="blue", lwd=2)
  lines(x, Median_sp, type="b", pch=19, col="red", lwd=2)
  
  legend("topleft",
         legend = c("MCFE", "ESF-NNMC"),
         col = c("blue", "red"),
         lwd = 2,
         pch = 19,
         bty = "n")
  
  dev.off()
}


### su matrice generata con AR(phi = 0.5)

rho_vals <- c(0, 0.2, 0.4, 0.6, 0.8)

for (rho in rho_vals) {

  start <- Sys.time()
  
  B = 200
  C = length(seq(2,20,2))
  Wspar = 0.3
  n = 500
  col = 10
  row = 50
  ar = 0.5 
  ns = sqrt(n)
  r = 3
  mape_missing = matrix(0, C, B)
  mape_missing_sp = matrix(0, C, B)
  
  set.seed(180326)
  
  for(i in 1:10) {
    for(j in 1:B) {
      
      W = matrix(rbinom(n,1,Wspar),nrow=row, ncol=row)
      W[lower.tri(W)] = t(W)[lower.tri(W)]
      diag(W) = 0
      
      lr = generate_low_rank_spatiotemp(n = row, p = col, rank = r, W = W, rho = rho, phi = ar, noise_sd = 1e-3, singular_scale = c(50, 15, 5), row_standardize = TRUE, check_stability = TRUE)
      mat = lr$M
      
      mask = matrix(rbinom(n,1,1-((i*2)/100)),row,col)
      
      sel = select_spatial_evec(W)
      A = sel$A
      
      res = mcnnm_cv_R(mat, mask, num_lam = 20,
                       to_estimate_u = TRUE, to_estimate_v = TRUE,
                       num_folds = 5, cv_ratio = 0.6,
                       niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred = res$L + matrix(rep(res$u,col),row,col) + 
             t(matrix(rep(res$v,row),col,row))
      
      res_sp = mcnnm_cv_R_with_A(mat, mask, A, num_lam = 20,
                                to_estimate_alpha = TRUE, to_estimate_v = TRUE,
                                num_folds = 5, cv_ratio = 0.6,
                                niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred_sp <- res_sp$L +
        (A %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L)) +
        matrix(1, nrow(res_sp$L), 1) %*% t(res_sp$v)
      
      mape_missing[i,j] = mean(abs((mat - pred) / mat)[mask == 0]) * 100
      mape_missing_sp[i,j] = mean(abs((mat - pred_sp) / mat)[mask == 0]) * 100
    }
  }
  
  end <- Sys.time()
  print(end - start)
  
  # --- summaries ---
  rownames(mape_missing) <- seq(2,20,2)
  colnames(mape_missing) <- 1:B
  rownames(mape_missing_sp) <- seq(2,20,2)
  colnames(mape_missing_sp) <- 1:B
  
  Median <- apply(mape_missing, 1, median, na.rm=TRUE)
  Median_sp <- apply(mape_missing_sp, 1, median, na.rm=TRUE)
  p25 <- apply(mape_missing, 1, quantile, probs = 0.25, na.rm = TRUE)
  p25_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.25, na.rm = TRUE)
  p75 <- apply(mape_missing, 1, quantile, probs = 0.75, na.rm = TRUE)
  p75_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.75, na.rm = TRUE)
  
  # --- dynamic filename ---
  rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
  # e.g. -0.5 -> "m0_5", 0.3 -> "0_3"
  
  setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")
  pdf(paste0("sim_missing_50_10_ar_rho_", rho_str, ".pdf"))
  
  x <- as.numeric(rownames(mape_missing))
  
  plot(x, Median, type="n",
       xlab="Percentage of missings", ylab="MAPE",
       main=paste("MAPE (rho =", rho, ")"),
       ylim = c(0,20))
  
  polygon(c(x, rev(x)), c(p25, rev(p75)),
          col = rgb(0, 0, 1, 0.15), border = NA)
  polygon(c(x, rev(x)), c(p25_sp, rev(p75_sp)),
          col = rgb(1, 0, 0, 0.15), border = NA)
  
  lines(x, Median, type="b", pch=19, col="blue", lwd=2)
  lines(x, Median_sp, type="b", pch=19, col="red", lwd=2)
  
  legend("topleft",
         legend = c("MCFE", "ESF-NNMC"),
         col = c("blue", "red"),
         lwd = 2,
         pch = 19,
         bty = "n")
  
  dev.off()
}



# Senza time fixed effects

rho_vals <- c(0, 0.2, 0.4, 0.6, 0.8)

for (rho in rho_vals) {

  start <- Sys.time()
  
  B = 200
  C = length(seq(2,20,2))
  Wspar = 0.3
  n = 500
  col = 10
  row = 50
  ns = sqrt(n)
  r = 3
  mape_missing = matrix(0, C, B)
  mape_missing_sp = matrix(0, C, B)
  
  set.seed(180326)
  
  for(i in 1:10) {
    for(j in 1:B) {
      
      W = matrix(rbinom(n,1,Wspar),nrow=row, ncol=row)
      W[lower.tri(W)] = t(W)[lower.tri(W)]
      diag(W) = 0
      
      lr = generate_low_rank_spatial(n = row, p = col, rank = r, W = W, rho = rho, noise_sd = 1e-3, singular_scale = c(50, 15, 5), row_standardize = TRUE, check_stability = TRUE)
      mat = lr$M
      
      mask = matrix(rbinom(n,1,1-((i*2)/100)),row,col)
      
      sel = select_spatial_evec(W)
      A = sel$A
      
      res = mcnnm_cv_R(mat, mask, num_lam = 20,
                       to_estimate_u = TRUE, to_estimate_v = FALSE,
                       num_folds = 5, cv_ratio = 0.6,
                       niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred = res$L + matrix(rep(res$u,col),row,col) 
           # + t(matrix(rep(res$v,row),col,row))
      
      res_sp = mcnnm_cv_R_with_A(mat, mask, A, num_lam = 20,
                                to_estimate_alpha = TRUE, to_estimate_v = FALSE,
                                num_folds = 5, cv_ratio = 0.6,
                                niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred_sp <- res_sp$L +
        (A %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L))
       # + matrix(1, nrow(res_sp$L), 1) %*% t(res_sp$v)
      
      mape_missing[i,j] = mean(abs((mat - pred) / mat)[mask == 0]) * 100
      mape_missing_sp[i,j] = mean(abs((mat - pred_sp) / mat)[mask == 0]) * 100
    }
  }
  
  end <- Sys.time()
  print(end - start)
  
  # --- summaries ---
  rownames(mape_missing) <- seq(2,20,2)
  colnames(mape_missing) <- 1:B
  rownames(mape_missing_sp) <- seq(2,20,2)
  colnames(mape_missing_sp) <- 1:B
  
  Median <- apply(mape_missing, 1, median, na.rm=TRUE)
  Median_sp <- apply(mape_missing_sp, 1, median, na.rm=TRUE)
  p25 <- apply(mape_missing, 1, quantile, probs = 0.25, na.rm = TRUE)
  p25_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.25, na.rm = TRUE)
  p75 <- apply(mape_missing, 1, quantile, probs = 0.75, na.rm = TRUE)
  p75_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.75, na.rm = TRUE)
  
  # --- dynamic filename ---
  rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
  # e.g. -0.5 -> "m0_5", 0.3 -> "0_3"
  
  setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")
  pdf(paste0("sim_missing_notime_50_10_rho_", rho_str, ".pdf"))
  
  x <- as.numeric(rownames(mape_missing))
  
  plot(x, Median, type="n",
       xlab="Percentage of missings", ylab="MAPE",
       main=paste("MAPE (rho =", rho, ")"),
       ylim = c(0,20))
  
  polygon(c(x, rev(x)), c(p25, rev(p75)),
          col = rgb(0, 0, 1, 0.15), border = NA)
  polygon(c(x, rev(x)), c(p25_sp, rev(p75_sp)),
          col = rgb(1, 0, 0, 0.15), border = NA)
  
  lines(x, Median, type="b", pch=19, col="blue", lwd=2)
  lines(x, Median_sp, type="b", pch=19, col="red", lwd=2)
  
  legend("topleft",
         legend = c("MCFE", "ESF-NNMC"),
         col = c("blue", "red"),
         lwd = 2,
         pch = 19,
         bty = "n")
  
  dev.off()
}

### su matrice generata con AR(phi = 0.5)

rho_vals <- c(0, 0.2, 0.4. 0.6, 0.8)

for (rho in rho_vals) {

  start <- Sys.time()
  
  B = 200
  C = length(seq(2,20,2))
  Wspar = 0.3
  n = 500
  col = 10
  row = 50
  ar = 0.5 
  ns = sqrt(n)
  r = 3
  mape_missing = matrix(0, C, B)
  mape_missing_sp = matrix(0, C, B)
  
  set.seed(180326)
  
  for(i in 1:10) {
    for(j in 1:B) {
      
      W = matrix(rbinom(n,1,Wspar),nrow=row, ncol=row)
      W[lower.tri(W)] = t(W)[lower.tri(W)]
      diag(W) = 0
      
      lr = generate_low_rank_spatiotemp(n = row, p = col, rank = r, W = W, rho = rho, phi = ar, noise_sd = 1e-3, singular_scale = c(50, 15, 5), row_standardize = TRUE, check_stability = TRUE)
      mat = lr$M
      
      mask = matrix(rbinom(n,1,1-((i*2)/100)),row,col)
      
      sel = select_spatial_evec(W)
      A = sel$A
      
      res = mcnnm_cv_R(mat, mask, num_lam = 20,
                       to_estimate_u = TRUE, to_estimate_v = FALSE,
                       num_folds = 5, cv_ratio = 0.6,
                       niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred = res$L + matrix(rep(res$u,col),row,col) 
      
      res_sp = mcnnm_cv_R_with_A(mat, mask, A, num_lam = 20,
                                to_estimate_alpha = TRUE, to_estimate_v = FALSE,
                                num_folds = 5, cv_ratio = 0.6,
                                niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred_sp <- res_sp$L +
        (A %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L))
      
      mape_missing[i,j] = mean(abs((mat - pred) / mat)[mask == 0]) * 100
      mape_missing_sp[i,j] = mean(abs((mat - pred_sp) / mat)[mask == 0]) * 100
    }
  }
  
  end <- Sys.time()
  print(end - start)
  
  # --- summaries ---
  rownames(mape_missing) <- seq(2,20,2)
  colnames(mape_missing) <- 1:B
  rownames(mape_missing_sp) <- seq(2,20,2)
  colnames(mape_missing_sp) <- 1:B
  
  Median <- apply(mape_missing, 1, median, na.rm=TRUE)
  Median_sp <- apply(mape_missing_sp, 1, median, na.rm=TRUE)
  p25 <- apply(mape_missing, 1, quantile, probs = 0.25, na.rm = TRUE)
  p25_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.25, na.rm = TRUE)
  p75 <- apply(mape_missing, 1, quantile, probs = 0.75, na.rm = TRUE)
  p75_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.75, na.rm = TRUE)
  
  # --- dynamic filename ---
  rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
  # e.g. -0.5 -> "m0_5", 0.3 -> "0_3"
  
  setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")
  pdf(paste0("sim_missing_notime_50_10_ar_rho_", rho_str, ".pdf"))
  
  x <- as.numeric(rownames(mape_missing))
  
  plot(x, Median, type="n",
       xlab="Percentage of missings", ylab="MAPE",
       main=paste("MAPE (rho =", rho, ")"),
       ylim = c(0,20))
  
  polygon(c(x, rev(x)), c(p25, rev(p75)),
          col = rgb(0, 0, 1, 0.15), border = NA)
  polygon(c(x, rev(x)), c(p25_sp, rev(p75_sp)),
          col = rgb(1, 0, 0, 0.15), border = NA)
  
  lines(x, Median, type="b", pch=19, col="blue", lwd=2)
  lines(x, Median_sp, type="b", pch=19, col="red", lwd=2)
  
  legend("topleft",
         legend = c("MCFE", "ESF-NNMC"),
         col = c("blue", "red"),
         lwd = 2,
         pch = 19,
         bty = "n")
  
  dev.off()
}


######### 28 03 2026: 50 righe e 50 colonne 

#rho_vals <- c(0, 0.2, 0.4, 0.6, 0.8)
rho_vals <- c(0.4, 0.6, 0.8)
for (rho in rho_vals) {

  start <- Sys.time()
  B = 200
  C = length(seq(2,20,2))
  Wspar = 0.3
  n = 500
  col = 50
  row = 50 
  ns = sqrt(n)
  r = 5
  mape_missing = matrix(0, C, B)
  mape_missing_sp = matrix(0, C, B)
  
  set.seed(180326)
  
  for(i in 1:10) {
    for(j in 1:B) {
      
      W = matrix(rbinom(n,1,Wspar),nrow=row, ncol=row)
      W[lower.tri(W)] = t(W)[lower.tri(W)]
      diag(W) = 0
      
      lr = generate_low_rank_spatial(n = row, p = col, rank = r, W = W, rho = rho, noise_sd = 1e-3, singular_scale = c(40, 16, 4, 3, 2), row_standardize = TRUE, check_stability = TRUE)
      mat = lr$M
      
      mask = matrix(rbinom(n,1,1-((i*2)/100)),row,col)
      
      sel = select_spatial_evec(W)
      A = sel$A
      
      res = mcnnm_cv_R(mat, mask, num_lam = 20,
                       to_estimate_u = TRUE, to_estimate_v = TRUE,
                       num_folds = 5, cv_ratio = 0.6,
                       niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred = res$L + matrix(rep(res$u,col),row,col) + 
             t(matrix(rep(res$v,row),col,row))
      
      res_sp = mcnnm_cv_R_with_A(mat, mask, A, num_lam = 20,
                                to_estimate_alpha = TRUE, to_estimate_v = TRUE,
                                num_folds = 5, cv_ratio = 0.6,
                                niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred_sp <- res_sp$L +
        (A %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L)) +
        matrix(1, nrow(res_sp$L), 1) %*% t(res_sp$v)
      
      mape_missing[i,j] = mean(abs((mat - pred) / mat)[mask == 0]) * 100
      mape_missing_sp[i,j] = mean(abs((mat - pred_sp) / mat)[mask == 0]) * 100
    }
  }
  
  end <- Sys.time()
  print(end - start)
  
  # --- summaries ---
  rownames(mape_missing) <- seq(2,20,2)
  colnames(mape_missing) <- 1:B
  rownames(mape_missing_sp) <- seq(2,20,2)
  colnames(mape_missing_sp) <- 1:B
  
  Median <- apply(mape_missing, 1, median, na.rm=TRUE)
  Median_sp <- apply(mape_missing_sp, 1, median, na.rm=TRUE)
  p25 <- apply(mape_missing, 1, quantile, probs = 0.25, na.rm = TRUE)
  p25_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.25, na.rm = TRUE)
  p75 <- apply(mape_missing, 1, quantile, probs = 0.75, na.rm = TRUE)
  p75_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.75, na.rm = TRUE)
  
  # --- dynamic filename ---
  rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
  # e.g. -0.5 -> "m0_5", 0.3 -> "0_3"
  
  setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")
  pdf(paste0("sim_missing_50_50_rho_", rho_str, ".pdf"))
  
  x <- as.numeric(rownames(mape_missing))
  
  plot(x, Median, type="n",
       xlab="Percentage of missings", ylab="MAPE",
       main=paste("MAPE (rho =", rho, ")"),
       ylim = c(0,20))
  
  polygon(c(x, rev(x)), c(p25, rev(p75)),
          col = rgb(0, 0, 1, 0.15), border = NA)
  polygon(c(x, rev(x)), c(p25_sp, rev(p75_sp)),
          col = rgb(1, 0, 0, 0.15), border = NA)
  
  lines(x, Median, type="b", pch=19, col="blue", lwd=2)
  lines(x, Median_sp, type="b", pch=19, col="red", lwd=2)
  
  legend("topleft",
         legend = c("MCFE", "ESF-NNMC"),
         col = c("blue", "red"),
         lwd = 2,
         pch = 19,
         bty = "n")
  
  dev.off()
}


### su matrice generata con AR(phi = 0.5)

rho_vals <- c(0, 0.2, 0.4, 0.6, 0.8)

for (rho in rho_vals) {

  start <- Sys.time()
  
  B = 200
  C = length(seq(2,20,2))
  Wspar = 0.3
  n = 500
  col = 50
  row = 50
  ar = 0.5 
  ns = sqrt(n)
  r = 5
  mape_missing = matrix(0, C, B)
  mape_missing_sp = matrix(0, C, B)
  
  set.seed(180326)
  
  for(i in 1:10) {
    for(j in 1:B) {
      
      W = matrix(rbinom(n,1,Wspar),nrow=row, ncol=row)
      W[lower.tri(W)] = t(W)[lower.tri(W)]
      diag(W) = 0
      
      lr = generate_low_rank_spatiotemp(n = row, p = col, rank = r, W = W, rho = rho, phi = ar, noise_sd = 1e-3, singular_scale = c(40, 16, 4, 3, 2), row_standardize = TRUE, check_stability = TRUE)
      mat = lr$M
      
      mask = matrix(rbinom(n,1,1-((i*2)/100)),row,col)
      
      sel = select_spatial_evec(W)
      A = sel$A
      
      res = mcnnm_cv_R(mat, mask, num_lam = 20,
                       to_estimate_u = TRUE, to_estimate_v = TRUE,
                       num_folds = 5, cv_ratio = 0.6,
                       niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred = res$L + matrix(rep(res$u,col),row,col) + 
             t(matrix(rep(res$v,row),col,row))
      
      res_sp = mcnnm_cv_R_with_A(mat, mask, A, num_lam = 20,
                                to_estimate_alpha = TRUE, to_estimate_v = TRUE,
                                num_folds = 5, cv_ratio = 0.6,
                                niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred_sp <- res_sp$L +
        (A %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L)) +
        matrix(1, nrow(res_sp$L), 1) %*% t(res_sp$v)
      
      mape_missing[i,j] = mean(abs((mat - pred) / mat)[mask == 0]) * 100
      mape_missing_sp[i,j] = mean(abs((mat - pred_sp) / mat)[mask == 0]) * 100
    }
  }
  
  end <- Sys.time()
  print(end - start)
  
  # --- summaries ---
  rownames(mape_missing) <- seq(2,20,2)
  colnames(mape_missing) <- 1:B
  rownames(mape_missing_sp) <- seq(2,20,2)
  colnames(mape_missing_sp) <- 1:B
  
  Median <- apply(mape_missing, 1, median, na.rm=TRUE)
  Median_sp <- apply(mape_missing_sp, 1, median, na.rm=TRUE)
  p25 <- apply(mape_missing, 1, quantile, probs = 0.25, na.rm = TRUE)
  p25_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.25, na.rm = TRUE)
  p75 <- apply(mape_missing, 1, quantile, probs = 0.75, na.rm = TRUE)
  p75_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.75, na.rm = TRUE)
  
  # --- dynamic filename ---
  rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
  # e.g. -0.5 -> "m0_5", 0.3 -> "0_3"
  
  setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")
  pdf(paste0("sim_missing_50_50_ar_rho_", rho_str, ".pdf"))
  
  x <- as.numeric(rownames(mape_missing))
  
  plot(x, Median, type="n",
       xlab="Percentage of missings", ylab="MAPE",
       main=paste("MAPE (rho =", rho, ")"),
       ylim = c(0,20))
  
  polygon(c(x, rev(x)), c(p25, rev(p75)),
          col = rgb(0, 0, 1, 0.15), border = NA)
  polygon(c(x, rev(x)), c(p25_sp, rev(p75_sp)),
          col = rgb(1, 0, 0, 0.15), border = NA)
  
  lines(x, Median, type="b", pch=19, col="blue", lwd=2)
  lines(x, Median_sp, type="b", pch=19, col="red", lwd=2)
  
  legend("topleft",
         legend = c("MCFE", "ESF-NNMC"),
         col = c("blue", "red"),
         lwd = 2,
         pch = 19,
         bty = "n")
  
  dev.off()
}



# Senza time fixed effects

rho_vals <- c(0, 0.2, 0.4, 0.6, 0.8)

for (rho in rho_vals) {

  start <- Sys.time()
  
  B = 200
  C = length(seq(2,20,2))
  Wspar = 0.3
  n = 500
  col = 50
  row = 50
  ns = sqrt(n)
  r = 5
  mape_missing = matrix(0, C, B)
  mape_missing_sp = matrix(0, C, B)
  
  set.seed(180326)
  
  for(i in 1:10) {
    for(j in 1:B) {
      
      W = matrix(rbinom(n,1,Wspar),nrow=row, ncol=row)
      W[lower.tri(W)] = t(W)[lower.tri(W)]
      diag(W) = 0
      
      lr = generate_low_rank_spatial(n = row, p = col, rank = r, W = W, rho = rho, noise_sd = 1e-3, singular_scale = c(40, 16, 4, 3, 2), row_standardize = TRUE, check_stability = TRUE)
      mat = lr$M
      
      mask = matrix(rbinom(n,1,1-((i*2)/100)),row,col)
      
      sel = select_spatial_evec(W)
      A = sel$A
      
      res = mcnnm_cv_R(mat, mask, num_lam = 20,
                       to_estimate_u = TRUE, to_estimate_v = FALSE,
                       num_folds = 5, cv_ratio = 0.6,
                       niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred = res$L + matrix(rep(res$u,col),row,col) 
           # + t(matrix(rep(res$v,row),col,row))
      
      res_sp = mcnnm_cv_R_with_A(mat, mask, A, num_lam = 20,
                                to_estimate_alpha = TRUE, to_estimate_v = FALSE,
                                num_folds = 5, cv_ratio = 0.6,
                                niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred_sp <- res_sp$L +
        (A %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L))
       # + matrix(1, nrow(res_sp$L), 1) %*% t(res_sp$v)
      
      mape_missing[i,j] = mean(abs((mat - pred) / mat)[mask == 0]) * 100
      mape_missing_sp[i,j] = mean(abs((mat - pred_sp) / mat)[mask == 0]) * 100
    }
  }
  
  end <- Sys.time()
  print(end - start)
  
  # --- summaries ---
  rownames(mape_missing) <- seq(2,20,2)
  colnames(mape_missing) <- 1:B
  rownames(mape_missing_sp) <- seq(2,20,2)
  colnames(mape_missing_sp) <- 1:B
  
  Median <- apply(mape_missing, 1, median, na.rm=TRUE)
  Median_sp <- apply(mape_missing_sp, 1, median, na.rm=TRUE)
  p25 <- apply(mape_missing, 1, quantile, probs = 0.25, na.rm = TRUE)
  p25_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.25, na.rm = TRUE)
  p75 <- apply(mape_missing, 1, quantile, probs = 0.75, na.rm = TRUE)
  p75_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.75, na.rm = TRUE)
  
  # --- dynamic filename ---
  rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
  # e.g. -0.5 -> "m0_5", 0.3 -> "0_3"
  
  setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")
  pdf(paste0("sim_missing_notime_50_50_rho_", rho_str, ".pdf"))
  
  x <- as.numeric(rownames(mape_missing))
  
  plot(x, Median, type="n",
       xlab="Percentage of missings", ylab="MAPE",
       main=paste("MAPE (rho =", rho, ")"),
       ylim = c(0,20))
  
  polygon(c(x, rev(x)), c(p25, rev(p75)),
          col = rgb(0, 0, 1, 0.15), border = NA)
  polygon(c(x, rev(x)), c(p25_sp, rev(p75_sp)),
          col = rgb(1, 0, 0, 0.15), border = NA)
  
  lines(x, Median, type="b", pch=19, col="blue", lwd=2)
  lines(x, Median_sp, type="b", pch=19, col="red", lwd=2)
  
  legend("topleft",
         legend = c("MCFE", "ESF-NNMC"),
         col = c("blue", "red"),
         lwd = 2,
         pch = 19,
         bty = "n")
  
  dev.off()
}

### su matrice generata con AR(phi = 0.5)

rho_vals <- c(0, 0.2, 0.4. 0.6, 0.8)

for (rho in rho_vals) {

  start <- Sys.time()
  
  B = 200
  C = length(seq(2,20,2))
  Wspar = 0.3
  col = 50
  row = 50
  ar = 0.5 
  ns = sqrt(n)
  r = 5
  mape_missing = matrix(0, C, B)
  mape_missing_sp = matrix(0, C, B)
  
  set.seed(180326)
  
  for(i in 1:10) {
    for(j in 1:B) {
      
      W = matrix(rbinom(n,1,Wspar),nrow=row, ncol=row)
      W[lower.tri(W)] = t(W)[lower.tri(W)]
      diag(W) = 0
      
      lr = generate_low_rank_spatiotemp(n = row, p = col, rank = r, W = W, rho = rho, phi = ar, noise_sd = 1e-3, singular_scale = c(40, 16, 4, 3, 2), row_standardize = TRUE, check_stability = TRUE)
      mat = lr$M
      
      mask = matrix(rbinom(n,1,1-((i*2)/100)),row,col)
      
      sel = select_spatial_evec(W)
      A = sel$A
      
      res = mcnnm_cv_R(mat, mask, num_lam = 20,
                       to_estimate_u = TRUE, to_estimate_v = FALSE,
                       num_folds = 5, cv_ratio = 0.6,
                       niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred = res$L + matrix(rep(res$u,col),row,col) 
      
      res_sp = mcnnm_cv_R_with_A(mat, mask, A, num_lam = 20,
                                to_estimate_alpha = TRUE, to_estimate_v = FALSE,
                                num_folds = 5, cv_ratio = 0.6,
                                niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
      
      pred_sp <- res_sp$L +
        (A %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L))
      
      mape_missing[i,j] = mean(abs((mat - pred) / mat)[mask == 0]) * 100
      mape_missing_sp[i,j] = mean(abs((mat - pred_sp) / mat)[mask == 0]) * 100
    }
  }
  
  end <- Sys.time()
  print(end - start)
  
  # --- summaries ---
  rownames(mape_missing) <- seq(2,20,2)
  colnames(mape_missing) <- 1:B
  rownames(mape_missing_sp) <- seq(2,20,2)
  colnames(mape_missing_sp) <- 1:B
  
  Median <- apply(mape_missing, 1, median, na.rm=TRUE)
  Median_sp <- apply(mape_missing_sp, 1, median, na.rm=TRUE)
  p25 <- apply(mape_missing, 1, quantile, probs = 0.25, na.rm = TRUE)
  p25_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.25, na.rm = TRUE)
  p75 <- apply(mape_missing, 1, quantile, probs = 0.75, na.rm = TRUE)
  p75_sp <- apply(mape_missing_sp, 1, quantile, probs = 0.75, na.rm = TRUE)
  
  # --- dynamic filename ---
  rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
  # e.g. -0.5 -> "m0_5", 0.3 -> "0_3"
  
  setwd("G:\\Il mio Drive\\00. PRIN_PNRR_2022\\02. Analisi\\Spatial matrix completion\\Dataset")
  pdf(paste0("sim_missing_notime_50_50_ar_rho_", rho_str, ".pdf"))
  
  x <- as.numeric(rownames(mape_missing))
  
  plot(x, Median, type="n",
       xlab="Percentage of missings", ylab="MAPE",
       main=paste("MAPE (rho =", rho, ")"),
       ylim = c(0,20))
  
  polygon(c(x, rev(x)), c(p25, rev(p75)),
          col = rgb(0, 0, 1, 0.15), border = NA)
  polygon(c(x, rev(x)), c(p25_sp, rev(p75_sp)),
          col = rgb(1, 0, 0, 0.15), border = NA)
  
  lines(x, Median, type="b", pch=19, col="blue", lwd=2)
  lines(x, Median_sp, type="b", pch=19, col="red", lwd=2)
  
  legend("topleft",
         legend = c("MCFE", "ESF-NNMC"),
         col = c("blue", "red"),
         lwd = 2,
         pch = 19,
         bty = "n")
  
  dev.off()
}




# stessa cosa ma con:
50 x 50




### al variare della sparsità di W (missing: 10%, quadrata, dim = 10 x 10, rank = 3, rho = 0.5)

start <- Sys.time()
rho = 0.5
B = 100
C= length(seq(0.1,0.55,0.05))
i = 10
n = 100 
col=10
row=10
ns= sqrt(n)
mape_Wspar = matrix(0,C,B)
mape_Wspar_sp = matrix(0,C,B)

for(Wspar in 1:C) {
for(j in 1:B) {
lr = generate_low_rank_spat(n = 10, rank = 3, singular_scale = c(50, 15, 5), noise_sd = 1e-4, rho=0.5, W = W, row_standardize = TRUE, check_stability = TRUE)
mat = lr$M
mask = matrix(rbinom(n,1,1-((i-1)/100)),row,col)
W = matrix(rbinom(n,1,Wspar/20 + 0.05),nrow=row, ncol=col)
sel = select_spatial_evec(W)
A = sel$A
res = mcnnm_cv_R(mat, mask, num_lam = 20,to_estimate_u = TRUE, to_estimate_v = TRUE, num_folds = 5, cv_ratio = 0.6, niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
pred = res$L + matrix(rep(res$u,col),row,col) + t(matrix(rep(res$v,row),col,row))
res_sp = mcnnm_cv_R_with_A(mat, mask, A, num_lam = 20,to_estimate_alpha = TRUE, to_estimate_v = TRUE, num_folds = 5, cv_ratio = 0.6, niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
pred_sp <- res_sp$L + (A %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L)) + matrix(1, nrow(res_sp$L), 1) %*% t(res_sp$v)
mape_Wspar[Wspar,j] = mean(abs((mat - pred) / mat)[mask == 0]) * 100
mape_Wspar_sp[Wspar,j] = mean(abs((mat - pred_sp) / mat)[mask == 0]) * 100
}
}
end <- Sys.time()
end - start

rownames(mape_Wspar) <- seq(0.1,0.55,0.05)
colnames(mape_Wspar) <- 1:B
rownames(mape_Wspar_sp) <- seq(0.1,0.55,0.05)
colnames(mape_Wspar_sp) <- 1:B

apply(mape_Wspar, 1, summary)
apply(mape_Wspar_sp, 1, summary)

Median <- apply(mape_Wspar, 1, median, na.rm=TRUE)
Median_sp <- apply(mape_Wspar_sp, 1, median, na.rm=TRUE)
p25 <- apply(mape_Wspar, 1, quantile, probs = 0.25, na.rm = TRUE)
p25_sp <- apply(mape_Wspar_sp, 1, quantile, probs = 0.25, na.rm = TRUE)
p75 <- apply(mape_Wspar, 1, quantile, probs = 0.75, na.rm = TRUE)
p75_sp <- apply(mape_Wspar_sp, 1, quantile, probs = 0.75, na.rm = TRUE)


# Plot
x <- as.numeric(rownames(mape_Wspar))
plot(x, Median, type="n",
     xlab="Percentage of valid connections in W", ylab="MAPE",
     main="MAPE at different levels of sparseness in W",
     ylim = c(0,50))

# Bands
polygon(c(x, rev(x)), c(p25, rev(p75)),
        col = rgb(0, 0, 1, 0.15), border = NA)
polygon(c(x, rev(x)), c(p25_sp, rev(p75_sp)),
        col = rgb(1, 0, 0, 0.15), border = NA)

# Median lines on top
lines(x, Median, type="b", pch=19, col="blue", lwd=2)
lines(x, Median_sp, type="b", pch=19, col="red", lwd=2)

legend("topleft",
       legend = c("MCFE", "ESF-NNMC"),
       col = c("blue", "red"),
       lwd = 2,
       pch = 19,
       bty = "n")

###########
#### al variare di rho ######

start <- Sys.time()
rho = 0.5
B = 
C= length(seq(0,0.9,0.1))
Wspar= 0.4
m = 10
n = 100 
col=10
row=10
ns= sqrt(n)
mape_rho = matrix(0,C,B)
mape_rho_sp = matrix(0,C,B)
set.seed(413)

for(i in 1:C) {
for(j in 1:B) {
lr = generate_low_rank_spat(n = 10, rank = 3, singular_scale = c(50, 15, 5), noise_sd = 1e-4, rho=(i-1)*0.1, W = W, row_standardize = TRUE, check_stability = TRUE)
mat = lr$M
mask = matrix(rbinom(n,1,1-((m-1)/100)),row,col)
W = matrix(rbinom(n,1,Wspar),nrow=row, ncol=col)
sel = select_spatial_evec(W)
A = sel$A
res = mcnnm_cv_R(mat, mask, num_lam = 20,to_estimate_u = TRUE, to_estimate_v = TRUE, num_folds = 5, cv_ratio = 0.6, niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
pred = res$L + matrix(rep(res$u,col),row,col) + t(matrix(rep(res$v,row),col,row))
res_sp = mcnnm_cv_R_with_A(mat, mask, A, num_lam = 20,to_estimate_alpha = TRUE, to_estimate_v = TRUE, num_folds = 5, cv_ratio = 0.6, niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
pred_sp <- res_sp$L + (A %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L)) + matrix(1, nrow(res_sp$L), 1) %*% t(res_sp$v)
mape_rho[i,j] = mean(abs((mat - pred) / mat)[mask == 0]) * 100
mape_rho_sp[i,j] = mean(abs((mat - pred_sp) / mat)[mask == 0]) * 100
}
}
end <- Sys.time()
end - start

rownames(mape_rho) <- seq(0.1,1,0.1)
colnames(mape_rho) <- 1:B
rownames(mape_rho_sp) <- seq(0.1,1,0.1)
colnames(mape_rho_sp) <- 1:B

apply(mape_rho, 1, summary)
apply(mape_rho_sp, 1, summary)

Median <- apply(mape_rho, 1, median, na.rm=TRUE)
Median_sp <- apply(mape_rho_sp, 1, median, na.rm=TRUE)
p25 <- apply(mape_rho, 1, quantile, probs = 0.25, na.rm = TRUE)
p25_sp <- apply(mape_rho_sp, 1, quantile, probs = 0.25, na.rm = TRUE)
p75 <- apply(mape_rho, 1, quantile, probs = 0.75, na.rm = TRUE)
p75_sp <- apply(mape_rho_sp, 1, quantile, probs = 0.75, na.rm = TRUE)


# Plot
x <- as.numeric(rownames(mape_rho))
plot(x, Median, type="n",
     xlab="Rho", ylab="MAPE",
     main="MAPE at different levels of parameter rho",
     ylim = c(0,120))

# Bands
polygon(c(x, rev(x)), c(p25, rev(p75)),
        col = rgb(0, 0, 1, 0.15), border = NA)
polygon(c(x, rev(x)), c(p25_sp, rev(p75_sp)),
        col = rgb(1, 0, 0, 0.15), border = NA)

# Median lines on top
lines(x, Median, type="b", pch=19, col="blue", lwd=2)
lines(x, Median_sp, type="b", pch=19, col="red", lwd=2)

legend("topleft",
       legend = c("MCFE", "ESF-NNMC"),
       col = c("blue", "red"),
       lwd = 2,
       pch = 19,
       bty = "n")

#############





### al variare della dimensione (missing: 10%, Wspar = 0.4, quadrata, rank = 5)

start <- Sys.time()
B = 100
C= length(c(10^2,20^2,30^2,40^2,50^2,60^2,70^2, 80^2,90^2, 100^2))
i = 10
Wspar = 0.4
mape_dim = matrix(0,C,B)
mape_dim_sp = matrix(0,C,B)

for(dim in 1:C) {
for(j in 1:B) {
n = (dim*10)^2
ns = sqrt(n)
lr = generate_low_rank(n = ns, rank = 5, singular_scale = c(40, 12, 8, 6, 4), noise_sd = 1e-4)
mat = lr$M
mask = matrix(rbinom(n,1,1-((i-1)/100)),ns,ns)
W = matrix(rbinom(n,1,Wspar),nrow=ns, ncol=ns)
sel = select_spatial_evec(W)
A = sel$A
res = mcnnm_cv_R(mat, mask, num_lam = 20,to_estimate_u = TRUE, to_estimate_v = TRUE, num_folds = 5, cv_ratio = 0.6, niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
pred = res$L + matrix(rep(res$u,ns),ns,ns) + t(matrix(rep(res$v,ns),ns,ns))
res_sp = mcnnm_cv_R_with_A(mat, mask, A, num_lam = 20,to_estimate_alpha = TRUE, to_estimate_v = TRUE, num_folds = 5, cv_ratio = 0.6, niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
pred_sp <- res_sp$L + (A %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L)) + matrix(1, nrow(res_sp$L), 1) %*% t(res_sp$v)
mape_dim[dim,j] = mean(abs((mat - pred) / mat)[mask == 0]) * 100
mape_dim_sp[dim,j] = mean(abs((mat - pred_sp) / mat)[mask == 0]) * 100
}
}
end <- Sys.time()
end - start

rownames(mape_missing) <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
colnames(mape_missing) <- 1:B
rownames(mape_missing_sp) <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
colnames(mape_missing_sp) <- 1:B

#par(mfrow=c(1,2))
#hist(mape_missing[)
#hist(mape_missing_sp)
apply(mape_missing, 1, summary)
apply(mape_missing_sp, 1, summary)

stats <- apply(mape_missing, 1, summary)
stats_sp <- apply(mape_missing_sp, 1, summary)

Median <- stats["Median", ]
Median_sp <- stats_sp["Median", ]

# Plot medians first
plot(Median, type="b", pch=19, col="blue",
     xaxt="n", xlab="Dimension of the square matrix (number of columns)", ylab="MAPE",
     main="MAPE at different dimensions of the square matrix", ylim = c(0,300))

axis(1, at = 1:length(Median), labels = colnames(stats))

# Add Q1 and Q3
lines(Median_sp, type="b", pch=19, col="red")

legend("topleft",
       legend = c("Median - MCFE", "Median - SPMCFE"),
       col = c("blue", "red"),
       pch = 19)



###al variare della forma della matrice

start <- Sys.time()
B = 100
Wspar= 0.4
n = 100
col = 25
row = 4
i=10
ns = sqrt(n)
mape_shape5 = matrix(0,1,B)
mape_shape5_sp = matrix(0,1,B)

for(j in 1:B) {
lr = generate_low_rank_nonsquared(m = row, n= col, rank = 3, singular_scale = c(50, 15, 5), noise_sd = 1e-4)
mat = lr$M
mask = matrix(rbinom(n,1,1-((i)/100)),row,col)
W = matrix(rbinom(n,1,Wspar),nrow=row, ncol=row)
sel = select_spatial_evec(W)
A = sel$A
res = mcnnm_cv_R(mat, mask, num_lam = 20,to_estimate_u = TRUE, to_estimate_v = TRUE, num_folds = 5, cv_ratio = 0.6, niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
pred = res$L + matrix(rep(res$u,col),row,col) + t(matrix(rep(res$v,row),col,row))
res_sp = mcnnm_cv_R_with_A(mat, mask, A, num_lam = 20,to_estimate_alpha = TRUE, to_estimate_v = TRUE, num_folds = 5, cv_ratio = 0.6, niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
pred_sp <- res_sp$L + (A %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L)) + matrix(1, nrow(res_sp$L), 1) %*% t(res_sp$v)
mape_shape5[,j] = mean(abs((mat - pred) / mat)[mask == 0]) * 100
mape_shape5_sp[,j] = mean(abs((mat - pred_sp) / mat)[mask == 0]) * 100
}

end <- Sys.time()
end - start

mape_shape = rbind(mape_shape1, mape_shape2, mape_shape3, mape_shape4, mape_shape5)
mape_shape_sp = rbind(mape_shape1_sp, mape_shape2_sp, mape_shape3_sp, mape_shape4_sp, mape_shape5_sp)

rownames(mape_shape) <- c(4,5,10,20,25)
colnames(mape_shape) <- 1:B
rownames(mape_shape_sp) <- c(4,5,10,20,25)
colnames(mape_shape_sp) <- 1:B

apply(mape_shape, 1, summary)
apply(mape_shape_sp, 1, summary)

stats <- apply(mape_shape[c(2:4),], 1, summary)
stats_sp <- apply(mape_shape_sp[c(2:4),], 1, summary)

Median <- stats["Median", ]
Median_sp <- stats_sp["Median", ]

# Plot medians first
plot(Median, type="b", pch=19, col="blue",
     xaxt="n", xlab="Number of columns", ylab="MAPE",
     main="MAPE at different number of columns (matrix dimension = 100)", ylim = c(0,80))

axis(1, at = 1:length(Median), labels = colnames(stats))

# Add Q1 and Q3
lines(Median_sp, type="b", pch=19, col="red")

legend("bottomleft",
       legend = c("Median - MCFE", "Median - SPMCFE"),
       col = c("blue", "red"),
       pch = 19)


