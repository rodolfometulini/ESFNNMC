## functions for the generation of synthetic data

## The processes are described in section 3.1

# low-rank spatio-temporal matrices
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


# low-rank matrices with spatial dependence and unit fixed effects
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
