# Set your directory
setwd("G:/Il mio Drive/00. PRIN_PNRR_2022/02. Analisi/Spatial matrix completion/GITHUB")

# source files
source("R/FENNMC.R")
source("R/ESFNNMC.R")
source("R/Weight_matrix.R")
source("R/Spatial_eigenvectors.R")


set.seed(123)


# --------------------------------------------------
# 1. Synthetic coordinates and spatial weights
# --------------------------------------------------
n <- 40
T <- 40
coords <- cbind(
  Longitude = runif(n, 8, 11),
  Latitude  = runif(n, 44, 47)
)

# the number of neighbours must be chosen
k = 5
W_knn <- build_knn_weights(coords = coords, k, symmetrize = TRUE)

# --------------------------------------------------
# 2. Spatial eigenvectors used by ESFNNMC
# --------------------------------------------------
tau = 0.90
sel_knn <- select_spatial_evec(W_knn, explained = tau)
A_knn <- sel_knn$A

# --------------------------------------------------
# 3. Generate a synthetic matrix
# --------------------------------------------------

M_obs <- matrix(rnorm(n * T, mean = 10, sd = 2), nrow = n, ncol = T)

# --------------------------------------------------
# 4. Introduce missing entries
# --------------------------------------------------

missing_rate <- 0.10

mask <- matrix(
  rbinom(n * T, 1, 1 - missing_rate),
  n,
  T
)

# The algorithm uses a numeric matrix and a
# separate observation mask.
M <- M_obs
M[mask == 0] <- 0

cat(
  "Percentage of missing entries:",
  round(100 * mean(mask == 0), 1),
  "%\n"
)

# --------------------------------------------------
# 5. Fit ESFNNMC and FENNMC 
# --------------------------------------------------

fit_FENNMC <- fennmc(
  M = M,
  mask = mask,
      num_lam = 20,
      to_estimate_u = TRUE,
      to_estimate_v = TRUE,
      num_folds = 5,
      cv_ratio = 0.6,
      niter = 200,
      rel_tol = 1e-5,
      is_quiet = TRUE
    )
    
fit_ESFNNMC <- esfnnmc(
  M = M,
  mask = mask,
  A = A_knn,
  num_lam = 20,
  to_estimate_alpha = TRUE,
  to_estimate_v = TRUE,
  num_folds = 5,
  cv_ratio = 0.6,
  niter = 200,
  rel_tol = 1e-5,
  is_quiet = TRUE
)


# --------------------------------------------------
# 6. Reconstruct the matrix
# --------------------------------------------------

M_hat_FENNMC <- compute_matrix(
  fit_FENNMC$L,
  fit_FENNMC$u,
  fit_FENNMC$v
)

M_hat_ESFNNMC <- compute_matrix_with_A(
  fit_ESFNNMC$L,
  A_knn,
  fit_ESFNNMC$alpha,
  fit_ESFNNMC$v
)



# --------------------------------------------------
# 7. Evaluate reconstruction on missing entries
# --------------------------------------------------

test <- mask == 0

# ESFNNMC
rmse <- sqrt( mean( (M_obs[test] - M_hat_ESFNNMC[test])^2 ) )

cat(
  "RMSE on artificially missing entries:",
  round(rmse, 4),
  "\n"
)

mape <- mean(
  abs((M_obs[test] - M_hat_ESFNNMC[test]) / M_obs[test])
) * 100

cat(
  "MAPE on artificially missing entries:",
  round(mape, 4),
  "\n"
)

cat(
  "Selected lambda:",
  fit_ESFNNMC$best_lambda,
  "\n"
)

# FENNMC
rmse <- sqrt( mean( (M_obs[test] - M_hat_FENNMC[test])^2 ) )

cat(
  "RMSE on artificially missing entries:",
  round(rmse, 4),
  "\n"
)

mape <- mean(
  abs((M_obs[test] - M_hat_FENNMC[test]) / M_obs[test])
) * 100

cat(
  "MAPE on artificially missing entries:",
  round(mape, 4),
  "\n"
)

cat(
  "Selected lambda:",
  fit_FENNMC$best_lambda,
  "\n"
)