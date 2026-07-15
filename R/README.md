# R source files

Files to include in the `R/` folder:

- `00_matrix_completion_helpers.R`: shared matrix-completion helper functions.
- `01_fennmc_benchmark.R`: fixed-effects nuclear-norm matrix completion benchmark.
- `02_esfnnmc.R`: proposed ESFNNMC estimator.
- `03_spatial_eigenvectors.R`: Moran eigenvector selection utilities.
- `04_spatial_weights.R`: utilities to build spatial weight matrices.

Recommended loading order:

```r
source("R/00_matrix_completion_helpers.R")
source("R/01_fennmc_benchmark.R")
source("R/02_esfnnmc.R")
source("R/03_spatial_eigenvectors.R")
source("R/04_spatial_weights.R")
```

Minimal ESFNNMC workflow with synthetic data:

```r
set.seed(123)

# 1. Synthetic coordinates and spatial weights
n <- 40
p <- 80
coords <- cbind(
  Longitude = runif(n, 8, 11),
  Latitude  = runif(n, 44, 47)
)

W_knn <- build_knn_weights(coords = coords, k = 5, symmetrize = TRUE)

# 2. Spatial eigenvectors used by ESFNNMC
sel_knn <- select_spatial_evec(W_knn, explained = 0.90)
A_knn <- sel_knn$A

# 3. Synthetic low-rank data matrix with spatial and time components
rank_true <- 3
U <- matrix(rnorm(n * rank_true), n, rank_true)
V <- matrix(rnorm(p * rank_true), p, rank_true)
L_true <- U %*% t(V)

alpha_true <- rnorm(ncol(A_knn))
v_true <- rnorm(p)

M_true <- L_true +
  as.vector(A_knn %*% alpha_true) %*% matrix(1, 1, p) +
  matrix(1, n, 1) %*% t(v_true)

M_obs <- M_true + matrix(rnorm(n * p, sd = 0.5), n, p)

# 4. Artificial missingness
mask <- matrix(rbinom(n * p, size = 1, prob = 0.80), n, p)
M_input <- M_obs
M_input[mask == 0] <- 0

# 5. ESFNNMC fit
fit <- mcnnm_cv_R_with_A(
  M = M_input,
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

# 6. Reconstructed matrix
M_hat <- fit$L +
  as.vector(A_knn %*% fit$alpha) %*% matrix(1, 1, p) +
  matrix(1, n, 1) %*% t(fit$v)
```
