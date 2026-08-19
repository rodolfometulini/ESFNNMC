# Eigenvector Spatial Filters Nuclear Norm Matrix Completion

This repository contains the R code accompanying the paper:

**Eigenvector Spatial Filters Nuclear Norm Matrix Completion with Application to Air Quality Data** (https://arxiv.org/abs/2606.05450)

The repository provides implementations of:

- **FENNMC** — Fixed Effects Nuclear Norm Matrix Completion, used as the main benchmark;
- **ESFNNMC** — Eigenvector Spatial Filters Nuclear Norm Matrix Completion, the spatial extension proposed in the paper.

ESFNNMC extends fixed-effects nuclear norm matrix completion by replacing unrestricted unit fixed effects with a parsimonious representation based on Moran eigenvectors derived from a spatial weights matrix.

The repository also contains reproducible examples illustrating how the methods can be used to reconstruct missing entries in partially observed panel and spatio-temporal matrices.

## Model overview

Let $M \in \mathbb{R}^{n \times T}$ be a partially observed matrix.

The FENNMC estimator (Athey et al., 2021) decomposes the matrix as $M = L + u\mathbf{1}_T^\top + \mathbf{1}_n v^\top + E$,

where:

- $L$ is a low-rank latent component;
- $u$ contains unit fixed effects;
- $v$ contains time fixed effects;
- $E$ is the residual component.

The estimator solves $\min_{L,u,v} \frac{1}{2} \sum_{(i,t)\in\Omega} (M_{it}-L_{it}-u_i-v_t)^2 +\lambda \|L\|_*$,

where $\|L\|_*$ denotes the nuclear norm.

Estimation is performed through a block-coordinate descent algorithm alternating between:

1. unit fixed effects;
2. time fixed effects;
3. the low-rank component, updated through singular-value soft thresholding.

The regularization parameter $\lambda$ is selected by cross-validation.

## Repository structure

```text
.
├── README.md
├── R/
│   ├── MCFE_soft_impute.R
│   └── ESFNNMC.R
├── examples/
│   ├── toy_FENNMC.R
│   └── toy_ESFNNMC.R
└── application/
    └── ...
```

The repository will be progressively updated with the code required to reproduce the simulation study and the empirical application presented in the paper.

# FENNMC in R

The implementation of Fixed Effects Nuclear Norm Matrix Completion is contained in:

```text
R/FENNMC.R
```

Load the functions using:

```r
source("R/FENNMC.R")
```

The main function is

```r
mcnnm_cv_R()
```

which estimates the model and selects the nuclear-norm regularization parameter by cross-validation.

Its basic usage is:

```r
fit <- mcnnm_cv_R(
  M,
  mask,
  to_estimate_u = TRUE,
  to_estimate_v = TRUE,
  num_lam = 20,
  num_folds = 5
)
```

where:

- `M` is the numeric matrix used by the algorithm;
- `mask` is a binary matrix with the same dimensions as `M`;
- `mask[i, j] = 1` indicates an observed entry;
- `mask[i, j] = 0` indicates an entry to be reconstructed.

Missing entries in `M` should be replaced by zero before fitting. Missingness is specified through `mask`, rather than through `NA` values.

The fitted matrix can be reconstructed as

```r
M_hat <- compute_matrix(
  fit$L,
  fit$u,
  fit$v
)
```

The returned object also contains the selected regularization parameter and cross-validation results:

```r
fit$best_lambda
fit$min_RMSE
fit$Avg_RMSE
```

# Toy example: Fixed Effects Nuclear Norm Matrix Completion

The following example generates a simple synthetic matrix, randomly removes 10% of its entries, and reconstructs them using FENNMC.

```r
source("R/FENNMC.R")

set.seed(123)

# --------------------------------------------------
# 1. Generate a synthetic matrix
# --------------------------------------------------
n <- 40
T <- 40
M_complete <- matrix(rnorm(n * T), nrow = n, ncol = T)

# --------------------------------------------------
# 2. Introduce missing entries
# --------------------------------------------------

missing_rate <- 0.20

mask <- matrix(
  rbinom(n * T, 1, 1 - missing_rate),
  n,
  T
)

# The algorithm uses a numeric matrix and a
# separate observation mask.
M <- M_complete
M[mask == 0] <- 0

cat(
  "Percentage of missing entries:",
  round(100 * mean(mask == 0), 1),
  "%\n"
)

# --------------------------------------------------
# 3. Fit FENNMC
# --------------------------------------------------

fit <- mcnnm_cv_R(
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

# --------------------------------------------------
# 4. Reconstruct the matrix
# --------------------------------------------------

M_hat <- compute_matrix(
  fit$L,
  fit$u,
  fit$v
)

# --------------------------------------------------
# 5. Evaluate reconstruction on missing entries
# --------------------------------------------------

test <- mask == 0

rmse_missing <- sqrt( mean( (M_complete[test] - M_hat[test])^2 ) )

cat(
  "RMSE on artificially missing entries:",
  round(rmse_missing, 4),
  "\n"
)

mape_missing <- mean(
  abs((M_complete[test] - M_hat[test]) / M_complete[test])
) * 100

cat(
  "MAPE on artificially missing entries:",
  round(mape_missing, 4),
  "\n"
)

cat(
  "Selected lambda:",
  fit$best_lambda,
  "\n"
)
```

The example illustrates the basic workflow:

```text
complete matrix
      ↓
artificial missingness
      ↓
observation mask
      ↓
cross-validation for lambda
      ↓
FENNMC estimation
      ↓
matrix reconstruction
      ↓
out-of-sample evaluation
```

# Disabling fixed effects

The implementation also allows unit or time fixed effects to be excluded.

For example, nuclear norm matrix completion without unit effects can be estimated with

```r
fit <- mcnnm_cv_R(
  M,
  mask,
  to_estimate_u = FALSE,
  to_estimate_v = TRUE
)
```

and a model without either set of fixed effects with

```r
fit <- mcnnm_cv_R(
  M,
  mask,
  to_estimate_u = FALSE,
  to_estimate_v = FALSE
)
```

# ESFNNMC

The proposed ESFNNMC method modifies the FENNMC specification by replacing the unrestricted unit effects $u$ with $u = A\alpha$, where the columns of \(A\) are selected Moran eigenvectors derived from a spatial weights matrix.

The resulting decomposition is $M = L + A\alpha\mathbf{1}_T^\top + \mathbf{1}_n v^\top + E$.

This provides a parsimonious and interpretable representation of spatially structured unit heterogeneity.

# ESFNNMC in R

The implementation of Eigenvector Spatial Filters Nuclear Norm Matrix Completion is contained in:

```text
R/ESFNNMC.R
```

Load the functions using (In the following order):

```r
source("R/FENNMC.R")
source("R/Weight_matrix.R")
source("R/Spatial_eigenvectors.R")

```

The main function is

```r
mcnnm_cv_R_with_A()
```

Its basic usage is:

```r
fit <- mcnnm_cv_R_with_A(
  M = M,
  mask = mask,
  A = A,
  num_lam = 20,
  to_estimate_alpha = TRUE,
  to_estimate_v = TRUE,
  num_folds = 5,
  cv_ratio = 0.6,
  niter = 200,
  rel_tol = 1e-5,
  is_quiet = TRUE
)

```

where:

- `M` is the numeric matrix used by the algorithm;
- `mask` is a binary matrix with the same dimensions as `M`;
- `mask[i, j] = 1` indicates an observed entry;
- `mask[i, j] = 0` indicates an entry to be reconstructed.
- `A` is the numeric matrix containing in columns the q selected eigenvectors

Missing entries in `M` should be replaced by zero before fitting. Missingness is specified through `mask`, rather than through `NA` values.

The fitted matrix can be reconstructed as

```r
M_hat <- fit$L +
  as.vector(A_knn %*% fit$alpha) %*% matrix(1, 1, p) +
  matrix(1, n, 1) %*% t(fit$v)
```

The returned object also contains the selected regularization parameter and cross-validation results:

```r
fit$best_lambda
fit$min_RMSE
fit$Avg_RMSE
```


# Toy example: ESF Nuclear Norm Matrix Completion

The following example generates a simple synthetic matrix, randomly removes 10% of its entries, and reconstructs them using ESFNNMC.

```r
set.seed(123)

# 1. Synthetic coordinates and spatial weights
n <- 40
T <- 40
coords <- cbind(
  Longitude = runif(n, 8, 11),
  Latitude  = runif(n, 44, 47)
)

# the number of neighbours must be chosen
W_knn <- build_knn_weights(coords = coords, k = 5, symmetrize = TRUE)

# 2. Spatial eigenvectors used by ESFNNMC
sel_knn <- select_spatial_evec(W_knn, explained = 0.90)
A_knn <- sel_knn$A

# 3. Synthetic matrix

M_obs <- matrix(rnorm(n * T), nrow = n, ncol = T)

# 4. Artificial missingness
mask <- matrix(rbinom(n * T, size = 1, prob = 0.90), n, T)
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
  as.vector(A_knn %*% fit$alpha) %*% matrix(1, 1, T) +
  matrix(1, n, 1) %*% t(fit$v)
  
# --------------------------------------------------
# 7. Evaluate reconstruction on missing entries
# --------------------------------------------------

test <- mask == 0

rmse_missing <- sqrt( mean( (M_input[test] - M_hat[test])^2 ) )

cat(
  "RMSE on artificially missing entries:",
  round(rmse_missing, 4),
  "\n"
)

mape_missing <- mean(
  abs((M_obs[test] - M_hat[test]) / M_obs[test])
) * 100

cat(
  "MAPE on artificially missing entries:",
  round(mape_missing, 4),
  "\n"
)

cat(
  "Selected lambda:",
  fit$best_lambda,
  "\n"
)
  
```

# Citation

If you use this code, please cite:

```text
Metulini, R. (2026). 
Eigenvector Spatial Filters Nuclear Norm Matrix Completion with Application to Air Quality Data. 
arXiv preprint arXiv:2606.05450.
```

Citation information will be updated after publication.

# License

TBD

# References

```text
Athey, S., Bayati, M., Doudchenko, N., Imbens, G., & Khosravi, K. (2021). 
Matrix completion methods for causal panel data models. 
Journal of the American Statistical Association, 116(536), 1716-1730.
```