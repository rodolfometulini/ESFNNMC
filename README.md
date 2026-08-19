# Eigenvector Spatial Filters Nuclear Norm Matrix Completion

This repository contains the R code accompanying the paper:

**Eigenvector Spatial Filters Nuclear Norm Matrix Completion with Application to Air Quality Data** (https://arxiv.org/abs/2606.05450)

The repository provides implementations of:

- **ESFNNMC** — Eigenvector Spatial Filters Nuclear Norm Matrix Completion, the spatial extension proposed in the paper.
- **FENNMC** — Fixed Effects Nuclear Norm Matrix Completion (Athey et al., 2021), used as the main benchmark;

ESFNNMC extends fixed-effects nuclear norm matrix completion by replacing unrestricted unit fixed effects with a parsimonious representation based on Moran eigenvectors derived from a spatial weights matrix.

The repository also contains reproducible examples illustrating how the methods can be used to reconstruct missing entries in partially observed panel and spatio-temporal matrices.

## Model overview

Let $M \in \mathbb{R}^{n \times T}$ be a partially observed matrix.

In short, the FENNMC estimator (Athey et al., 2021) decomposes the matrix as $M = L + u\mathbf{1}_T^\top + \mathbf{1}_n v^\top + E$,

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

**ESFNNMC**

The proposed ESFNNMC method modifies the FENNMC specification by replacing the unrestricted unit effects $u$ with $u = A\alpha$, where the columns of \(A\) are selected Moran eigenvectors derived from a spatial weights matrix.

The resulting decomposition is $M = L + A\alpha\mathbf{1}_T^\top + \mathbf{1}_n v^\top + E$.

This provides a parsimonious and interpretable representation of spatially structured unit heterogeneity.


## Repository structure

```text
.
├── README.md
├── R/
│   ├── FENNMC.R
│   └── ESFNNMC.R
├── examples/
│   ├── toy_FENNMC.R
│   └── toy_ESFNNMC.R
└── application/
|__ simulation/
    └── ...
```

The repository will be progressively updated with the code required to reproduce the simulation study and the empirical application presented in the paper.

# ESFNNMC and FENNMC in R

The implementation of Eigenvector Spatial Filters Nuclear Norm Matrix Completion (ESFNNMC) and Fixed Effects Nuclear Norm Matrix Completion (FENNMC) is contained in:

```text
R/FENNMC.R
R/ESFNNMC.R
R/Weight_matrix.R
R/Spatial_eigenvectors.R
```

Load the functions using (in this order):

```r
source("R/FENNMC.R")
source("R/FENNMC.R")
source("R/Weight_matrix.R")
source("R/Spatial_eigenvectors.R")
```

The main functions are

```r
fennmc()
esfnnmc()
```

which estimate the model and selects the nuclear-norm regularization parameter by cross-validation.

The basic usage is:

```r
fit <- fennmc(
  M,
  mask,
  to_estimate_u = TRUE,
  to_estimate_v = TRUE,
  num_lam = 20,
  num_folds = 5,
  cv_ratio = 0.6,
  niter = 200,
  rel_tol = 1e-5,
  is_quiet = FALSE
)


fit <- esfnnmc(
  M,
  mask,
  A,
  to_estimate_alpha = TRUE,
  to_estimate_v = TRUE,
  num_lam = 20,
  num_folds = 5,
  cv_ratio = 0.6,
  niter = 200,
  rel_tol = 1e-5,
  is_quiet = FALSE
)


```

where:

- `M` is the numeric matrix used by the algorithm;
- `mask` is a binary matrix with the same dimensions as `M`;
- `mask[i, j] = 1` indicates an observed entry;
- `mask[i, j] = 0` indicates an entry to be reconstructed.
- `A` is the numeric matrix containing in columns the q selected eigenvectors

Missing entries in `M` should be replaced by zero before fitting. Missingness is specified through `mask`, rather than through `NA` values.

The fitted matrix for FENNMC can be reconstructed as

```r
M_hat <- compute_matrix(
  fit$L,
  fit$u,
  fit$v
)
```
while for ESFNNMC can be reconstructed as

```r
M_hat = compute_matrix_with_A(
  fit$L,
  A,
  fit$alpha,
  fit$v
)

```

The returned object also contains the selected regularization parameter and cross-validation results:

```r
fit$best_lambda
fit$min_RMSE
fit$Avg_RMSE
```

# Toy example

The following example generates a simple synthetic matrix, randomly removes 10% of its entries, and reconstructs them using ESFNNMC and FENNMC.

```r
source("R/FENNMC.R")
source("R/FENNMC.R")
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

M_obs <- matrix(rnorm(n * T), nrow = n, ncol = T)

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
  fit$L,
  fit$u,
  fit$v
)

M_hat_ESFNNMC <- compute_matrix_with_A(
  fit$L,
  A,
  fit$alpha,
  fit$v
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
```

# Disabling fixed effects

The implementation also allows unit or time fixed effects to be excluded.

For example, FE nuclear norm matrix completion without time effects can be estimated with

```r
fit_FENNMC_notime <- fennmc(
  M = M,
  mask = mask,
      num_lam = 20,
      to_estimate_u = TRUE,
      to_estimate_v = FALSE,
      num_folds = 5,
      cv_ratio = 0.6,
      niter = 200,
      rel_tol = 1e-5,
      is_quiet = TRUE
    )
```
and ESF nuclear norm matrix completion without unit effects can be estimated with

```r
fit_ESFNNMC_notime <- esfnnmc(
M = M,
  mask = mask,
  A = A_knn,
  num_lam = 20,
  to_estimate_alpha = TRUE,
  to_estimate_v = FALSE,
  num_folds = 5,
  cv_ratio = 0.6,
  niter = 200,
  rel_tol = 1e-5,
  is_quiet = TRUE
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

ARXiv CC BY 4.0

# References

```text
Athey, S., Bayati, M., Doudchenko, N., Imbens, G., & Khosravi, K. (2021). 
Matrix completion methods for causal panel data models. 
Journal of the American Statistical Association, 116(536), 1716-1730.
```