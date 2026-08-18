# Eigenvector Spatial Filters Nuclear Norm Matrix Completion

This repository contains the R code accompanying the paper:

**Eigenvector Spatial Filters Nuclear Norm Matrix Completion with Application to Air Quality Data**

The repository provides implementations of:

- **FENNMC** — Fixed Effects Nuclear Norm Matrix Completion, used as the main benchmark;
- **ESFNNMC** — Eigenvector Spatial Filters Nuclear Norm Matrix Completion, the spatial extension proposed in the paper.

ESFNNMC extends fixed-effects nuclear norm matrix completion by replacing unrestricted unit fixed effects with a parsimonious representation based on Moran eigenvectors derived from a spatial weights matrix.

The repository also contains reproducible examples illustrating how the methods can be used to reconstruct missing entries in partially observed panel and spatio-temporal matrices.

## Model overview

Let

$M \in \mathbb{R}^{n \times T}$

be a partially observed matrix.

The FENNMC estimator decomposes the matrix as

$M = L + u\mathbf{1}_T^\top + \mathbf{1}_n v^\top + E$,


where:

- \(L\) is a low-rank latent component;
- \(u\) contains unit fixed effects;
- \(v\) contains time fixed effects;
- \(E\) is the residual component.

The estimator solves

$\min_{L,u,v} \frac{1}{2} \sum_{(i,t)\in\Omega} (M_{it}-L_{it}-u_i-v_t)^2 +\lambda \|L\|_*$,

where \(\|L\|_*\) denotes the nuclear norm.

Estimation is performed through a block-coordinate descent algorithm alternating between:

1. unit fixed effects;
2. time fixed effects;
3. the low-rank component, updated through singular-value soft thresholding.

The regularization parameter \(\lambda\) can be selected by cross-validation.

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
R/MCFE_soft_impute.R
```

Load the functions using:

```r
source("R/MCFE_soft_impute.R")
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

The following example generates a low-rank panel matrix with unit and time effects, randomly removes 20% of its entries, and reconstructs them using FENNMC.

```r
source("R/MCFE_soft_impute.R")

set.seed(123)

# --------------------------------------------------
# 1. Generate a synthetic low-rank panel
# --------------------------------------------------

n <- 20
T <- 30
r <- 3

U <- matrix(rnorm(n * r), n, r)
V <- matrix(rnorm(T * r), T, r)

L_true <- U %*% t(V)

# Unit fixed effects
u_true <- rnorm(n, sd = 2)

# Time fixed effects
v_true <- sin(seq(0, 2 * pi, length.out = T)) * 2

# Complete matrix
M_complete <-
  L_true +
  outer(u_true, rep(1, T)) +
  outer(rep(1, n), v_true) +
  matrix(rnorm(n * T, sd = 0.2), n, T)

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
  to_estimate_u = TRUE,
  to_estimate_v = TRUE,
  num_lam = 20,
  num_folds = 5,
  cv_ratio = 0.8,
  is_quiet = FALSE
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

rmse_missing <- sqrt(
  mean(
    (M_complete[test] - M_hat[test])^2
  )
)

cat(
  "RMSE on artificially missing entries:",
  round(rmse_missing, 4),
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

The proposed ESFNNMC method modifies the FENNMC specification by replacing the unrestricted unit effects

$u$

with

$u = A\alpha$,

where the columns of \(A\) are selected Moran eigenvectors derived from a spatial weights matrix.

The resulting decomposition is

$ M = L + A\alpha\mathbf{1}_T^\top + \mathbf{1}_n v^\top + E$.

This provides a parsimonious and interpretable representation of spatially structured unit heterogeneity.

Code and examples for ESFNNMC are provided in the corresponding sections of this repository.

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