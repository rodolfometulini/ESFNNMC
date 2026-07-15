# R source files

- `00_matrix_completion_helpers.R`: shared matrix-completion helper functions.
- `01_fennmc_benchmark.R`: fixed-effects nuclear-norm matrix completion benchmark.
- `02_esfnnmc.R`: proposed ESFNNMC estimator.
- `03_spatial_eigenvectors.R`: Moran eigenvector selection utilities.
- `04_spatial_weights.R`: utilities to build spatial weight matrices from coordinates.

Recommended loading order:

```r
source("R/00_matrix_completion_helpers.R")
source("R/01_fennmc_benchmark.R")
source("R/02_esfnnmc.R")
source("R/03_spatial_eigenvectors.R")
source("R/04_spatial_weights.R")
```

Minimal ESFNNMC workflow:

```r
coords <- stations[, c("Longitude", "Latitude")]

spatial_filters <- build_knn_spatial_filters(
  coords = coords,
  k = 10,
  station_ids = stations$IDStations,
  explained = 0.90
)

W_knn <- spatial_filters$W
A_knn <- spatial_filters$A

fit <- mcnnm_cv_R_with_A(
  M = mat_r,
  mask = mask,
  A = A_knn,
  num_lam = 20,
  to_estimate_alpha = TRUE,
  to_estimate_v = TRUE,
  num_folds = 5,
  cv_ratio = 0.6,
  niter = 200,
  rel_tol = 1e-5
)
```
