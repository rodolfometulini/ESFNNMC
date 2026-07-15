# R source files

- `00_matrix_completion_helpers.R`: shared helper functions.
- `01_fennmc_benchmark.R`: fixed-effects nuclear-norm matrix completion benchmark.
- `02_esfnnmc.R`: proposed ESFNNMC estimator.
- `03_spatial_eigenvectors.R`: Moran eigenvector selection utilities.

Recommended loading order:

```r
source("R/00_matrix_completion_helpers.R")
source("R/01_fennmc_benchmark.R")
source("R/03_spatial_eigenvectors.R")
source("R/02_esfnnmc.R")
```
