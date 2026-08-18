# R source files

Files to include in the `R/` folder:

- `FENNMC.R`: fixed-effects nuclear-norm matrix completion benchmark (Athey et al., 2021).
- `ESFNNMC.R`: Eigenvector spatial filters nuclear-norm matrix completion.
- `Weight_matrix.R`: Function to create a k-nn spatial weight matrix.
- `Spatial_eigenvectors.R`: Function to select candidate spatial filters.

Recommended loading order:

```r
source("R/FENNMC.R")
source("R/ESFNNMC.R")
source("R/Weight_matrix.R")
source("R/Spatial_eigenvectors.R")
```

