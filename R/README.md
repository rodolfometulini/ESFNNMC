# R source files

This codes allows to perform fixed-effects nuclear-norm matrix completion via Soft-Impute algorithm, as described in Athey et al. (2021), 
and the new eigenvector spatial filters nuclear-norm matrix completion. 

Athey et al. codes in the R package MCPanel (https://github.com/susanathey/MCPanel) are translated from C++ to R, 
and modified by replacing unit fixed effects with a set of eigenvectors spatial filters.

Files included in the `R/` folder:

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

# References

```text
Athey, S., Bayati, M., Doudchenko, N., Imbens, G., & Khosravi, K. (2021). 
Matrix completion methods for causal panel data models. 
Journal of the American Statistical Association, 116(536), 1716-1730.
```

