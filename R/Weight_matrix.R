# Spatial weight matrix utilities
#
# These functions build spatial weight matrices from unit coordinates. The
# resulting matrix W can be passed to `select_spatial_evec()` to obtain the
# spatial eigenvector matrix A used by ESFNNMC.

# Build a binary k-nearest-neighbor spatial weight matrix.
#
# Args:
#   coords: numeric matrix or data frame with two columns: longitude and latitude,
#     or any two-dimensional coordinate system.
#   k: number of nearest neighbors.
#   station_ids: optional vector of unit identifiers used as row and column names.
#   symmetrize: if TRUE, make the matrix symmetric by setting W[i, j] = W[j, i] = 1
#     whenever either i is a neighbor of j or j is a neighbor of i.
#
# Returns:
#   A binary spatial weight matrix W.

# NOTE: the function requires spdep

build_knn_weights <- function(coords, k = 10, station_ids = NULL, symmetrize = FALSE) {
  if (!requireNamespace("spdep", quietly = TRUE)) {
    stop("Package 'spdep' is required. Install it with install.packages('spdep').")
  }

  coords <- as.matrix(coords)

  if (!is.numeric(coords) || ncol(coords) != 2) {
    stop("coords must be a numeric matrix or data frame with two coordinate columns.")
  }

  n <- nrow(coords)

  if (k < 1 || k >= n) {
    stop("k must be between 1 and nrow(coords) - 1.")
  }

  knn <- spdep::knearneigh(coords, k = k)
  nb <- spdep::knn2nb(knn)
  W <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)

  if (symmetrize) {
    W <- 1 * ((W + t(W)) > 0)
    diag(W) <- 0
  }

  if (!is.null(station_ids)) {
    if (length(station_ids) != n) {
      stop("station_ids must have length equal to nrow(coords).")
    }
    rownames(W) <- station_ids
    colnames(W) <- station_ids
  }

  W
}