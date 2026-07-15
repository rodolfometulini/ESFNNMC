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

# Build a binary distance-band spatial weight matrix.
#
# Args:
#   coords: numeric matrix or data frame with two columns: longitude and latitude.
#   distance_threshold_km: maximum distance, in kilometers, for two units to be
#     considered neighbors.
#   station_ids: optional vector of unit identifiers used as row and column names.
#   distance_fun: distance function passed to geosphere::distm().
#
# Returns:
#   A binary spatial weight matrix W.
build_distance_band_weights <- function(coords,
                                        distance_threshold_km,
                                        station_ids = NULL,
                                        distance_fun = NULL) {
  if (!requireNamespace("geosphere", quietly = TRUE)) {
    stop("Package 'geosphere' is required. Install it with install.packages('geosphere').")
  }

  coords <- as.matrix(coords)

  if (!is.numeric(coords) || ncol(coords) != 2) {
    stop("coords must be a numeric matrix or data frame with two columns: longitude and latitude.")
  }

  if (missing(distance_threshold_km) || distance_threshold_km <= 0) {
    stop("distance_threshold_km must be a positive number.")
  }

  if (is.null(distance_fun)) {
    distance_fun <- geosphere::distHaversine
  }

  dist_m <- geosphere::distm(coords, coords, fun = distance_fun)
  dist_km <- dist_m / 1000

  W <- ifelse(dist_km < distance_threshold_km & dist_km > 0, 1, 0)
  diag(W) <- 0

  if (!is.null(station_ids)) {
    if (length(station_ids) != nrow(coords)) {
      stop("station_ids must have length equal to nrow(coords).")
    }
    rownames(W) <- station_ids
    colnames(W) <- station_ids
  }

  W
}
