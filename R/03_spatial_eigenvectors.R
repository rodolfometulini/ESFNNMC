# Spatial eigenvector selection utilities
#
# These functions construct Moran eigenvectors from a spatial weight matrix W,
# rank them by Moran's I, and select the leading eigenvectors to be used in the
# ESFNNMC spatial component.

# Compute Moran's I for each column of A.
moran_I_evec <- function(A, W) {
  if (!is.matrix(A)) A <- as.matrix(A)
  if (!is.matrix(W)) W <- as.matrix(W)

  if (nrow(A) != nrow(W) || nrow(W) != ncol(W)) {
    stop("A and W have incompatible dimensions.")
  }

  n <- nrow(A)

  # Row-standardize W.
  rs <- rowSums(W)
  Wrs <- W / ifelse(rs == 0, 1, rs)
  S0 <- sum(Wrs)

  if (S0 == 0) {
    stop("The spatial weight matrix has no positive links.")
  }

  sapply(seq_len(ncol(A)), function(j) {
    a <- A[, j]
    a <- a - mean(a)

    num <- t(a) %*% Wrs %*% a
    den <- sum(a^2)

    as.numeric((n / S0) * (num / den))
  })
}

# Select Moran eigenvectors from a spatial weight matrix.
#
# Args:
#   W: spatial weight matrix.
#   q: optional number of eigenvectors to retain. If NULL, q is selected using
#     the cumulative share of positive Moran's I signal.
#   only_positive: if TRUE, retain only eigenvectors with positive Moran's I.
#   threshold: minimum eigenvalue used to retain candidate eigenvectors.
#   explained: cumulative share of Moran's I signal retained when q is NULL.
#   symmetrize: if TRUE, replace W by (W + t(W)) / 2 before decomposition.
#
# Returns:
#   A list with the selected eigenvector matrix A, the corresponding Moran's I
#   values, retained eigenvalue indices, and q.
select_spatial_evec <- function(W,
                                q = NULL,
                                only_positive = TRUE,
                                threshold = 1e-6,
                                explained = 0.90,
                                symmetrize = TRUE) {
  if (!is.matrix(W)) W <- as.matrix(W)

  if (nrow(W) != ncol(W)) {
    stop("W must be a square matrix.")
  }

  if (explained <= 0 || explained > 1) {
    stop("explained must be in the interval (0, 1].")
  }

  if (symmetrize) {
    W <- (W + t(W)) / 2
  }

  n <- nrow(W)
  H <- diag(n) - matrix(1, n, n) / n
  C <- H %*% W %*% H

  eig <- eigen(C, symmetric = TRUE)
  vals <- eig$values
  vecs <- eig$vectors

  # Retain eigenvectors associated with sufficiently positive eigenvalues.
  idx <- which(vals > threshold)

  if (length(idx) == 0L) {
    stop("No eigenvectors retained. Consider lowering `threshold` or checking W.")
  }

  Afull <- vecs[, idx, drop = FALSE]

  # Compute Moran's I for each candidate eigenvector.
  Ivals <- moran_I_evec(Afull, W)

  # Optionally retain only eigenvectors with positive Moran's I.
  if (only_positive) {
    keep <- which(Ivals > 0)

    if (length(keep) == 0L) {
      stop("No eigenvectors with positive Moran's I were retained.")
    }

    Afull <- Afull[, keep, drop = FALSE]
    Ivals <- Ivals[keep]
    idx <- idx[keep]
  }

  # Rank eigenvectors by Moran's I in decreasing order.
  ord <- order(Ivals, decreasing = TRUE)
  Afull <- Afull[, ord, drop = FALSE]
  Ivals <- Ivals[ord]
  idx <- idx[ord]

  # Choose q; by default retain enough eigenvectors to explain the requested
  # share of the positive spatial autocorrelation signal.
  if (is.null(q)) {
    cum <- cumsum(Ivals) / sum(Ivals)
    q <- which(cum >= explained)[1]
  } else {
    q <- min(q, ncol(Afull))
  }

  A <- Afull[, seq_len(q), drop = FALSE]

  list(
    A = A,
    MoranI = Ivals[seq_len(q)],
    idx = idx[seq_len(q)],
    q = q
  )
}

# Backward-compatible wrappers for fixed cumulative-selection thresholds.
select_spatial_evec80 <- function(W, q = NULL, only_positive = TRUE, threshold = 1e-6) {
  select_spatial_evec(W, q = q, only_positive = only_positive, threshold = threshold, explained = 0.80)
}

select_spatial_evec95 <- function(W, q = NULL, only_positive = TRUE, threshold = 1e-6) {
  select_spatial_evec(W, q = q, only_positive = only_positive, threshold = threshold, explained = 0.95)
}
