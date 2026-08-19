# Spatial eigenvector selection utilities
#
# These functions construct Moran eigenvectors from a spatial weight matrix W,
# rank them by Moran's I, and select the leading eigenvectors to be used in the
# ESFNNMC spatial component.

## compute Moran's I for each column of A.
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

# Select Moran eigenvectors from a spatial weight matrix W.
#
# Args:
#   W: spatial weight matrix.
#   q: optional number of eigenvectors to retain. If NULL, q is selected using
#   the cumulative share of positive Moran's I signal (default is NULL).
#   only_positive: if TRUE, retain only eigenvectors with positive Moran's I (default is TRUE).
#   threshold: minimum eigenvalue used to retain candidate eigenvectors.
#   explained: cumulative share of Moran's I signal retained when q is NULL (default is 0.90).
#   symmetrize: if TRUE, replace W by (W + t(W)) / 2 before decomposition (default is TRUE).
#
# Returns:
#   A list with the selected eigenvector matrix A, the corresponding Moran's I
#   values, retained eigenvalue indices, and q.

select_spatial_evec <- function(W, q = NULL, only_positive = TRUE, explained = 0.90, threshold = 1e-6, symmetrize = TRUE) {
  if (symmetrize) {
    W <- (W + t(W)) / 2
  }
  
  n <- nrow(W)
  H <- diag(n) - matrix(1, n, n) / n
  C <- H %*% W %*% H
  
  eig <- eigen(C, symmetric = TRUE)
  vals <- eig$values
  vecs <- eig$vectors
  
  # filter eigenvectors by eigenvalue sign
  idx <- which(vals > threshold)
  Afull <- vecs[, idx, drop = FALSE]
  
  # compute Moran's I for each
  Ivals <- moran_I_evec(Afull, W)
  
  # keep only positive Moran I if requested
  if (only_positive) {
    keep <- which(Ivals > 0)
    Afull <- Afull[, keep, drop = FALSE]
    Ivals <- Ivals[keep]
  }
  
  # rank by Moran I magnitude (highest first)
  ord <- order(Ivals, decreasing = TRUE)
  Afull <- Afull[, ord, drop = FALSE]
  Ivals <- Ivals[ord]
  
  # choose q (default: retain enough to explain 90% autocorrelation signal)
  if (is.null(q)) {
    cum <- cumsum(Ivals) / sum(Ivals)
    q <- max(1, which(cum <= explained))
  } else {
    q <- min(q, ncol(Afull))
  }
  
  A <- Afull[, seq_len(q), drop = FALSE]
  
  list(A = A, MoranI = Ivals[seq_len(q)], idx = idx[seq_len(q)], q = q)
}