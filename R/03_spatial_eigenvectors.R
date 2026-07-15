# Spatial eigenvector selection utilities
#
# These functions construct Moran eigenvectors from a spatial weight matrix W,
# rank them by Moran's I, and select the leading eigenvectors to be used in the
# ESFNNMC spatial component.

select_spatial_evec80 <- function(W, q = NULL, only_positive = TRUE, threshold = 1e-6) {
  # Symmetrize W for numerical stability.
  W <- (W + t(W)) / 2
  n <- nrow(W)
  H <- diag(n) - matrix(1, n, n) / n
  C <- H %*% W %*% H
  
  eig <- eigen(C, symmetric = TRUE)
  vals <- eig$values
  vecs <- eig$vectors

  # Retain eigenvectors associated with positive eigenvalues.
  idx <- which(vals > threshold)
  Afull <- vecs[, idx, drop = FALSE]

  # Compute Moran's I for each candidate eigenvector.
  Ivals <- moran_I_evec(Afull, W)

  # Optionally retain only eigenvectors with positive Moran's I.
  if (only_positive) {
    keep <- which(Ivals > 0)
    Afull <- Afull[, keep, drop = FALSE]
    Ivals <- Ivals[keep]
  }

  # Rank eigenvectors by Moran's I in decreasing order.
  ord <- order(Ivals, decreasing = TRUE)
  Afull <- Afull[, ord, drop = FALSE]
  Ivals <- Ivals[ord]

  # Choose q; by default retain eigenvectors explaining about 80% of the positive spatial autocorrelation signal.
  if (is.null(q)) {
    cum <- cumsum(Ivals) / sum(Ivals)
    q <- max(1, which(cum <= 0.80))
  } else {
    q <- min(q, ncol(Afull))
  }

  A <- Afull[, seq_len(q), drop = FALSE]

  list(A = A, MoranI = Ivals[seq_len(q)], idx = idx[seq_len(q)], q = q)
}

select_spatial_evec <- function(W, q = NULL, only_positive = TRUE, threshold = 1e-6) {
  # Symmetrize W for numerical stability.
  W <- (W + t(W)) / 2
  n <- nrow(W)
  H <- diag(n) - matrix(1, n, n) / n
  C <- H %*% W %*% H
  
  eig <- eigen(C, symmetric = TRUE)
  vals <- eig$values
  vecs <- eig$vectors

  # Retain eigenvectors associated with positive eigenvalues.
  idx <- which(vals > threshold)
  Afull <- vecs[, idx, drop = FALSE]

  # Compute Moran's I for each candidate eigenvector.
  Ivals <- moran_I_evec(Afull, W)

  # Optionally retain only eigenvectors with positive Moran's I.
  if (only_positive) {
    keep <- which(Ivals > 0)
    Afull <- Afull[, keep, drop = FALSE]
    Ivals <- Ivals[keep]
  }

  # Rank eigenvectors by Moran's I in decreasing order.
  ord <- order(Ivals, decreasing = TRUE)
  Afull <- Afull[, ord, drop = FALSE]
  Ivals <- Ivals[ord]

  # Choose q; by default retain eigenvectors explaining about 90% of the positive spatial autocorrelation signal.
  if (is.null(q)) {
    cum <- cumsum(Ivals) / sum(Ivals)
    q <- max(1, which(cum <= 0.90))
  } else {
    q <- min(q, ncol(Afull))
  }

  A <- Afull[, seq_len(q), drop = FALSE]

  list(A = A, MoranI = Ivals[seq_len(q)], idx = idx[seq_len(q)], q = q)
}

select_spatial_evec95 <- function(W, q = NULL, only_positive = TRUE, threshold = 1e-6) {
  # Symmetrize W for numerical stability.
  W <- (W + t(W)) / 2
  n <- nrow(W)
  H <- diag(n) - matrix(1, n, n) / n
  C <- H %*% W %*% H
  
  eig <- eigen(C, symmetric = TRUE)
  vals <- eig$values
  vecs <- eig$vectors

  # Retain eigenvectors associated with positive eigenvalues.
  idx <- which(vals > threshold)
  Afull <- vecs[, idx, drop = FALSE]

  # Compute Moran's I for each candidate eigenvector.
  Ivals <- moran_I_evec(Afull, W)

  # Optionally retain only eigenvectors with positive Moran's I.
  if (only_positive) {
    keep <- which(Ivals > 0)
    Afull <- Afull[, keep, drop = FALSE]
    Ivals <- Ivals[keep]
  }

  # Rank eigenvectors by Moran's I in decreasing order.
  ord <- order(Ivals, decreasing = TRUE)
  Afull <- Afull[, ord, drop = FALSE]
  Ivals <- Ivals[ord]

  # Choose q; by default retain eigenvectors explaining about 95% of the positive spatial autocorrelation signal.
  if (is.null(q)) {
    cum <- cumsum(Ivals) / sum(Ivals)
    q <- max(1, which(cum <= 0.95))
  } else {
    q <- min(q, ncol(Afull))
  }

  A <- Afull[, seq_len(q), drop = FALSE]

  list(A = A, MoranI = Ivals[seq_len(q)], idx = idx[seq_len(q)], q = q)
}



# Compute Moran's I for each column of A.
moran_I_evec <- function(A, W) {
  n <- nrow(A)
  # Row-standardize W.
  rs <- rowSums(W)
  Wrs <- W / ifelse(rs == 0, 1, rs)

  S0 <- sum(Wrs)

  Ivals <- sapply(seq_len(ncol(A)), function(j) {
    a <- A[, j]
    a <- a - mean(a)
    num <- t(a) %*% Wrs %*% a
    den <- sum(a^2)
    as.numeric((n / S0) * (num / den))
  })
  Ivals
}


