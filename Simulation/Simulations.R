# Set your directory
#setwd("")

# sources
source("R/FENNMC.R")
source("R/ESFNNMC.R")
source("R/Weight_matrix.R")
source("R/Spatial_eigenvectors.R")
source("Simulation/DGP_functions.R")


### Simulation to obtain the data for Table 1, Table 2 and Table A1.
# MAPE across methods for the 10×10 matrix with rank 5 and three regimes.
# The following code prints the summary file "simulation_summary_10_10_rank5.txt" (and a pdf figure for each rho and phi combination) in your current directory.
# Min, Max, Q1, Q3 and standard deviations are provided besides median values.

rho_vals <- c(0, 0.4, 0.8)
phi_vals <- c(0, 0.4, 0.8)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

log_file <- "simulation_summary_10_10_rank5.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 100
    nrow_mat <- 10
    ncol_mat <- 10
    r <- 5
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = ncol_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(50, 25, 15, 10, 5), # harder to estimate
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(3, 6),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- fennmc(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- esfnnmc(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- fennmc(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) Spatial MC without time fixed effects
        res_sp_ntf <- esfnnmc(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss]) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss]) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_allmodels_rank5_rho_", rho_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}

sink()


### Simulation to obtain the data for Table 3, Table A2 and Table A3.
# MAPE across methods for the 10×50 matrix with rank 5 and 3 regimes.
# The following code prints the summary file "simulation_summary_10_50_rank5.txt" (and a pdf figure for each rho and phi combination) in your current directory.
# Min, Max, Q1, Q3, and standard deviations are provided besides median values.

rho_vals <- c(0, 0.4, 0.8)
phi_vals <- c(0, 0.4, 0.8)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

log_file <- "simulation_summary_10_50_rank5.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 500
    nrow_mat <- 10
    ncol_mat <- 50
    r <- 5
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = nrow_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(50, 25, 15, 10, 5), # harder to estimate
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(15, 30),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- fennmc(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- esfnnmc(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- fennmc(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) Spatial MC without time fixed effects
        res_sp_ntf <- esfnnmc(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss]) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss]) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_allmodels_10_50_rank5_rho_", rho_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}

sink()

### Simulation to obtain the data for Table 4, Table A4 and Table A5.
# MAPE across methods for the 50×10 matrix with rank 10 and 3 regimes.
# The following code prints the summary file "simulation_summary_50_10_rank10.txt" (and a pdf figure for each rho and phi combination) in your current directory.
# Min, Max, Q1, Q3, and standard deviations are provided besides median values.

rho_vals <- c(0, 0.4, 0.8)
phi_vals <- c(0, 0.4, 0.8)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

log_file <- "simulation_summary_50_10_rank10.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 500
    nrow_mat <- 50
    ncol_mat <- 10
    r <- 10
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = nrow_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(40, 30, 25, 20, 15, 12, 10, 7, 5, 2),
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(3, 6),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- fennmc(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- esfnnmc(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- fennmc(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) Spatial MC without time fixed effects
        res_sp_ntf <- esfnnmc(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss]) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss]) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_allmodels_50_10_rank10_rho_", rho_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}

sink()



### Simulation to obtain the data for Table 5, Table A6 and Table A7.
# MAPE across methods for the 30×30 matrix with rank 10 and 3 regimes.
# The following code prints the summary file "simulation_summary_30_30_rank10.txt" (and a pdf figure for each rho and phi combination) in your current directory.
# Min, Max, Q1, Q3, and standard deviations are provided besides median values.

rho_vals <- c(0, 0.4, 0.8)
phi_vals <- c(0, 0.4, 0.8)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}


log_file <- "simulation_summary_30_30_rank10.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 900
    nrow_mat <- 30
    ncol_mat <- 30
    r <- 10
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = nrow_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(40, 30, 25, 20, 15, 12, 10, 7, 5, 2),
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(10, 20),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- fennmc(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- esfnnmc(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- fennmc(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) Spatial MC without time fixed effects
        res_sp_ntf <- esfnnmc(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss]) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss]) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_allmodels_30_30_rank10_rho_", rho_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}

sink()


# Codes for producing data for Table 6
## stress test
# Median MAPE for FENNMC and ESFNNMC for different strengths of unit effects and temporal dependence. 
10% and 20% of missing values
# 10 x 10 matrix, rank = 5,  3 regimes
# it generates the file "simulation_summary_10_10_rank5_nonspatial_FE.txt"

unit_fe_sd_vals <- c(0, 5, 15)
phi_vals <- c(0, 0.4, 0.8)

extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

log_file <- file.path("simulation_summary_10_10_rank5_nonspatial_FE.txt")
sink(log_file, split = TRUE)

for (unit_fe_sd in unit_fe_sd_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 100
    nrow_mat <- 10
    ncol_mat <- 10
    r <- 5
    
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)
    
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix used only to construct ESF filters ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = ncol_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data with non-spatial unit fixed effects ---
        lr <- generate_low_rank_fixed_effects(
          n = nrow_mat,
          p = ncol_mat,
          rank = r,
          phi = ar,
          unit_fe_sd = unit_fe_sd,
          noise_sd = 1e-3,
          singular_scale = c(50, 25, 15, 10, 5),
          enforce_positive = TRUE,
          regime_breaks = c(3, 6),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- fennmc(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) ESFNNMC with time fixed effects
        res_sp_tf <- esfnnmc(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- fennmc(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) ESFNNMC without time fixed effects
        res_sp_ntf <- esfnnmc(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss], na.rm = TRUE) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss], na.rm = TRUE) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss], na.rm = TRUE) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss], na.rm = TRUE) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("unit_fe_sd =", unit_fe_sd, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(
      nfilters,
      miss_seq,
      paste0("Number of spatial filters (ncol(A)) - unit_fe_sd=", unit_fe_sd, ", phi=", ar)
    )
    
    # --- dynamic filename ---
    u_str <- gsub("-", "m", gsub("\\.", "_", as.character(unit_fe_sd)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
#    pdf(paste0("sim_missing_allmodels_rank5_nonspatialFE_unitSD_", u_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (unit FE sd =", unit_fe_sd, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}


sink()




# Simulation to obtain the data for Table 7.

# MAPE across methods for the 10×10 matrix with rank 5 and 3 regimes. Parameter tau = 0.80.
# The following code prints the summary file "simulation_summary_10_10_rank5_tau80.txt" 
# Min, Max, Q1, Q3, and standard deviations are provided besides median values.

rho_vals <- c(0.4)
phi_vals <- c(0.4)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

log_file <- "simulation_summary_10_10_rank5_tau80.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 100
    nrow_mat <- 10
    ncol_mat <- 10
    r <- 5
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = ncol_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(50, 25, 15, 10, 5), # harder to estimate
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(3, 6),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W, explained = 0.80)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- fennmc(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- esfnnmc(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- fennmc(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) Spatial MC without time fixed effects
        res_sp_ntf <- esfnnmc(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss]) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss]) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_allmodels_rank5_tau80_rho_", rho_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}

sink()

# MAPE across methods for the 10×10 matrix with rank 5 and 3 regimes. Parameter tau = 0.90.
# The following code prints the summary file "simulation_summary_10_10_rank5_tau90.txt" 
# Min, Max, Q1, Q3, and standard deviations are provided besides median values.

rho_vals <- c(0.4)
phi_vals <- c(0.4)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

log_file <- "simulation_summary_10_10_rank5_tau90.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 100
    nrow_mat <- 10
    ncol_mat <- 10
    r <- 5
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = ncol_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(50, 25, 15, 10, 5), # harder to estimate
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(3, 6),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W, explained = 0.90)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- fennmc(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- esfnnmc(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- fennmc(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) Spatial MC without time fixed effects
        res_sp_ntf <- esfnnmc(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss]) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss]) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_allmodels_rank5_tau80_rho_", rho_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}

sink()


# MAPE across methods for the 10×10 matrix with rank 5 and 3 regimes. Parameter tau = 0.95.
# The following code prints the summary file "simulation_summary_10_10_rank5_tau95.txt" 
# Min, Max, Q1, Q3, and standard deviations are provided besides median values.

rho_vals <- c(0.4)
phi_vals <- c(0.4)

# helper to extract selected lambda from model output
extract_lambda_L <- function(res) {
  candidate_names <- c("best_lambda", "lambda_L", "lam_L", "lambda", "opt_lambda", "lambda_opt")
  for (nm in candidate_names) {
    if (!is.null(res[[nm]]) && length(res[[nm]]) == 1 && is.numeric(res[[nm]])) {
      return(res[[nm]])
    }
  }
  return(NA_real_)
}

# helper to print summary stats by missingness level
print_summary_table <- function(mat, miss_seq, label) {
  summary_mat <- t(apply(mat, 1, function(x) {
    c(
      Min    = min(x, na.rm = TRUE),
      Q1     = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean   = mean(x, na.rm = TRUE),
      Q3     = quantile(x, 0.75, na.rm = TRUE),
      Max    = max(x, na.rm = TRUE),
      SD     = sd(x, na.rm = TRUE)
    )
  }))
  
  summary_df <- data.frame(
    Missing = miss_seq,
    summary_mat,
    row.names = NULL
  )
  
  cat("\n", label, "\n", sep = "")
  print(summary_df, digits = 4)
}

log_file <- "simulation_summary_10_10_rank5_tau95.txt"
sink(log_file, split = TRUE)

for (rho in rho_vals) {
  for (ar in phi_vals) {
    
    start <- Sys.time()
    
    B <- 200
    miss_seq <- c(0, 2, 4, 6, 8, 10, 15, 20, 25)
    C <- length(miss_seq)
    Wspar <- 0.3
    n <- 100
    nrow_mat <- 10
    ncol_mat <- 10
    r <- 5
    
    # MAPE containers
    mape_mcfe_tf   <- matrix(NA, C, B)
    mape_sp_tf     <- matrix(NA, C, B)
    mape_mcfe_ntf  <- matrix(NA, C, B)
    mape_sp_ntf    <- matrix(NA, C, B)
    
    # lambda containers
    lambda_mcfe_tf   <- matrix(NA, C, B)
    lambda_sp_tf     <- matrix(NA, C, B)
    lambda_mcfe_ntf  <- matrix(NA, C, B)
    lambda_sp_ntf    <- matrix(NA, C, B)

    # filters container
    nfilters <- matrix(NA, C, B)
    
    set.seed(180326)
    
    for (i in seq_along(miss_seq)) {
      miss_rate <- miss_seq[i] / 100
      
      for (j in 1:B) {
        
        # --- spatial weights matrix ---
        W <- matrix(rbinom(n, 1, Wspar), nrow = nrow_mat, ncol = ncol_mat)
        W[lower.tri(W)] <- t(W)[lower.tri(W)]
        diag(W) <- 0
        
        # --- generate data ---
        lr <- generate_low_rank_spatiotemp(
          n = nrow_mat, p = ncol_mat, rank = r, W = W,
          rho = rho, phi = ar, noise_sd = 1e-3,
          singular_scale = c(50, 25, 15, 10, 5), # harder to estimate
          row_standardize = TRUE, check_stability = TRUE,
          enforce_positive = TRUE,
          regime_breaks = c(3, 6),
          regime_means = c(0, 15, 30),
          regime_sd = 1
        )
        mat <- lr$M
        
        # --- missingness mask ---
        mask <- matrix(rbinom(n, 1, 1 - miss_rate), nrow_mat, ncol_mat)
        
        # --- spatial eigenvectors ---
        sel <- select_spatial_evec(W, explained = 0.95)
        A <- sel$A
        nfilters[i, j] <- ncol(A)
        
        # 1) MCFE with time fixed effects
        res_tf <- fennmc(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_tf <- res_tf$L +
          matrix(rep(res_tf$u, ncol_mat), nrow_mat, ncol_mat) +
          t(matrix(rep(res_tf$v, nrow_mat), ncol_mat, nrow_mat))
        
        lambda_mcfe_tf[i, j] <- extract_lambda_L(res_tf)
        
        # 2) Spatial MC with time fixed effects
        res_sp_tf <- esfnnmc(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = TRUE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_tf <- res_sp_tf$L +
          (A %*% res_sp_tf$alpha) %*% matrix(1, 1, ncol(res_sp_tf$L)) +
          matrix(1, nrow(res_sp_tf$L), 1) %*% t(res_sp_tf$v)
        
        lambda_sp_tf[i, j] <- extract_lambda_L(res_sp_tf)
        
        # 3) MCFE without time fixed effects
        res_ntf <- fennmc(
          mat, mask,
          num_lam = 20,
          to_estimate_u = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_ntf <- res_ntf$L +
          matrix(rep(res_ntf$u, ncol_mat), nrow_mat, ncol_mat)
        
        lambda_mcfe_ntf[i, j] <- extract_lambda_L(res_ntf)
        
        # 4) Spatial MC without time fixed effects
        res_sp_ntf <- esfnnmc(
          mat, mask, A,
          num_lam = 20,
          to_estimate_alpha = TRUE,
          to_estimate_v = FALSE,
          num_folds = 5,
          cv_ratio = 0.6,
          niter = 200,
          rel_tol = 1e-5,
          is_quiet = FALSE
        )
        
        pred_sp_ntf <- res_sp_ntf$L +
          (A %*% res_sp_ntf$alpha) %*% matrix(1, 1, ncol(res_sp_ntf$L))
        
        lambda_sp_ntf[i, j] <- extract_lambda_L(res_sp_ntf)
        
        # --- performance ---
        idx_miss <- (mask == 0)
        mape_mcfe_tf[i, j]  <- mean(abs((mat - pred_tf)     / mat)[idx_miss]) * 100
        mape_sp_tf[i, j]    <- mean(abs((mat - pred_sp_tf)  / mat)[idx_miss]) * 100
        mape_mcfe_ntf[i, j] <- mean(abs((mat - pred_ntf)    / mat)[idx_miss]) * 100
        mape_sp_ntf[i, j]   <- mean(abs((mat - pred_sp_ntf) / mat)[idx_miss]) * 100
      }
    }
    
    end <- Sys.time()
    elapsed_sec <- as.numeric(difftime(end, start, units = "secs"))
    hours <- floor(elapsed_sec / 3600)
    mins  <- floor((elapsed_sec %% 3600) / 60)
    secs  <- round(elapsed_sec %% 60, 2)
    
    cat("\n====================================================\n")
    cat("rho =", rho, "| phi =", ar, "\n")
    cat("Computation time:", hours, "h", mins, "min", secs, "sec\n")
    cat("====================================================\n")
    
    # --- names ---
    rownames(mape_mcfe_tf)  <- miss_seq
    rownames(mape_sp_tf)    <- miss_seq
    rownames(mape_mcfe_ntf) <- miss_seq
    rownames(mape_sp_ntf)   <- miss_seq
    
    colnames(mape_mcfe_tf)  <- 1:B
    colnames(mape_sp_tf)    <- 1:B
    colnames(mape_mcfe_ntf) <- 1:B
    colnames(mape_sp_ntf)   <- 1:B
    
    # --- summaries for median lines ---
    med_mcfe_tf  <- apply(mape_mcfe_tf, 1, median, na.rm = TRUE)
    med_sp_tf    <- apply(mape_sp_tf, 1, median, na.rm = TRUE)
    med_mcfe_ntf <- apply(mape_mcfe_ntf, 1, median, na.rm = TRUE)
    med_sp_ntf   <- apply(mape_sp_ntf, 1, median, na.rm = TRUE)
    
    # --- confidence bands (25th and 75th percentiles) ---
    q25_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_tf  <- apply(mape_mcfe_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_tf    <- apply(mape_sp_tf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_mcfe_ntf <- apply(mape_mcfe_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    q25_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.25, na.rm = TRUE)
    q75_sp_ntf   <- apply(mape_sp_ntf, 1, quantile, probs = 0.75, na.rm = TRUE)
    
    # --- summaries written to txt file ---
    print_summary_table(mape_mcfe_tf,  miss_seq, "MAPE - MCFE with time FE")
    print_summary_table(mape_sp_tf,    miss_seq, "MAPE - ESF-NNMC with time FE")
    print_summary_table(mape_mcfe_ntf, miss_seq, "MAPE - MCFE without time FE")
    print_summary_table(mape_sp_ntf,   miss_seq, "MAPE - ESF-NNMC without time FE")
    
    print_summary_table(lambda_mcfe_tf,  miss_seq, "lambda_L - MCFE with time FE")
    print_summary_table(lambda_sp_tf,    miss_seq, "lambda_L - ESF-NNMC with time FE")
    print_summary_table(lambda_mcfe_ntf, miss_seq, "lambda_L - MCFE without time FE")
    print_summary_table(lambda_sp_ntf,   miss_seq, "lambda_L - ESF-NNMC without time FE")
    
    print_summary_table(nfilters, miss_seq, paste0("Number of spatial filters (ncol(A)) - rho=", rho, ", phi=", ar))
    # --- dynamic filename ---
    rho_str <- gsub("-", "m", gsub("\\.", "_", as.character(rho)))
    phi_str <- gsub("-", "m", gsub("\\.", "_", as.character(ar)))
    
    pdf(paste0("sim_missing_allmodels_rank5_tau80_rho_", rho_str, "_phi_", phi_str, ".pdf"))
    
    x <- miss_seq
    
    ymax <- max(
      q75_mcfe_tf, q75_sp_tf, q75_mcfe_ntf, q75_sp_ntf,
      med_mcfe_tf, med_sp_tf, med_mcfe_ntf, med_sp_ntf,
      na.rm = TRUE
    )
    
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    plot(
      x, med_mcfe_tf, type = "n",
      xlab = "Percentage of missings",
      ylab = "MAPE",
      main = paste("MAPE comparison (rho =", rho, ", phi =", ar, ")"),
      ylim = c(0, ymax * 1.02)
    )
    
    # --- dotted confidence bands ---
    lines(x, q25_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_tf,  col = "blue",      lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    lines(x, q75_sp_tf,    col = "red",       lty = 3, lwd = 1.2)
    
    lines(x, q25_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    lines(x, q75_mcfe_ntf, col = "darkgreen", lty = 3, lwd = 1.2)
    
    lines(x, q25_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    lines(x, q75_sp_ntf,   col = "purple",    lty = 3, lwd = 1.2)
    
    # --- all median lines solid ---
    lines(x, med_mcfe_tf,  type = "b", pch = 19, col = "blue",      lwd = 2, lty = 1)
    lines(x, med_sp_tf,    type = "b", pch = 19, col = "red",       lwd = 2, lty = 1)
    lines(x, med_mcfe_ntf, type = "b", pch = 17, col = "darkgreen", lwd = 2, lty = 1)
    lines(x, med_sp_ntf,   type = "b", pch = 17, col = "purple",    lwd = 2, lty = 1)
    
    legend(
      "topleft",
      legend = c(
        "MCFE with time FE",
        "ESF-NNMC with time FE",
        "MCFE without time FE",
        "ESF-NNMC without time FE",
        "25%-75% band"
      ),
      col = c("blue", "red", "darkgreen", "purple", "black"),
      lwd = c(2, 2, 2, 2, 1.2),
      lty = c(1, 1, 1, 1, 3),
      pch = c(19, 19, 17, 17, NA),
      bty = "n"
    )
    
    dev.off()
  }
}

sink()



