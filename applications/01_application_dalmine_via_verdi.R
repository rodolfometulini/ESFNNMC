# ============================================================
# ESFNNMC application: PM10 imputation for Dalmine - via Verdi
# ============================================================
#
# Purpose:
#   Apply FENNMC and ESFNNMC to the Agrimonia PM10 dataset and
#   compare fitted/imputed time series for the monitoring station
#   Dalmine - via Verdi.
#
# Notes:
#   - Core model functions are assumed to be available from the
#     package/source files in the R/ folder.
#   - Spatial weights are based on 10-nearest neighbours.
#   - This script does not run out-of-sample validation.
#
# Expected functions from the package:
#   - mcnnm_cv_R()
#   - mcnnm_cv_R_with_A()
#   - select_spatial_evec()
#
# Optional recommended helpers, if included in the package:
#   - create_pm10_matrix()
#   - build_knn_weights()
#   - reconstruct_mcfe()
#   - reconstruct_esfnnmc()
# ============================================================


# ---- 0. Packages ----

required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "spdep"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Please install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(spdep)


# ---- 1. User settings ----

# Path to the Agrimonia .RData file.
# The file is expected to contain an object named `a`.
data_path <- "data/Agrimonia_stations.RData"

# Output folder.
out_dir <- "outputs/application_dalmine_via_verdi"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Application settings.
pollutant_var <- "AQ_pm10"
time_var <- "Time"
station_id_var <- "IDStations"
year_selected <- "2021"

# Dalmine - via Verdi.
# In the original application script this corresponds to row 4 / station ID 1269.
target_station_id <- "1269"
target_station_label <- "Dalmine - via Verdi"

# Spatial weights.
k_neighbours <- 10

# Model settings.
num_lam <- 20
num_folds <- 5
cv_ratio <- 0.6
niter <- 200
rel_tol <- 1e-5
is_quiet <- FALSE


# ---- 2. Load data ----

if (!file.exists(data_path)) {
  stop("Data file not found: ", data_path)
}

load(data_path)

if (!exists("a")) {
  stop("The loaded .RData file must contain an object named `a`.")
}

aq_data <- a %>%
  mutate(year = substr(.data[[time_var]], 1, 4)) %>%
  filter(.data$year == year_selected)


# ---- 3. Select PM10 stations ----

# Stations used in the original application.
pm10_station_ids <- c(
  "1264", "1265", "1266", "1269", "1297", "1374",
  "517", "528", "531", "542", "546", "548", "554", "558", "560",
  "561", "564", "565", "569", "571", "572", "574", "576", "583",
  "584", "592", "595", "596", "598", "600", "604", "608", "609",
  "627", "629", "633", "642", "643", "649", "654", "655", "659",
  "661", "663", "664", "669", "670", "673", "674", "677", "679",
  "681", "683", "685", "687", "690", "693", "695", "697", "703",
  "705", "706", "708", "709"
)

aq_pm10 <- aq_data %>%
  filter(.data[[station_id_var]] %in% pm10_station_ids)

# Keep one coordinate pair per station.
station_coords <- aq_pm10 %>%
  distinct(
    IDStations = .data[[station_id_var]],
    Latitude,
    Longitude
  ) %>%
  mutate(IDStations = as.character(IDStations)) %>%
  arrange(match(IDStations, pm10_station_ids))

if (!target_station_id %in% station_coords$IDStations) {
  stop("Target station ID not found in selected stations: ", target_station_id)
}


# ---- 4. Build station-by-time matrix ----

mat_df <- aq_pm10 %>%
  select(
    IDStations = all_of(station_id_var),
    Time = all_of(time_var),
    value = all_of(pollutant_var)
  ) %>%
  mutate(IDStations = as.character(IDStations)) %>%
  pivot_wider(
    names_from = Time,
    values_from = value
  ) %>%
  arrange(match(IDStations, pm10_station_ids))

time_names <- setdiff(names(mat_df), "IDStations")

mat <- as.matrix(mat_df[, time_names])
storage.mode(mat) <- "numeric"
rownames(mat) <- mat_df$IDStations

row <- nrow(mat)
col <- ncol(mat)

message("Matrix size: ", row, " stations x ", col, " time points")
message("Observed cells: ", sum(!is.na(mat)))
message("Missing cells: ", sum(is.na(mat)))


# ---- 5. Create mask and numerical matrix ----

mask <- matrix(as.integer(!is.na(mat)), row, col)

mat_r <- mat
mat_r[is.na(mat_r)] <- 0


# ---- 6. KNN spatial weights, k = 10 ----

coords <- station_coords %>%
  filter(IDStations %in% rownames(mat)) %>%
  arrange(match(IDStations, rownames(mat)))

coords_sp <- cbind(coords$Longitude, coords$Latitude)

knn <- spdep::knearneigh(coords_sp, k = k_neighbours)
nb <- spdep::knn2nb(knn)

W_knn <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
rownames(W_knn) <- coords$IDStations
colnames(W_knn) <- coords$IDStations

# Select spatial eigenvectors.
sel_knn <- select_spatial_evec(W_knn)
A_knn <- sel_knn$A

message("Number of selected spatial eigenvectors: ", ncol(A_knn))
message("Selected Moran's I values:")
print(sel_knn$MoranI)


# ---- 7. Fit FENNMC ----

start_mcfe <- Sys.time()

res_mcfe <- mcnnm_cv_R(
  mat_r,
  mask,
  num_lam = num_lam,
  to_estimate_u = TRUE,
  to_estimate_v = TRUE,
  num_folds = num_folds,
  cv_ratio = cv_ratio,
  niter = niter,
  rel_tol = rel_tol,
  is_quiet = is_quiet
)

pred_mcfe <- res_mcfe$L +
  matrix(rep(res_mcfe$u, col), row, col) +
  t(matrix(rep(res_mcfe$v, row), col, row))

end_mcfe <- Sys.time()
message("FENNMC elapsed time:")
print(end_mcfe - start_mcfe)


# ---- 8. Fit ESFNNMC ----

start_esf <- Sys.time()

res_esf <- mcnnm_cv_R_with_A(
  mat_r,
  mask,
  A_knn,
  num_lam = num_lam,
  to_estimate_alpha = TRUE,
  to_estimate_v = TRUE,
  num_folds = num_folds,
  cv_ratio = cv_ratio,
  niter = niter,
  rel_tol = rel_tol,
  is_quiet = is_quiet
)

pred_esf <- res_esf$L +
  (A_knn %*% res_esf$alpha) %*% matrix(1, 1, ncol(res_esf$L)) +
  matrix(1, nrow(res_esf$L), 1) %*% t(res_esf$v)

end_esf <- Sys.time()
message("ESFNNMC elapsed time:")
print(end_esf - start_esf)


# ---- 9. In-sample fit on observed cells ----

mape_observed <- function(y, yhat, observed_mask) {
  err <- abs((y - yhat) / y)[observed_mask == 1]
  mean(err[is.finite(err)], na.rm = TRUE) * 100
}

mape_mcfe <- mape_observed(mat, pred_mcfe, mask)
mape_esf <- mape_observed(mat, pred_esf, mask)

fit_summary <- data.frame(
  Method = c("FENNMC", "ESFNNMC"),
  MAPE_observed = c(mape_mcfe, mape_esf),
  Best_lambda = c(res_mcfe$best_lambda, res_esf$best_lambda)
)

print(fit_summary)

write.table(
  fit_summary,
  file = file.path(out_dir, "fit_summary_dalmine_via_verdi.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ---- 10. Extract Dalmine - via Verdi time series ----

target_row <- which(rownames(mat) == target_station_id)

if (length(target_row) != 1) {
  stop("Target station ID should match exactly one row.")
}

date_seq <- as.Date(time_names)

series_dalmine <- data.frame(
  Date = date_seq,
  Observed = as.numeric(mat[target_row, ]),
  FENNMC = as.numeric(pred_mcfe[target_row, ]),
  ESFNNMC = as.numeric(pred_esf[target_row, ]),
  Missing = as.integer(is.na(mat[target_row, ]))
)

write.table(
  series_dalmine,
  file = file.path(out_dir, "predictions_dalmine_via_verdi.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ---- 11. Plot Dalmine - via Verdi predictions ----

highlight <- series_dalmine %>%
  mutate(group = cumsum(Missing != lag(Missing, default = 0))) %>%
  filter(Missing == 1) %>%
  group_by(group) %>%
  summarise(
    start = min(Date),
    end = max(Date),
    .groups = "drop"
  )

plot_file <- file.path(out_dir, "PM10_predictions_dalmine_via_verdi.pdf")

pdf(plot_file, width = 12, height = 6.5)

ggplot(series_dalmine, aes(x = Date)) +
  geom_rect(
    data = highlight,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
    fill = "lightblue",
    alpha = 0.4,
    inherit.aes = FALSE
  ) +
  geom_line(
    aes(y = FENNMC, color = "FENNMC"),
    linewidth = 0.5,
    linetype = "solid"
  ) +
  geom_line(
    aes(y = ESFNNMC, color = "ESFNNMC"),
    linewidth = 0.5,
    linetype = "solid"
  ) +
  geom_line(
    aes(y = Observed, color = "Observed"),
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c(
      "FENNMC" = "darkblue",
      "ESFNNMC" = "darkorange3",
      "Observed" = "gray30"
    )
  ) +
  labs(
    title = paste("PM10 predictions -", target_station_label),
    x = NULL,
    y = "PM10 level",
    color = NULL
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = "top"
  )

dev.off()

message("Saved plot: ", plot_file)


# ---- 12. Optional: Moran's I over time ----

listw_knn <- spdep::mat2listw(W_knn, style = "W", zero.policy = TRUE)

moran_values <- apply(mat_r, 2, function(z) {
  spdep::moran.test(z, listw_knn, zero.policy = TRUE)$estimate["Moran I statistic"]
})

moran_df <- data.frame(
  Date = date_seq,
  Moran_I = as.numeric(moran_values)
)

write.table(
  moran_df,
  file = file.path(out_dir, "moran_I_by_day.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

pdf(file.path(out_dir, "moran_I_by_day.pdf"), width = 12, height = 5)

ggplot(moran_df, aes(x = Date, y = Moran_I)) +
  geom_line(linewidth = 0.5) +
  labs(
    title = "Moran's I by day",
    x = NULL,
    y = "Moran's I"
  ) +
  theme_classic(base_size = 16)

dev.off()

message("Application completed.")
