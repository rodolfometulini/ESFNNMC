# Set your directory
#setwd("")

# packages
install.packages("dplyr")
install.packages("tidyr")
install.packages("tibble")
install.packages("sp")
install.packages("sf")
install.packages("mapview")
install.packages("ggplot2")
install.packages("spdep")
library(dplyr)
library(tidyr)
library(tibble)
library(sp)
library(mapview)
library(sf)
library(ggplot2)
library(spdep)


### load the dataset
load("Data/Agrimonia_stations.RData")

### load codes
source("R/FENNMC.R")
source("R/ESFNNMC.R")
source("R/Weight_matrix.R")
source("R/Spatial_eigenvectors.R")


# check data
head(a)

# extract 2021 daily data
a$year = substr(a$Time,1,4)
ds21 = a[a$year=="2021",]

# Select stations in Lombardy that do not have all PM10 values missing.
options(tibble.print_max = 1000)
ds21 %>%
  group_by(IDStations) %>%
 summarise(missing_pm10 = sum(is.na(AQ_pm10)))

ds21_red = ds21[ds21$IDStations %in% c("1264","1265","1266","1269","1297","1374","517","528","531","542","546","548","554","558","560",
"561","564","565","569","571","572","574","576","583","584","592","595","596","598","600","604","608",
"609","627","629","633","642","643","649","654","655","659","661","663","664","669","670","673",
"674","677","679","681","683","685","687","690","693","695","697","703","705","706","708","709"),1:8]

head(ds21_red, 20)

# from vector to matrix
mat <- ds21_red %>%
  pivot_wider(
    names_from = Time,   # variable whose categories become column names
    values_from = AQ_pm10      # variable whose values fill the cells
  )

### CREATE THE MAP

# Convert tibble to SpatialPointsDataFrame
PM10 = mat
coordinates(PM10) <- ~Longitude + Latitude   # specify coordinate columns
proj4string(PM10) <- CRS("+proj=longlat +datum=WGS84")  # assign CRS
plot(PM10, pch = 19, cex = 1.2)

italy = st_read("Data/georef-italy-regione-millesime.shp")
milan <- st_sfc(st_point(c(9.19, 45.46)),crs = 4326)
milan <- st_transform(milan, st_crs(italy))
lombardy <- italy[st_intersects(italy, milan, sparse = FALSE), ]

# FIGURE 1: map with the localization of the stations and the level of PM10 on January 31 and July 31
# to be visualized in a browser
mapview(PM10, zcol="2021-07-31") + mapview(lombardy, col.regions="transparent", alpha.regions = 0, lwd=3, col="blue", legend=FALSE)
mapview(PM10, zcol="2021-01-31") + mapview(lombardy, col.regions="transparent", alpha.regions = 0, lwd=3, col="blue", legend=FALSE)
# END FIGURE 1

# create matrix to be used
mat = as.matrix(mat[,7:371])
rownames(mat) = c("1264","1265","1266","1269","1297","1374","517","528","531","542","546","548","554","558","560",
"561","564","565","569","571","572","574","576","583","584","592","595","596","598","600","604","608",
"609","627","629","633","642","643","649","654","655","659","661","663","664","669","670","673",
"674","677","679","681","683","685","687","690","693","695","697","703","705","706","708","709")

# draw the missing values
# convert matrix to long format
df <- as.data.frame(as.table(mat))

# FIGURE 2: heatmap
pdf("Data/heatmap.pdf", width = 14, height = 6)
ggplot(df, aes(Var2, Var1, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "green", high = "red", na.value = "black") +
  theme_minimal() +
  labs(x = "Day", y = "Station", fill = "PM10")
dev.off()
# END FIGURE 2


## Models

# ESF

### create the spatial weight matrix
# k-nn, k = 10

df <- data.frame(
  IDStations = c(1264,1265,1266,1269,1297,1374,517,528,531,542,546,548,554,558,560,561,564,565,569,571,572,574,576,583,584,592,595,596,598,600,604,608,609,627,629,633,642,643,649,654,655,659,661,663,664,669,670,673,674,677,679,681,683,685,687,690,693,695,697,703,705,706,708,709),
  Latitude = c(46.2,45.3,45.2,45.6,45.2,45.6,45.5,45.5,45.5,45.7,45.5,45.5,45.6,45.6,45.8,45.8,45.8,45.7,46.2,46.5,46.1,45.9,45.7,45.7,45.7,45.5,45.6,45.6,45.1,45.3,45.3,45.2,45.5,45.1,45.4,45.3,45.2,45.2,45.5,45.6,45.9,45.6,45.5,45.1,45.2,45.5,45.2,45.0,45.6,45.1,45.8,45.9,45.5,45.7,45.6,45.1,45.1,45.4,45.0,45.0,45.5,45.9,45.3,45.3),
  Longitude = c(9.88,9.50,9.67,9.60,9.93,9.28,8.74,9.20,9.33,9.16,8.88,9.20,9.03,8.84,8.82,9.07,9.22,9.13,9.87,10.4,9.57,9.40,9.41,9.64,9.66,9.59,9.56,9.61,9.70,9.49,9.41,9.70,9.56,10.0,9.70,9.86,9.16,9.15,10.2,10.2,10.2,10.4,10.3,10.8,10.8,10.2,10.8,9.01,9.27,10.0,9.35,9.50,9.52,9.48,8.76,8.87,8.90,10.7,11.2,11.1,9.24,9.40,8.75,8.84)
)

# Coordinates
coords_sp <- cbind(df$Longitude, df$Latitude)

# 10-nearest neighbors
knn <- knearneigh(coords_sp, k = 10)
nb  <- knn2nb(knn)

# Convert to weight matrix
W_knn <- nb2mat(nb, style = "B")  # "B" = binary weights (0/1)
W_knn

# Convert df to sf object
stations_sf <- st_as_sf(df, coords = c("Longitude", "Latitude"), crs = 4326)

# Create line geometry from neighbor list
edges <- lapply(1:length(nb), function(i) {
  if(length(nb[[i]]) == 0) return(NULL)
  do.call(st_sfc, lapply(nb[[i]], function(j) {
    st_linestring(rbind(coords_sp[i,], coords_sp[j,]))
  }))
})

# Remove NULLs and combine
edges_sf <- st_as_sf(data.frame(id = 1:length(edges)), geometry = st_sfc(do.call(c, edges)), crs = 4326)

# FIGURE 3: map with valid connections
mapview(stations_sf, col.region = "red", cex = 4, layer.name = "Stations") +
mapview(edges_sf, color = "blue", layer.name = "KNN Links", lwd=1) +
mapview(lombardy, col.regions="transparent", alpha.regions = 0, lwd=3, col="blue", legend=FALSE)
# END FIGURE 3

# Select eigenvectors 
sel <- select_spatial_evec(W_knn)
A <- sel$A
print(sel$MoranI)


#FIGURE 4: plot eigenvectors
#to be visualized using the browser
df$E1 = A[,1]; df$E2 = A[,2]; df$E3 = A[,3];df$E4 = A[,4];df$E5 = A[,5];df$E6 = A[,6];df$E7 = A[,7]
stations_sf <- st_as_sf(df, coords = c("Longitude", "Latitude"), crs = 4326)
# Map first spatial dimension
mapview(stations_sf, zcol="E1", cex = 6, layer.name = "Stations") + mapview(lombardy, col.regions="transparent", alpha.regions = 0, lwd=3, col="blue", legend=FALSE)
# Map second spatial dimension
mapview(stations_sf, zcol="E2", cex = 6, layer.name = "Stations") + mapview(lombardy, col.regions="transparent", alpha.regions = 0, lwd=3, col="blue", legend=FALSE)
# Map third spatial dimension
mapview(stations_sf, zcol="E3", cex = 6, layer.name = "Stations") + mapview(lombardy, col.regions="transparent", alpha.regions = 0, lwd=3, col="blue", legend=FALSE)
# Map fourth spatial dimension
mapview(stations_sf, zcol="E4", cex = 6, layer.name = "Stations") + mapview(lombardy, col.regions="transparent", alpha.regions = 0, lwd=3, col="blue", legend=FALSE)
# Map fifth spatial dimension
mapview(stations_sf, zcol="E5", cex = 6, layer.name = "Stations") + mapview(lombardy, col.regions="transparent", alpha.regions = 0, lwd=3, col="blue", legend=FALSE)
# Map six spatial dimension
mapview(stations_sf, zcol="E6", cex = 6, layer.name = "Stations") + mapview(lombardy, col.regions="transparent", alpha.regions = 0, lwd=3, col="blue", legend=FALSE)
# Map six spatial dimension
mapview(stations_sf, zcol="E7", cex = 6, layer.name = "Stations") + mapview(lombardy, col.regions="transparent", alpha.regions = 0, lwd=3, col="blue", legend=FALSE)

# END FIGURE 4

# compute Moran'I of the original matrix
# 31 January 2021
mat1 = mat[,31]
# replace the missing data with an average of the two near values.
mat1[is.na(mat1)]=33
lw = nb2listw(nb)
moran.test(mat1, lw, zero.policy = FALSE)

# 31 July 2021
mat1 = mat[,212]
# replace the missing data with an average of the two near values.
mat1[16]=15.5
mat1[23]=14.5
mat1[52]=19
mat1[53]=19
lw = nb2listw(nb)
moran.test(mat1, lw, zero.policy = FALSE)


# estimate the models

col = dim(mat)[2]
row = dim(mat)[1]
n = col*row
ns = sqrt(n)
mask =  matrix(as.integer(ifelse(is.na(mat), 0, 1)), row, col)
mat_r = mat
mat_r[is.na(mat_r)] = 999

# ESFNNMC
start <- Sys.time()
fit_esf = esfnnmc(mat_r, mask, A, num_lam = 20,to_estimate_alpha = TRUE, to_estimate_v = TRUE, num_folds = 5, cv_ratio = 0.6, niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
mhat_esf <- compute_matrix_with_A(fit_esf$L, A, fit_esf$alpha, fit_esf$v)
end <- Sys.time()
end - start

# estimated alpha coefficients
print(fit_esf$alpha)

# FENNMC
start <- Sys.time()
fit_fe = fennmc(mat_r, mask, num_lam = 20,to_estimate_u = TRUE, to_estimate_v = TRUE, num_folds = 5, cv_ratio = 0.6, niter = 200, rel_tol = 1e-5, is_quiet = FALSE)
mhat_fe = compute_matrix(fit_fe$L, fit_fe$u, fit_fe$v)
end <- Sys.time()
end - start


# FIGURE 5: Plot of the estimate of the ESFNNMC and FENNMC versus true values for Dalmine via Verdi station.
mhat_fe_dalmine = mhat_fe[4,]
mhat_esf_dalmine = mhat_esf[4,]
true = mat[4,] 
Date <- seq(as.Date("2021-01-01"), by="day", length.out=365)
flag = ifelse(is.na(mat[4,]),1,0)
data <- data.frame(Date, true, mhat_fe_dalmine , mhat_esf_dalmine, flag)

highlight <- data %>%
  mutate(group = cumsum(flag != lag(flag, default = 0))) %>% 
  filter(flag == 1) %>%
  group_by(group) %>%
  summarize(start = min(Date), end = max(Date))

pdf(
  "PM10_predictions2.pdf",
  width = 12,   # inches
  height = 6.5
)

ggplot(data, aes(x = Date)) +
  geom_rect(
    data = highlight,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
    fill = "lightblue", alpha = 0.4, inherit.aes = FALSE
  ) +
geom_line(
  aes(y = mhat_fe_dalmine, color = "FENNMC"),
  linewidth = 1,
  alpha = 0.7
) +
geom_line(
  aes(y = mhat_esf_dalmine, color = "ESFNNMC"),
  linewidth = 0.4
) +
  geom_line(
    aes(y = true, color = "True"),
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c(
      "FENNMC" = "darkblue",
      "ESFNNMC" = "darkorange3",
      "True" = "gray30"
    )
  ) +
  labs(
    color = "Predictions",
    y = "PM10 level"
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.position = "top"
  )
dev.off()

# END FIGURE 5

### Validation strategy
# TABLE 8
# do not run: Large computational time.
set.seed(123)
B <- 200
mask_levels <- c(0.05, 0.10, 0.20)

# Original observed mask
mask_obs <- matrix(as.integer(!is.na(mat)), row, col)

# Replace NA only for numerical stability; mask controls observed cells
mat_r <- mat
mat_r[is.na(mat_r)] <- 0

# Spatial filters
sel_knn <- select_spatial_evec(W_knn)
A_knn <- sel_knn$A
print(sel_knn$MoranI)

# Storage list
results_list <- list()

for (mlev in mask_levels) {
  
  cat("\nArtificial masking level:", mlev, "\n")
  
  mape_esf <- numeric(B)
  mape_mcfe <- numeric(B)
  
  for (b in 1:B) {
    
    cat("Replication", b, "of", B, "\n")
    
    # Observed cells eligible for artificial masking
    obs_id <- which(mask_obs == 1)
    
    # Randomly select validation cells
    n_val <- floor(length(obs_id) * mlev)
    val_id <- sample(obs_id, n_val, replace = FALSE)
    
    # Training mask: original observed cells minus validation cells
    mask_train <- mask_obs
    mask_train[val_id] <- 0
    
    # Validation mask
    mask_val <- matrix(0, row, col)
    mask_val[val_id] <- 1
    
    # -----------------------------
    # ESFNNMC
    # -----------------------------
    
    res_sp <- esfnnmc(
      mat_r, mask_train, A_knn,
      num_lam = 20,
      to_estimate_alpha = TRUE,
      to_estimate_v = TRUE,
      num_folds = 5,
      cv_ratio = 0.6,
      niter = 200,
      rel_tol = 1e-5,
      is_quiet = TRUE
    )
    
    pred_sp <- res_sp$L +
      (A_knn %*% res_sp$alpha) %*% matrix(1, 1, ncol(res_sp$L)) +
      matrix(1, nrow(res_sp$L), 1) %*% t(res_sp$v)
    
    err_sp <- abs((mat - pred_sp) / mat)[mask_val == 1]
    mape_esf[b] <- mean(err_sp[is.finite(err_sp)], na.rm = TRUE) * 100
    
    # -----------------------------
    # FENNMC
    # -----------------------------
    
    res <- fennmc(
      mat_r, mask_train,
      num_lam = 20,
      to_estimate_u = TRUE,
      to_estimate_v = TRUE,
      num_folds = 5,
      cv_ratio = 0.6,
      niter = 200,
      rel_tol = 1e-5,
      is_quiet = TRUE
    )
    
    pred_mcfe <- res$L +
      matrix(rep(res$u, col), row, col) +
      t(matrix(rep(res$v, row), col, row))
    
    err_mcfe <- abs((mat - pred_mcfe) / mat)[mask_val == 1]
    mape_mcfe[b] <- mean(err_mcfe[is.finite(err_mcfe)], na.rm = TRUE) * 100
  }
  
  # Summary table
  summary_tab <- data.frame(
    Masking = paste0(mlev * 100, "%"),
    Method = c("ESFNNMC", "FENNMC"),
    Q1 = c(
      quantile(mape_esf, 0.25, na.rm = TRUE),
      quantile(mape_mcfe, 0.25, na.rm = TRUE)
    ),
    Median = c(
      median(mape_esf, na.rm = TRUE),
      median(mape_mcfe, na.rm = TRUE)
    ),
    Q3 = c(
      quantile(mape_esf, 0.75, na.rm = TRUE),
      quantile(mape_mcfe, 0.75, na.rm = TRUE)
    )
  )
  
  results_list[[paste0("mask_", mlev * 100)]] <- summary_tab
}

# Combine results
final_results <- do.call(rbind, results_list)

# Print to console
print(final_results)

# END TABLE 8

# relative contribution of components (section 4.4)

L = var(as.vector(fit_esf$L)) 
aA = (A %*% fit_esf$alpha) %*% matrix(1, 1, ncol(fit_esf$L))
spat = var(as.vector(aA))
t = var(as.vector(matrix(1, nrow(fit_esf$L), 1) %*% t(fit_esf$v)))
M = L + spat + t
esf = spat/M
temp = t/M
lr = L/M
contrib = c(lr*100,esf*100,temp*100)
names(contrib) <- c("Low-rank", "ESF", "Time")
# relative contributions in ESFNNMC
print(contrib)

L = var(as.vector(fit_fe$L)) 
u = matrix(rep(fit_fe$u,col),row,col)
fe = var(as.vector(u))
t = var(as.vector(t(matrix(rep(fit_fe$v,row),col,row))))
M = L + fe + t
fe = fe/M
temp = t/M
lr = L/M
contrib = c(lr*100,fe*100,temp*100)
names(contrib) <- c("Low-rank", "Unit", "Time")
# relative contributions in FEFNNMC
print(contrib)











