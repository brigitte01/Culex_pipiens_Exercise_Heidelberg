# ============================================================
# SECTION 0: PACKAGES AND LIBRARIES
# ============================================================
library(dplyr)
library(sf)
library(terra)
library(ks)
library(readr)

# ============================================================
# SECTION 1: KDE-WEIGHTED BACKGROUND SAMPLING
# ============================================================


set.seed(123)

# Number of background points
n_bg <- 2000

# Environmental predictors
env <- environmental_stack_italy_aligned

# ------------------------------------------------------------
# 1.1 Prepare all surveyed locations
# ------------------------------------------------------------

surveyed_cells_italy <- culex_cells_prepared %>%
  select(
    cell_id,
    longitude,
    latitude
  ) %>%
  distinct(cell_id, .keep_all = TRUE)

# Check dataset
cat(
  "Total surveyed grid cells:",
  nrow(surveyed_cells_italy),
  "\n"
)

# Convert to spatial points
surveyed_sf <- sf::st_as_sf(
  surveyed_cells_italy,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

# Project to a metric coordinate system
surveyed_laea <- sf::st_transform(
  surveyed_sf,
  crs = 3035
)

# Extract projected coordinates
surveyed_xy_m <- sf::st_coordinates(surveyed_laea)

head(surveyed_xy_m)

#--- Exporting the 728 surveyed locations as a separate prepared CSV.
readr::write_csv(
  surveyed_cells_italy,
  file.path(
    "Data",
    "Italy",
    "Presence_culex",
    "Culex_pipiens_Surveyed_Cells_Italy.csv"
  )
)
# ------------------------------------------------------------
# 1.2 Identify candidate background cells
# ------------------------------------------------------------

# Identify complete environmental cells
valid_env <- terra::app(
  env,
  fun = function(x) as.integer(all(!is.na(x)))
)

# Extract valid raster cell IDs
valid_cells <- which(
  terra::values(valid_env, mat = FALSE) == 1
)

# Identify presence raster cells
presence_cells <- terra::cellFromXY(
  env[[1]],
  as.matrix(
    culex_presences_italy[, c("longitude", "latitude")]
  )
)

# Exclude presence cells
candidate_cells <- setdiff(
  valid_cells,
  presence_cells
)

# Summary
cat("Valid environmental cells:", length(valid_cells), "\n")
cat("Presence cells:", length(unique(presence_cells)), "\n")
cat("Background candidate cells:", length(candidate_cells), "\n")

# ============================================================
# SECTION 2: ESTIMATE KDE BANDWIDTH
# ============================================================

library(ks)

# Estimate bandwidth matrix using all surveyed locations
H_kde <- ks::Hpi(
  x = surveyed_xy_m
)

# Display bandwidth matrix
print(H_kde)

# Approximate bandwidth along principal axes (in km)
bandwidth_km <- sqrt(
  eigen(H_kde, symmetric = TRUE)$values
) / 1000

cat(
  "KDE principal-axis bandwidths (km):",
  round(bandwidth_km, 2),
  "\n"
)

# ============================================================
# SECTION 3: CALCULATE KDE SAMPLING-BIAS SURFACE
# ============================================================


# ------------------------------------------------------------
# 3.1 Extract centres of valid environmental cells
# ------------------------------------------------------------

valid_xy <- terra::xyFromCell(
  env,
  valid_cells
)

# Convert cell centres to sf points
valid_sf <- sf::st_as_sf(
  data.frame(
    longitude = valid_xy[, 1],
    latitude = valid_xy[, 2]
  ),
  coords = c("longitude", "latitude"),
  crs = 4326
)

# Project to EPSG:3035
valid_laea <- sf::st_transform(
  valid_sf,
  crs = 3035
)

# Extract projected coordinates
valid_xy_m <- sf::st_coordinates(valid_laea)

# ------------------------------------------------------------
# 3.2 Evaluate KDE at valid environmental cells
# ------------------------------------------------------------

kde_result <- ks::kde(
  x = surveyed_xy_m,
  H = H_kde,
  eval.points = valid_xy_m
)

# Extract KDE density values
kde_weights <- as.numeric(kde_result$estimate)

# Inspect KDE values
summary(kde_weights)

# Check for invalid values
stopifnot(
  length(kde_weights) == length(valid_cells),
  all(is.finite(kde_weights)),
  all(kde_weights >= 0)
)

# ============================================================
# SECTION 4: CREATE KDE BIAS RASTER
# ============================================================

# Create empty raster using the environmental template
kde_bias_italy <- env[[1]]

# Set all cells to NA
terra::values(kde_bias_italy) <- NA_real_

# Assign KDE values to valid environmental cells
kde_bias_italy[valid_cells] <- kde_weights

# Rename layer
names(kde_bias_italy) <- "KDE_sampling_bias"

# Inspect
print(kde_bias_italy)

# Plot KDE bias surface
plot(
  kde_bias_italy,
  main = "KDE sampling-bias surface — Italy"
)

# Overlay surveyed locations
points(
  surveyed_cells_italy$longitude,
  surveyed_cells_italy$latitude,
  pch = 16,
  cex = 0.3,
  col = "black"
)

# ============================================================
# SECTION 5: KDE-WEIGHTED BACKGROUND SAMPLING
# ============================================================

# Extract KDE weights for candidate background cells
candidate_weights <- terra::values(
  kde_bias_italy,
  mat = FALSE
)[candidate_cells]

# Verify sampling weights
stopifnot(
  length(candidate_weights) == length(candidate_cells),
  all(is.finite(candidate_weights)),
  all(candidate_weights >= 0),
  sum(candidate_weights) > 0,
  n_bg <= length(candidate_cells)
)

# Normalize weights
candidate_weights <- candidate_weights / sum(candidate_weights)

# Sample background cells
set.seed(123)

bg_cells <- sample(
  candidate_cells,
  size = n_bg,
  replace = FALSE,
  prob = candidate_weights
)

# Extract geographical coordinates
bg_xy <- terra::xyFromCell(
  env,
  bg_cells
)

# Create background dataset
background_italy <- data.frame(
  raster_cell = bg_cells,
  longitude = bg_xy[, 1],
  latitude = bg_xy[, 2]
)

# Verify background sample
cat("Number of background points:", nrow(background_italy), "\n")

cat(
  "Number of unique background cells:",
  length(unique(background_italy$raster_cell)),
  "\n"
)

cat(
  "Overlap with presence cells:",
  sum(background_italy$raster_cell %in% presence_cells),
  "\n"
)

head(background_italy)

# ============================================================
# SECTION 6: VISUALIZE PRESENCE AND BACKGROUND LOCATIONS
# ============================================================

plot(
  kde_bias_italy,
  main = "KDE-weighted background sampling — Italy"
)

# Background points
points(
  background_italy$longitude,
  background_italy$latitude,
  pch = 16,
  cex = 0.25,
  col = "blue"
)

# Presence points
points(
  culex_presences_italy$longitude,
  culex_presences_italy$latitude,
  pch = 16,
  cex = 0.4,
  col = "red"
)

legend(
  "bottomleft",
  legend = c("Presence", "Background"),
  col = c("red", "blue"),
  pch = 16,
  bty = "n"
)

# ============================================================
# SECTION 8: EXPORT BACKGROUND POINTS AND KDE BIAS RASTER
# ============================================================

# ------------------------------------------------------------
# 8.1 Save background points
# ------------------------------------------------------------

background_output_dir <- file.path(
  "Data", "Italy", "Background"
)


readr::write_csv(
  background_italy,
  file.path(
    background_output_dir,
    "Culex_pipiens_Background_Italy.csv"
  )
)

# ------------------------------------------------------------
# 8.2 Save KDE sampling-bias raster
# ------------------------------------------------------------

terra::writeRaster(
  kde_bias_italy,
  file.path(
    background_output_dir,
    "KDE_Sampling_Bias_Italy.tif"
  ),
  overwrite = TRUE,
  datatype = "FLT8S",
  gdal = "COMPRESS=LZW"
)

#####################################################################################
############################ END ####################################################