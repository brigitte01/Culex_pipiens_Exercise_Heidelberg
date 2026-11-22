
# ============================================================
# SECTION 0: PACKAGES AND LIBRARIES
# ============================================================


library(terra)
library(maxnet)
library(terra)
library(dplyr)

# ------------------------------------------------------------
# Set directory
# ------------------------------------------------------------

covariates_dir <- "Data/Europe/Covariates"

# ============================================================
# SECTION 1: LOAD EUROPEAN ENVIRONMENTAL COVARIATES
# ============================================================

# ------------------------------------------------------------
# 1.1 Load the five predictors
# ------------------------------------------------------------

builtup_europe <- rast(
  file.path(covariates_dir,
            "Built_up_2024_CERRA.tif")
)

cropland_europe <- rast(
  file.path(covariates_dir,
            "Cropland_2024_CERRA.tif")
)


temperature_europe <- rast(
  file.path(covariates_dir,
            "cerra_mean_2024-2025_2m_temperature.tif")
)

humidity_europe <- rast(
  file.path(covariates_dir,
            "cerra_mean_2024-2025_2m_relative_humidity.tif")
)

precipitation_europe <- rast(
  file.path(covariates_dir,
            "cerra_mean_2024-2025_total_precipitation.tif")
)
# ============================================================
# SECTION 2: INSPECT RASTER PROPERTIES
# ============================================================

europe_rasters <- list(
  Built_up = builtup_europe,
  Cropland = cropland_europe,
  Temperature = temperature_europe,
  Humidity = humidity_europe,
  Precipitation = precipitation_europe
)

for (variable in names(europe_rasters)) {
  
  cat("\n====================================\n")
  cat("Variable:", variable, "\n")
  cat("====================================\n")
  
  print(europe_rasters[[variable]])
  
  cat(
    "Missing cells:",
    terra::global(
      is.na(europe_rasters[[variable]]),
      "sum"
    )[1, 1],
    "\n"
  )
}
# ============================================================
# SECTION 3: CHECK SPATIAL ALIGNMENT
# ============================================================

for (variable in names(europe_rasters)[-1]) {
  
  cat(
    variable,
    "matches Built-up geometry:",
    terra::compareGeom(
      builtup_europe,
      europe_rasters[[variable]],
      stopOnError = FALSE
    ),
    "\n"
  )
}

# ============================================================
# SECTION 4: VERIFY GRID ALIGNMENT
# ============================================================

# Check climate rasters against one another
terra::compareGeom(
  temperature_europe,
  humidity_europe,
  precipitation_europe,
  stopOnError = FALSE
)

# Check that the climate grid is aligned with CERRA
terra::res(temperature_europe)
terra::res(builtup_europe)

terra::crs(
  temperature_europe
) == terra::crs(
  builtup_europe
)

# Check raster grid origin
terra::origin(temperature_europe)
terra::origin(builtup_europe)

# Check whether climate extent lies within land-cover extent
terra::xmin(temperature_europe) >= terra::xmin(builtup_europe) &&
  terra::xmax(temperature_europe) <= terra::xmax(builtup_europe) &&
  terra::ymin(temperature_europe) >= terra::ymin(builtup_europe) &&
  terra::ymax(temperature_europe) <= terra::ymax(builtup_europe)

# ============================================================
# SECTION 5: CROP LAND COVER TO CLIMATE EXTENT
# ============================================================

builtup_europe_aligned <- terra::crop(
  builtup_europe,
  temperature_europe,
  snap = "near"
)

cropland_europe_aligned <- terra::crop(
  cropland_europe,
  temperature_europe,
  snap = "near"
)

# Verify geometry
cat(
  "Built-up matches climate:",
  terra::compareGeom(
    builtup_europe_aligned,
    temperature_europe,
    stopOnError = FALSE
  ),
  "\n"
)

cat(
  "Cropland matches climate:",
  terra::compareGeom(
    cropland_europe_aligned,
    temperature_europe,
    stopOnError = FALSE
  ),
  "\n"
)

# ============================================================
# SECTION 6: BUILD EUROPEAN PREDICTOR STACK
# ============================================================

europe_stack <- c(
  builtup_europe_aligned,
  cropland_europe_aligned,
  temperature_europe,
  humidity_europe,
  precipitation_europe
)

names(europe_stack) <- c(
  "Built_up_mean",
  "Cropland_mean",
  "Temperature",
  "Humidity",
  "Precipitation"
)

print(europe_stack)
plot(europe_stack)

# ============================================================
# SECTION 7: SUMMARY OF EUROPEAN PREDICTORS
# ============================================================

# Minimum, mean and maximum for each predictor
terra::global(
  europe_stack,
  fun = c("min", "mean", "max"),
  na.rm = TRUE
)

# Number of missing cells for each predictor
missing_summary <- data.frame(
  Variable = names(europe_stack),
  
  Missing_cells = sapply(
    1:terra::nlyr(europe_stack),
    function(i) {
      terra::global(
        is.na(europe_stack[[i]]),
        "sum"
      )[1, 1]
    }
  )
)

missing_summary$Missing_percent <- round(
  100 * missing_summary$Missing_cells /
    terra::ncell(europe_stack),
  2
)

print(missing_summary)

# ============================================================
# SECTION 8: COMPLETE-DATA MASK
# ============================================================

# Count available predictors per cell
n_available <- terra::app(
  !is.na(europe_stack),
  fun = "sum"
)

# Cells with all five predictors
complete_mask <- terra::ifel(
  n_available == 5,
  1,
  NA
)

names(complete_mask) <- "Complete_data"

# Number of complete cells
n_complete <- terra::global(
  !is.na(complete_mask),
  "sum",
  na.rm = TRUE
)[1, 1]

cat(
  "Total grid cells:",
  terra::ncell(europe_stack),
  "\n"
)

cat(
  "Cells with all five predictors:",
  n_complete,
  "\n"
)

cat(
  "Percentage of grid with complete data:",
  round(
    100 * n_complete / terra::ncell(europe_stack),
    2
  ),
  "%\n"
)

# ============================================================
# SECTION 9: VISUALIZE COMPLETE-DATA COVERAGE
# ============================================================

plot(
  complete_mask,
  col = "darkgreen",
  legend = FALSE,
  main = "European Cells with All Five Predictors"
)

# ============================================================
# SECTION 10: EUROPEAN ENVIRONMENTAL RANGES
# ============================================================

# Keep only cells with all five predictors
europe_stack_complete <- terra::mask(
  europe_stack,
  complete_mask
)

# Summary statistics
europe_ranges <- terra::global(
  europe_stack_complete,
  fun = c("min", "mean", "max"),
  na.rm = TRUE
)

print(europe_ranges)

# ============================================================
# SECTION 11: ITALIAN ENVIRONMENTAL RANGES
# ============================================================
# Load Italian environmental predictors

environmental_stack_italy_aligned <- terra::rast(
  file.path(
    "Data",
    "Italy",
    "Covariates",
    "Environmental_Predictors_Italy.tif"
  )
)
# Use the same stack for the Italian range summaries
environmental_stack_italy <- environmental_stack_italy_aligned


italy_ranges <- terra::global(
  environmental_stack_italy,
  fun = c("min", "mean", "max"),
  na.rm = TRUE
)

print(italy_ranges)

# ============================================================
# SECTION 12: ENVIRONMENTAL EXTRAPOLATION CHECK
# ============================================================

# Italian environmental ranges
italy_ranges <- terra::global(
  environmental_stack_italy_aligned,
  c("min", "max"),
  na.rm = TRUE
)

# Initialize a raster for counting out-of-range predictors
novelty_count <- europe_stack_complete[[1]]
terra::values(novelty_count) <- 0

# Count predictors outside Italian training ranges
for (i in seq_len(terra::nlyr(europe_stack_complete))) {
  
  predictor <- europe_stack_complete[[i]]
  
  lower <- italy_ranges$min[i]
  upper <- italy_ranges$max[i]
  
  outside_range <- terra::ifel(
    predictor < lower | predictor > upper,
    1,
    0
  )
  
  novelty_count <- novelty_count + outside_range
}

# Name output
names(novelty_count) <- "Number_outside_training_range"

# Visualize
plot(
  novelty_count,
  main = "Environmental Extrapolation: Europe"
)

# ============================================================
# SECTION 13: EXTRAPOLATION BY PREDICTOR
# ============================================================

# Italian training ranges
italy_ranges <- terra::global(
  environmental_stack_italy_aligned,
  c("min", "max"),
  na.rm = TRUE
)

# Number of complete European cells
n_complete <- terra::global(
  !is.na(europe_stack_complete[[1]]),
  "sum",
  na.rm = TRUE
)[1, 1]

# Calculate extrapolation for each predictor
extrapolation_summary <- lapply(
  seq_len(terra::nlyr(europe_stack_complete)),
  function(i) {
    
    predictor <- europe_stack_complete[[i]]
    
    lower <- italy_ranges$min[i]
    upper <- italy_ranges$max[i]
    
    # Count cells below and above Italian training range
    below <- terra::global(
      predictor < lower,
      "sum",
      na.rm = TRUE
    )[1, 1]
    
    above <- terra::global(
      predictor > upper,
      "sum",
      na.rm = TRUE
    )[1, 1]
    
    data.frame(
      Variable = names(europe_stack_complete)[i],
      Italy_min = lower,
      Italy_max = upper,
      Cells_below = below,
      Cells_above = above,
      Total_outside = below + above,
      Percent_outside = round(
        100 * (below + above) / n_complete,
        2
      )
    )
  }
)

extrapolation_summary <- bind_rows(
  extrapolation_summary
)

print(extrapolation_summary)

# ============================================================
# SECTION 14: EUROPEAN MAXENT PREDICTION
# ============================================================
# ============================================================
# SECTION 14: EUROPEAN MAXENT PREDICTION
# ============================================================

# ------------------------------------------------------------
# 14.0 Load the selected Italian MaxEnt model
# ------------------------------------------------------------

model_path <- file.path(
  "Outputs",
  "Italy",
  "Maxent",
  "MaxEnt_Best_Model_Italy.rds"
)

# Check that the saved model exists
stopifnot(file.exists(model_path))

# Load the selected model
best_model_italy <- readRDS(model_path)

# Verify that the model was loaded
stopifnot(inherits(best_model_italy, "maxnet"))

cat("Selected Italian MaxEnt model loaded successfully.\n")
# ------------------------------------------------------------
# 14.1 Verify predictor compatibility
# ------------------------------------------------------------

# Predictors used to train the Italian model
italy_predictors <- names(
  environmental_stack_italy_aligned
)

# Predictors available for Europe
europe_predictors <- names(
  europe_stack_complete
)

# Check that predictor names match
stopifnot(
  identical(
    italy_predictors,
    europe_predictors
  )
)

cat("Predictor names match successfully.\n")

# ------------------------------------------------------------
# 14.2 Define MaxEnt prediction function
# ------------------------------------------------------------

predict_maxnet_europe <- function(model, data, ...) {
  
  predict(
    model,
    newdata = as.data.frame(data),
    type = "cloglog",
    clamp = TRUE
  )
  
}

# ------------------------------------------------------------
# 14.3 Project selected Italian model to Europe
# ------------------------------------------------------------

maxent_europe <- terra::predict(
  europe_stack_complete,
  best_model_italy,
  fun = predict_maxnet_europe,
  na.rm = TRUE
)

# Rename output
names(maxent_europe) <- "Culex_pipiens_suitability"

# ------------------------------------------------------------
# 14.4 Inspect predictions
# ------------------------------------------------------------

print(maxent_europe)

terra::global(
  maxent_europe,
  c("min", "mean", "max"),
  na.rm = TRUE
)

# Verify prediction range
prediction_range <- terra::global(
  maxent_europe,
  c("min", "max"),
  na.rm = TRUE
)

stopifnot(
  prediction_range$min >= 0,
  prediction_range$max <= 1
)

# ============================================================
# SECTION 15: EUROPEAN SUITABILITY MAP WITH SEPARATE LEGEND
# ============================================================

# Suitability categories
breaks <- c(0, 0.25, 0.50, 0.75, 1.000001)

# Colours
suitability_colors <- c(
  "#2166AC",  # Very low
  "#FEE08B",  # Low
  "#F46D43",  # Moderate
  "#D73027"   # High
)

# Create separate panels for map and legend
layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1.5))

# ---------------------------
# Panel 1: Suitability map
# ---------------------------

par(mar = c(2, 2, 4, 1))

plot(
  maxent_europe,
  col = suitability_colors,
  breaks = breaks,
  main = expression(
    italic("Culex pipiens") ~ "Habitat Suitability - Europe"
  ),
  axes = FALSE,
  legend = FALSE
)

# ---------------------------
# Panel 2: Legend
# ---------------------------

par(mar = c(2, 0, 4, 1))

plot.new()

legend(
  "center",
  legend = c(
    "Very low (0–0.25)",
    "Low (0.25–0.50)",
    "Moderate (0.50–0.75)",
    "High (0.75–1.00)"
  ),
  fill = suitability_colors,
  title = "Habitat suitability",
  bty = "n",
  cex = 0.9
)

# Restore normal plotting layout
layout(1)

# ============================================================
# SECTION 16: EXPORT EUROPEAN SUITABILITY
# ============================================================

prediction_dir <- "Outputs"

dir.create(
  prediction_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

terra::writeRaster(
  maxent_europe,
  filename = file.path(
    prediction_dir,
    "Culex_pipiens_MaxEnt_Europe_CERRA_2024_2025.tif"
  ),
  overwrite = TRUE,
  datatype = "FLT4S",
  gdal = c("COMPRESS=LZW")
)

cat("European MaxEnt suitability raster saved.\n")

##############################################################################
########################### END
