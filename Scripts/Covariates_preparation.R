# ============================================================
# SECTION 0: PACKAGES AND LIBRARIES
# ============================================================
library(dplyr)
library(terra)
library(corrplot)

# ============================================================
# SECTION 1: IDENTIFY ITALIAN ENVIRONMENTAL INPUT FILES
# ============================================================

# Relative to the Italy_SDM_Analysis project directory
covariate_dir <- file.path(
  "Data", "Italy", "Covariates"
)

# Check directory exists
stopifnot(dir.exists(covariate_dir))

# List files, including extensions
covariate_files <- list.files(
  covariate_dir,
  full.names = TRUE
)

print(basename(covariate_files))

# ============================================================
# SECTION 2: ITALIAN ENVIRONMENTAL PREDICTORS
# ============================================================

# ------------------------------------------------------------
# 2.1 Load prepared environmental datasets
# ------------------------------------------------------------

covariate_dir <- file.path("Data", "Italy", "Covariates")

# Land-cover proportions (2008–2020)
landcover_italy <- terra::rast(
  file.path(covariate_dir, "Landcover_Mean_2008_2020.tif")
)

# ERA5-Land mean temperature (2014–2022)
temp_italy <- terra::rast(
  file.path(covariate_dir,
            "era5land_mean_2014-2022_2m_temperature.tif")
)

# ERA5-Land mean relative humidity (2014–2022)
humidity_italy <- terra::rast(
  file.path(covariate_dir,
            "era5land_mean_2014-2022_2m_relative_humidity.tif")
)

# ERA5-Land mean annual precipitation (2014–2022)
precip_italy <- terra::rast(
  file.path(covariate_dir,
            "era5land_mean_2014-2022_total_precipitation.tif")
)

# ------------------------------------------------------------
# 2.2 Inspect environmental raster properties
# ------------------------------------------------------------

# Land cover
print(landcover_italy)
names(landcover_italy)

# Climate variables
print(temp_italy)
print(humidity_italy)
print(precip_italy)

# ------------------------------------------------------------
# 2.3 Check environmental variable ranges
# ------------------------------------------------------------

terra::global(
  temp_italy,
  c("min", "mean", "max"),
  na.rm = TRUE
)

terra::global(
  humidity_italy,
  c("min", "mean", "max"),
  na.rm = TRUE
)

terra::global(
  precip_italy,
  c("min", "mean", "max"),
  na.rm = TRUE
)

# ============================================================
# SECTION 3: HARMONIZE ITALIAN ENVIRONMENTAL PREDICTORS
# ============================================================

# ------------------------------------------------------------
# 3.1 Select the two land-cover predictors
# ------------------------------------------------------------

landcover_selected <- landcover_italy[[
  c("Built_up_mean", "Cropland_mean")
]]

# Use land cover as the reference grid
template_italy <- landcover_selected[[1]]

# ------------------------------------------------------------
# 3.2 Align climate predictors to the reference grid
# ------------------------------------------------------------

temp_aligned_italy <- terra::project(
  temp_italy,
  template_italy,
  method = "bilinear"
)

humidity_aligned_italy <- terra::project(
  humidity_italy,
  template_italy,
  method = "bilinear"
)

precip_aligned_italy <- terra::project(
  precip_italy,
  template_italy,
  method = "bilinear"
)

# Rename climate layers
names(temp_aligned_italy) <- "Temperature"
names(humidity_aligned_italy) <- "Humidity"
names(precip_aligned_italy) <- "Precipitation"

# ------------------------------------------------------------
# 3.3 Combine all five predictors
# ------------------------------------------------------------

environmental_stack_italy_aligned <- c(
  landcover_selected,
  temp_aligned_italy,
  humidity_aligned_italy,
  precip_aligned_italy
)

# Inspect
print(environmental_stack_italy_aligned)

names(environmental_stack_italy_aligned)

# ============================================================
# SECTION 4: CHECK ENVIRONMENTAL RASTER ALIGNMENT
# ============================================================

geometry_check <- sapply(
  2:terra::nlyr(environmental_stack_italy_aligned),
  function(i) {
    terra::compareGeom(
      environmental_stack_italy_aligned[[1]],
      environmental_stack_italy_aligned[[i]],
      stopOnError = FALSE
    )
  }
)

print(geometry_check)

cat(
  "All environmental layers aligned:",
  all(geometry_check),
  "\n"
)

# ============================================================
# SECTION 5: SAVE VALIDATED ENVIRONMENTAL PREDICTORS
# ============================================================

# Ensure correct layer names
names(environmental_stack_italy_aligned) <- c(
  "Built_up_mean",
  "Cropland_mean",
  "Temperature",
  "Humidity",
  "Precipitation"
)

# Save harmonized predictors
terra::writeRaster(
  environmental_stack_italy_aligned,
  filename = file.path(
    covariate_dir,
    "Environmental_Predictors_Italy.tif"
  ),
  overwrite = TRUE,
  datatype = "FLT4S",
  gdal = "COMPRESS=LZW"
)

# Reload saved environmental stack
environmental_stack_italy <- terra::rast(
  file.path(
    covariate_dir,
    "Environmental_Predictors_Italy.tif"
  )
)

# Inspect metadata
print(environmental_stack_italy)

names(environmental_stack_italy)

terra::varnames(environmental_stack_italy)

# ============================================================
# SECTION 6: CORRELATION ANALYSIS
# ============================================================

# Extract environmental predictor values
env_values <- terra::values(
  environmental_stack_italy_aligned,
  na.rm = TRUE
)

# Convert to data frame
env_df <- as.data.frame(env_values)

# Check dimensions
dim(env_df)

# Inspect
head(env_df)

# ------------------------------------------------------------
# 6.1 Pearson correlation matrix
# ------------------------------------------------------------

cor_matrix <- cor(
  env_df,
  method = "pearson",
  use = "complete.obs"
)

# Display correlation matrix
round(cor_matrix, 3)

# ------------------------------------------------------------
# 6.2 Correlation plot
# ------------------------------------------------------------

corrplot(
  cor_matrix,
  method = "circle",
  type = "upper",
  order = "original",
  addCoef.col = "black",
  number.cex = 0.8,
  tl.col = "black",
  tl.cex = 0.9,
  tl.srt = 45,
  diag = FALSE
)

# ------------------------------------------------------------
# 6.3 Identify highly correlated predictors
# ------------------------------------------------------------

# Extract upper triangle to avoid duplicate pairs
cor_pairs <- which(
  abs(cor_matrix) >= 0.70 &
    upper.tri(cor_matrix),
  arr.ind = TRUE
)

# Create summary table
high_correlations <- data.frame(
  Variable_1 = rownames(cor_matrix)[cor_pairs[, 1]],
  Variable_2 = colnames(cor_matrix)[cor_pairs[, 2]],
  Correlation = cor_matrix[cor_pairs]
)

print(high_correlations)

############################ END ####################################################