# ============================================================
# SECTION 0: PACKAGES AND LIBRARIES
# ============================================================

library(dplyr)
library(terra)

# ============================================================
# SECTION 1: MOSQUITO OCCURRENCE DATA PREPARATION
# ============================================================

# ------------------------------------------------------------
# 1.1 Identify available occurrence datasets
# ------------------------------------------------------------

presence_dir <- file.path(
  "Data", "Italy", "Presence_culex"
)

stopifnot(dir.exists(presence_dir))

presence_files <- list.files(
  presence_dir,
  full.names = TRUE,
  recursive = TRUE
)

print(basename(presence_files))

# ============================================================
# SECTION 2: IMPORT AND INSPECT MOSQUITO SURVEILLANCE DATA
# ============================================================

# Import surveillance data
abundance_file <- file.path(
  presence_dir,
  "longo_et_al_db_aggregated.csv"
)

abundance_data <- readr::read_csv(
  abundance_file,
  show_col_types = FALSE
)

# Dataset dimensions
dim(abundance_data)

# Variable names
names(abundance_data)

# Preview
head(abundance_data)

# Check structure
glimpse(abundance_data)

# ============================================================
# SECTION 3: CHECK SURVEILLANCE DATA
# ============================================================

# 1. Species represented in the dataset
abundance_data %>%
  count(genus, species, sort = TRUE)

# 2. Surveillance characteristics
abundance_data %>%
  count(trap_type, life_stage, sex, sort = TRUE)

# 3. Abundance values
abundance_data %>%
  summarise(
    total_records = n(),
    positive_records = sum(value > 0, na.rm = TRUE),
    zero_records = sum(value == 0, na.rm = TRUE),
    missing_values = sum(is.na(value)),
    unique_cells = n_distinct(cell_id),
    first_year = min(year, na.rm = TRUE),
    last_year = max(year, na.rm = TRUE)
  )

# 4. Number of unique cells with positive Culex pipiens records
abundance_data %>%
  filter(
    genus == "Culex",
    species == "pipiens",
    value > 0
  ) %>%
  summarise(
    positive_records = n(),
    unique_positive_cells = n_distinct(cell_id)
  )

# ============================================================
# SECTION 4: CHECK SURVEILLANCE GRID COORDINATES
# ============================================================


# Check coordinate consistency within each cell
coordinate_check <- abundance_data %>%
  group_by(cell_id) %>%
  summarise(
    n_longitudes = n_distinct(longitude),
    n_latitudes = n_distinct(latitude),
    .groups = "drop"
  )

# Summary
coordinate_check %>%
  summarise(
    total_cells = n(),
    cells_multiple_longitudes = sum(n_longitudes > 1),
    cells_multiple_latitudes = sum(n_latitudes > 1)
  )

# Check missing coordinates and cell IDs
abundance_data %>%
  summarise(
    missing_cell_ids = sum(is.na(cell_id)),
    missing_longitude = sum(is.na(longitude)),
    missing_latitude = sum(is.na(latitude))
  )

# ============================================================
# SECTION 5: PREPARE GRID-LEVEL OCCURRENCE DATA
# ============================================================

culex_cells_prepared <- abundance_data %>%
  group_by(cell_id) %>%
  summarise(
    longitude = first(longitude),
    latitude = first(latitude),
    n_records = n(),
    n_positive = sum(value > 0, na.rm = TRUE),
    n_zero = sum(value == 0, na.rm = TRUE),
    first_year = min(year, na.rm = TRUE),
    last_year = max(year, na.rm = TRUE),
    presence = as.integer(any(value > 0, na.rm = TRUE)),
    .groups = "drop"
  )

# Inspect
glimpse(culex_cells_prepared)

# Count surveyed cells by presence status
culex_cells_prepared %>%
  count(presence)

# ============================================================
# SECTION 6: EXTRACT ENVIRONMENTAL VALUES
# ============================================================

# Extract predictors at surveillance coordinates
env_at_cells <- terra::extract(
  environmental_stack_italy_aligned,
  culex_cells_prepared[, c("longitude", "latitude")]
)

# Add environmental completeness indicator
culex_cells_prepared$complete_env <- complete.cases(
  env_at_cells[, -1]
)

# Summarise
culex_cells_prepared %>%
  count(presence, complete_env)

# Count positive cells with complete environmental data
culex_cells_prepared %>%
  filter(presence == 1, complete_env) %>%
  summarise(n_presence_cells = n())

# ============================================================
# SECTION 7: CHECK SPATIAL DUPLICATION
# ============================================================


# Select presence locations with complete environmental data
culex_presences <- culex_cells_prepared %>%
  filter(
    presence == 1,
    complete_env
  )

# Identify environmental raster cell for each presence
culex_presences$raster_cell <- terra::cellFromXY(
  environmental_stack_italy_aligned[[1]],
  as.matrix(
    culex_presences[, c("longitude", "latitude")]
  )
)

# Summarise spatial duplication
culex_presences %>%
  summarise(
    n_presence_locations = n(),
    n_unique_raster_cells = n_distinct(raster_cell),
    n_duplicate_locations = n() - n_distinct(raster_cell)
  )

# ============================================================
# SECTION 8: PREPARE FINAL MAXENT PRESENCE DATA
# ============================================================

# Retain presences with complete environmental data
culex_presences_italy <- culex_presences %>%
  distinct(raster_cell, .keep_all = TRUE) %>%
  select(
    cell_id,
    longitude,
    latitude,
    raster_cell
  )

# Verify final presence count
cat(
  "Final number of MaxEnt presence locations:",
  nrow(culex_presences_italy),
  "\n"
)

# Verify spatial uniqueness
stopifnot(
  nrow(culex_presences_italy) == 626,
  !anyDuplicated(culex_presences_italy$raster_cell)
)

head(culex_presences_italy)

# ============================================================
# SECTION 9: EXPORT FINAL CULEX PIPIENS PRESENCE DATA
# ============================================================

# Prepare final presence dataset
culex_presences_italy_final <- culex_presences_italy %>%
  mutate(
    species = "Culex pipiens"
  ) %>%
  select(
    species,
    longitude,
    latitude,
    cell_id
  )

# Verify dataset
print(head(culex_presences_italy_final))

# Export CSV
readr::write_csv(
  culex_presences_italy_final,
  file.path(
    presence_dir,
    "Culex_pipiens_Presence_Italy_Final.csv"
  )
)
#----------------------------------------------------
# Reload exported CSV
culex_presences_italy_final_reload <- readr::read_csv(
  file.path(
    presence_dir,
    "Culex_pipiens_Presence_Italy_Final.csv"
  )
)
names(culex_presences_italy_final_reload)
# Verify
stopifnot(
  nrow(culex_presences_italy_final_reload) == 626,
  all(culex_presences_italy_final_reload$species == "Culex pipiens"),
  !anyNA(culex_presences_italy_final_reload)
)


print(head(culex_presences_italy_final_reload))

#####################################################################################
############################ END ####################################################