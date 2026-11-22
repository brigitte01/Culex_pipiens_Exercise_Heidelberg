# ============================================================
# SECTION 0: PACKAGES AND LIBRARIES
# ============================================================
library(terra)
library(ENMeval)
library(maxnet)
library(dplyr)
library(ggplot2)
library(tidyr)

# ============================================================
# SECTION 1: MAXENT MODEL TUNING — ITALY
# ============================================================

set.seed(123)

# ------------------------------------------------------------
# 1.1 Prepare modelling inputs
# ------------------------------------------------------------

# Environmental predictors
env <- environmental_stack_italy_aligned

# Presence coordinates
occ_coords <- culex_presences_italy_final %>%
  select(longitude, latitude) %>%
  as.data.frame()

# Background coordinates
bg_coords <- background_italy %>%
  select(longitude, latitude) %>%
  as.data.frame()

# ------------------------------------------------------------
# Verify modelling inputs
# ------------------------------------------------------------

cat("Presence locations:", nrow(occ_coords), "\n")
cat("Background locations:", nrow(bg_coords), "\n")
cat("Environmental predictors:", terra::nlyr(env), "\n")

print(names(env))

# Check missing coordinates
stopifnot(
  nrow(occ_coords) == 626,
  nrow(bg_coords) == 2000,
  !anyNA(occ_coords),
  !anyNA(bg_coords),
  terra::nlyr(env) == 5
)

# Verify environmental values at modelling locations
occ_env <- terra::extract(env, occ_coords)
bg_env  <- terra::extract(env, bg_coords)

stopifnot(
  all(complete.cases(occ_env)),
  all(complete.cases(bg_env))
)

# ------------------------------------------------------------
# 1.2 Define MaxEnt tuning parameters
# ------------------------------------------------------------

tune_args_italy <- list(
  fc = c("L", "Q", "LQ", "LQH", "LQHP"),
  rm = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)
)

# Number of candidate models
n_models <- length(tune_args_italy$fc) *
  length(tune_args_italy$rm)

cat("Total candidate MaxEnt models:", n_models, "\n")

# ------------------------------------------------------------
# 1.3 Tune MaxEnt models using ENMeval
# ------------------------------------------------------------

maxent_tuning_italy <- ENMeval::ENMevaluate(
  occs = occ_coords,
  envs = env,
  bg = bg_coords,
  algorithm = "maxnet",
  partitions = "block",
  tune.args = tune_args_italy,
  raster.preds = FALSE,
  parallel = FALSE
)

# ============================================================
# SECTION 2: SELECT BEST MAXENT MODEL
# ============================================================


# Extract tuning results
tuning_results_italy <- ENMeval::eval.results(
  maxent_tuning_italy
)

# ------------------------------------------------------------
# 2.1 Retain models with delta AICc <= 10
# ------------------------------------------------------------

candidate_models_italy <- tuning_results_italy %>%
  filter(
    !is.na(delta.AICc),
    delta.AICc <= 10,
    !is.na(or.10p.avg)
  ) %>%
  arrange(or.10p.avg, delta.AICc)

# Display candidate models
print(candidate_models_italy)

cat(
  "Number of candidate models:",
  nrow(candidate_models_italy),
  "\n"
)

# ------------------------------------------------------------
# 2.2 Select model with lowest omission rate
# ------------------------------------------------------------

stopifnot(nrow(candidate_models_italy) > 0)

best_model_result_italy <- candidate_models_italy %>%
  slice(1)

print(best_model_result_italy)

# Selected model summary
best_model_result_italy %>%
  select(
    fc,
    rm,
    AICc,
    delta.AICc,
    auc.val.avg,
    or.10p.avg,
    cbi.val.avg
  )

# ============================================================
# SECTION 3: EXTRACT SELECTED MAXENT MODEL
# ============================================================

# Selected model identifier
best_tune_id_italy <- as.character(
  best_model_result_italy$tune.args[1]
)

cat("Selected model:", best_tune_id_italy, "\n")

# Extract fitted model
best_model_italy <- ENMeval::eval.models(
  maxent_tuning_italy
)[[best_tune_id_italy]]

# Inspect model
print(best_model_italy)

# Confirm extraction
stopifnot(!is.null(best_model_italy))

# ============================================================
# SECTION 4: SAVE MODEL AND TUNING RESULTS
# ============================================================

model_dir <- file.path("Outputs", "Italy", "MaxEnt")

# Save fitted model
saveRDS(
  best_model_italy,
  file.path(model_dir, "MaxEnt_Best_Model_Italy.rds")
)

# Save full tuning object
saveRDS(
  maxent_tuning_italy,
  file.path(model_dir, "MaxEnt_Tuning_Italy.rds")
)

# Save evaluation table
readr::write_csv(
  tuning_results_italy,
  file.path(model_dir, "MaxEnt_Tuning_Results_Italy.csv")
)

# Save selected configuration
readr::write_csv(
  best_model_result_italy,
  file.path(model_dir, "MaxEnt_Selected_Model_Italy.csv")
)

# ============================================================
# SECTION 5: FINAL MODEL EVALUATION SUMMARY
# ============================================================

# Extract evaluation metrics for the selected model
model_evaluation_italy <- tuning_results_italy %>%
  filter(tune.args == best_tune_id_italy) %>%
  select(
    fc,
    rm,
    auc.train,
    auc.val.avg,
    auc.val.sd,
    auc.diff.avg,
    auc.diff.sd,
    cbi.train,
    cbi.val.avg,
    cbi.val.sd,
    or.10p.avg,
    or.10p.sd,
    or.mtp.avg,
    or.mtp.sd,
    AICc,
    delta.AICc,
    ncoef
  )

# Display results
print(model_evaluation_italy)

# Save evaluation summary
readr::write_csv(
  model_evaluation_italy,
  file.path(
    model_dir,
    "MaxEnt_Final_Model_Evaluation_Italy.csv"
  )
)

# ============================================================
# SECTION 6: INSPECT SPATIAL CROSS-VALIDATION RESULTS
# ============================================================

# Inspect the ENMevaluate object
class(maxent_tuning_italy)

# List available components
slotNames(maxent_tuning_italy)

# Check available partition information
ENMeval::eval.partition.settings(
  maxent_tuning_italy
)

# Check partition assignments
partition_groups_italy <- ENMeval::eval.occs.grp(
  maxent_tuning_italy
)

print(table(partition_groups_italy))

# ============================================================
# SECTION 7: EVALUATION BY SPATIAL FOLD
# ============================================================

# Extract partition-level evaluation results
partition_results_italy <- ENMeval::eval.results.partitions(
  maxent_tuning_italy
)

# Inspect structure
dim(partition_results_italy)

# Check available columns
names(partition_results_italy)

# Display first rows
head(partition_results_italy)

# Inspect model identifier columns
names(partition_results_italy)[
  grepl(
    "tune|fc|rm|partition|fold",
    names(partition_results_italy),
    ignore.case = TRUE
  )
]

# ============================================================
# SECTION 8: SPATIAL CROSS-VALIDATION PERFORMANCE
# ============================================================

# Extract results for the selected model
best_model_folds_italy <- partition_results_italy %>%
  filter(tune.args == best_tune_id_italy) %>%
  arrange(fold)

# Display results
print(best_model_folds_italy)

# Confirm that all four folds are present
stopifnot(nrow(best_model_folds_italy) == 4)

# ------------------------------------------------------------
# Prepare data for plotting
# ------------------------------------------------------------

fold_plot_data <- best_model_folds_italy %>%
  select(
    fold,
    auc.val,
    cbi.val,
    or.10p
  ) %>%
  pivot_longer(
    cols = c(auc.val, cbi.val, or.10p),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    fold = factor(fold),
    Metric = recode(
      Metric,
      "auc.val" = "Validation AUC",
      "cbi.val" = "Validation Boyce index",
      "or.10p" = "10th-percentile omission rate"
    )
  )

# ------------------------------------------------------------
# Plot fold-level evaluation metrics
# ------------------------------------------------------------

fold_evaluation_plot <- ggplot(
  fold_plot_data,
  aes(x = fold, y = Value)
) +
  geom_col(
    fill = "steelblue",
    width = 0.65
  ) +
  geom_text(
    aes(label = sprintf("%.3f", Value)),
    vjust = -0.4,
    size = 3.5
  ) +
  facet_wrap(
    ~ Metric,
    scales = "free_y",
    ncol = 3
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "Spatial cross-validation performance",
    subtitle = "Culex pipiens MaxEnt — Italy (LQHP, RM = 2)",
    x = "Spatial validation fold",
    y = "Metric value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(fold_evaluation_plot)

# Save fold-level results
readr::write_csv(
  best_model_folds_italy,
  file.path(
    model_dir,
    "MaxEnt_Selected_Model_Fold_Evaluation_Italy.csv"
  )
)

# Save evaluation plot
ggsave(
  filename = file.path(
    model_dir,
    "MaxEnt_Spatial_CrossValidation_Italy.png"
  ),
  plot = fold_evaluation_plot,
  width = 11,
  height = 4.5,
  dpi = 300
)

# ============================================================
# SECTION 9: JACKKNIFE DELTA VALIDATION BOYCE INDEX
# ============================================================

# Full-model validation Boyce index
full_model_cbi <- best_model_result_italy$cbi.val.avg[1]

cat("Full-model validation Boyce:", full_model_cbi, "\n")

# Extract the original spatial partition assignments
occ_groups <- ENMeval::eval.occs.grp(maxent_tuning_italy)
bg_groups  <- ENMeval::eval.bg.grp(maxent_tuning_italy)

# Check group assignments
table(occ_groups)
table(bg_groups)

# Environmental predictors
predictor_names <- names(environmental_stack_italy_aligned)

print(predictor_names)

# ------------------------------------------------------------
# 9.1 Leave-one-variable-out jackknife
# ------------------------------------------------------------

jackknife_results <- list()

for (var in predictor_names) {
  
  cat("\n====================================\n")
  cat("Omitting predictor:", var, "\n")
  cat("====================================\n")
  
  # Keep all predictors except the omitted one
  vars_remaining <- setdiff(predictor_names, var)
  
  env_reduced <- environmental_stack_italy_aligned[[
    vars_remaining
  ]]
  
  # Fit model using the same spatial partitions
  jk_eval <- ENMeval::ENMevaluate(
    occs = occ_coords,
    envs = env_reduced,
    bg = bg_coords,
    algorithm = "maxnet",
    partitions = "user",
    user.grp = list(
      occs.grp = occ_groups,
      bg.grp = bg_groups
    ),
    tune.args = list(
      fc = "LQHP",
      rm = 2
    ),
    raster.preds = FALSE,
    parallel = FALSE
  )
  
  # Extract evaluation results
  jk_results <- ENMeval::eval.results(jk_eval)
  
  # Save summary
  jackknife_results[[var]] <- data.frame(
    Variable = var,
    CBI_without = jk_results$cbi.val.avg[1],
    CBI_without_sd = jk_results$cbi.val.sd[1],
    Delta_CBI = full_model_cbi - jk_results$cbi.val.avg[1],
    AUC_without = jk_results$auc.val.avg[1],
    Omission_without = jk_results$or.10p.avg[1]
  )
  
}

# Combine results
jackknife_summary <- bind_rows(jackknife_results)

# Rank predictors by delta Boyce
jackknife_summary <- jackknife_summary %>%
  arrange(desc(Delta_CBI))

print(jackknife_summary)

# ------------------------------------------------------------
# 9.2 Jackknife Delta Boyce Index Plot
# ------------------------------------------------------------

jackknife_plot <- ggplot(
  jackknife_summary,
  aes(
    x = reorder(Variable, Delta_CBI),
    y = Delta_CBI,
    fill = Delta_CBI > 0
  )
) +
  geom_col(width = 0.7) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "TRUE" = "steelblue",
      "FALSE" = "tomato"
    ),
    guide = "none"
  ) +
  labs(
    title = "Jackknife variable importance",
    subtitle = "Change in validation Boyce index when each predictor is omitted",
    x = NULL,
    y = expression(Delta ~ "Validation Boyce index")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

print(jackknife_plot)

# ============================================================
# SECTION 10: ENVIRONMENTAL RESPONSE CURVES
# ============================================================

# Environmental values from complete raster cells
env_df <- as.data.frame(
  terra::values(
    environmental_stack_italy_aligned,
    na.rm = TRUE
  )
)

# Mean values for all predictors
env_means <- colMeans(env_df)

print(env_means)

# ------------------------------------------------------------
# 10.1 Generate marginal response curves
# ------------------------------------------------------------

response_curves_italy <- lapply(
  names(env_df),
  function(var) {
    
    # Sequence covering the observed range
    x_values <- seq(
      min(env_df[[var]], na.rm = TRUE),
      max(env_df[[var]], na.rm = TRUE),
      length.out = 200
    )
    
    # Reference dataset: other predictors at their means
    newdata <- as.data.frame(
      matrix(
        rep(env_means, each = length(x_values)),
        nrow = length(x_values)
      )
    )
    
    names(newdata) <- names(env_means)
    
    # Vary the predictor of interest
    newdata[[var]] <- x_values
    
    # Predict using the selected MaxEnt model
    predictions <- predict(
      best_model_italy,
      newdata = newdata,
      type = "cloglog",
      clamp = TRUE
    )
    
    # Return response-curve data
    data.frame(
      Variable = var,
      Value = x_values,
      Suitability = as.numeric(predictions)
    )
  }
)

# Combine all curves
response_curves_italy <- bind_rows(response_curves_italy)

head(response_curves_italy)

# ------------------------------------------------------------
# 10.2 Plot environmental response curves
# ------------------------------------------------------------

response_plot_italy <- ggplot(
  response_curves_italy,
  aes(x = Value, y = Suitability)
) +
  geom_line(
    color = "steelblue",
    linewidth = 1.1
  ) +
  facet_wrap(
    ~ Variable,
    scales = "free_x",
    ncol = 2
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  labs(
    title = "Environmental response curves",
    subtitle = "Culex pipiens MaxEnt — Italy (LQHP, RM = 2)",
    x = "Environmental predictor value",
    y = "Predicted suitability (cloglog)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(response_plot_italy)

# ------------------------------------------------------------
# 10.3 Export response curves
# ------------------------------------------------------------

# Save underlying data
readr::write_csv(
  response_curves_italy,
  file.path(
    model_dir,
    "MaxEnt_Response_Curves_Italy.csv"
  )
)

# Save figure
ggsave(
  filename = file.path(
    model_dir,
    "MaxEnt_Response_Curves_Italy.png"
  ),
  plot = response_plot_italy,
  width = 10,
  height = 9,
  dpi = 300
)

##############################################################################
########################### END

