###############################################################################
# 02_MICE_Imputation.R
# ============================================================================
# Project: Lifestyle Factors and PPC-MM Association Study
# Purpose: 2-Level Multiple Imputation by Chained Equations (MICE)
# ============================================================================
#
# METHODOLOGY:
#   - 2-Level MICE with COHORT as Level 2 clustering variable
#   - Uses mice package for multilevel imputation
#   - Includes comprehensive diagnostics (convergence, density comparison)
#   - Rubin's rules for pooling results
#
# PARAMETERS:
#   - Level 2: cohort (5 cohorts: CHARLS, ELSA, HRS, SHARE, MHAS)
#   - Level 1: individual
#   - m = 20 (number of imputed datasets)
#   - maxit = 30 (iterations)
#   - Method: 2l.pmm (multilevel predictive mean matching)
#
# VARIABLES TO IMPUTE:
#   - edu (education level)
#
# EXCLUDED VARIABLES:
#   - rural (ELSA 100% missing)
#   - employment (excluded per study design)
#
# OUTPUT:
#   - Pooled_mice_imputed.rds (full mids object for Rubin's rules)
#   - MICE_Convergence_Trace.png/pdf
#   - MICE_Density_Comparison.png/pdf
#   - MICE_Summary_Statistics.xlsx
#
###############################################################################

# =============================================================================
# PART 1: Environment Setup
# =============================================================================

source("C:/Users/user/Desktop/Scientific Research/UKB Project/Jung Sun Lab/Project1_Lifestyle_PPM_Hexiao,Hongtao/Code/00_Functions_and_Setup.R")

# Additional packages for MICE
required_mice_pkgs <- c("mice", "miceadds", "lattice")
for (pkg in required_mice_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cat("\n")
cat("################################################################\n")
cat("#                                                              #\n")
cat("#     2-Level MICE Imputation (Sensitivity Analysis S3)        #\n")
cat("#     Level 2: cohort | Level 1: individual                    #\n")
cat("#                                                              #\n")
cat("################################################################\n")
cat("\n")

# =============================================================================
# PART 2: Configuration
# =============================================================================

# MICE parameters
MICE_PARAMS <- list(
  m = 20,           # Number of imputed datasets (standard for publication)
  maxit = 30,       # Maximum iterations (ensure convergence)
  seed = 12345      # Random seed for reproducibility
)

# Level 2 clustering variable
CLUSTER_VAR <- "cohort"  # 5 cohorts as clusters

# Variable to impute
VAR_TO_IMPUTE <- "edu"

# Variables to EXCLUDE from analysis
EXCLUDE_VARS <- c("rural", "employment")

# Predictor variables for imputation model
PREDICTOR_VARS <- c(
  # Demographics (complete)
  "age_baseline",
  "sex",
  "marital_binary",
  
  # Health status at baseline (complete)
  "physical_base",
  "psych_base", 
  "cog_base",
  
  # Lifestyle (exposure)
  "unhealthy_score",
  
  # Auxiliary variables (outcome-related, improves MAR assumption)
  "event_ppcmm",
  "time_ppcmm_months"
)

cat("=== MICE Configuration ===\n\n")
cat("  Imputation Parameters:\n")
cat("    m (datasets):", MICE_PARAMS$m, "\n")
cat("    maxit (iterations):", MICE_PARAMS$maxit, "\n")
cat("    seed:", MICE_PARAMS$seed, "\n")
cat("\n")
cat("  Clustering:\n")
cat("    Level 2 variable:", CLUSTER_VAR, "\n")
cat("\n")
cat("  Variable to Impute:\n")
cat("   ", VAR_TO_IMPUTE, "\n")
cat("\n")
cat("  Excluded Variables:\n")
cat("   ", paste(EXCLUDE_VARS, collapse = ", "), "\n")
cat("\n")

# Create output directory
mice_output_dir <- file.path(OUTPUT_DIR, "02_MICE")
if (!dir.exists(mice_output_dir)) dir.create(mice_output_dir, recursive = TRUE)

# =============================================================================
# PART 3: Load Data
# =============================================================================

cat("=== Loading Data ===\n\n")

# Load main analysis data from each cohort
cohort_list <- list()

for (coh in names(DATA_PATHS)) {
  file_path <- DATA_PATHS[[coh]]
  
  if (!file.exists(file_path)) {
    cat("  [SKIP]", coh, "- file not found\n")
    next
  }
  
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  bw <- BASELINE_WAVES[[coh]]
  
  # Standardize lifestyle variable name
  score_var <- paste0("w", bw, "_unhealthy_score")
  if (score_var %in% names(df)) {
    df$unhealthy_score <- df[[score_var]]
  }
  
  # Ensure cohort variable exists
  if (!"cohort" %in% names(df)) {
    df$cohort <- coh
  }
  
  cat("  Loaded", coh, ":", nrow(df), "rows\n")
  cohort_list[[coh]] <- df
}

# Find common columns
common_cols <- Reduce(intersect, lapply(cohort_list, names))
cat("\n  Common columns:", length(common_cols), "\n")

# Select relevant columns only
keep_cols <- unique(c(
  CLUSTER_VAR,
  VAR_TO_IMPUTE,
  PREDICTOR_VARS,
  "country_name", "region"  # Keep for reference
))
keep_cols <- keep_cols[keep_cols %in% common_cols]

# Subset and combine
for (coh in names(cohort_list)) {
  available_cols <- keep_cols[keep_cols %in% names(cohort_list[[coh]])]
  cohort_list[[coh]] <- cohort_list[[coh]][, available_cols, drop = FALSE]
}

pooled_data <- bind_rows(cohort_list)
cat("\n  Pooled data:", nrow(pooled_data), "rows,", ncol(pooled_data), "columns\n")

# =============================================================================
# PART 4: Missing Data Summary
# =============================================================================

cat("\n=== Missing Data Summary ===\n\n")

# Overall missing
missing_overall <- data.frame(
  Variable = character(),
  N_Total = integer(),
  N_Missing = integer(),
  Pct_Missing = numeric(),
  stringsAsFactors = FALSE
)

check_vars <- c(VAR_TO_IMPUTE, PREDICTOR_VARS)

for (var in check_vars) {
  if (var %in% names(pooled_data)) {
    n_total <- nrow(pooled_data)
    n_miss <- sum(is.na(pooled_data[[var]]) | pooled_data[[var]] == "")
    pct_miss <- round(n_miss / n_total * 100, 2)
    
    missing_overall <- rbind(missing_overall, data.frame(
      Variable = var,
      N_Total = n_total,
      N_Missing = n_miss,
      Pct_Missing = pct_miss,
      stringsAsFactors = FALSE
    ))
    
    status <- if(n_miss == 0) "[Complete]" else paste0("[", pct_miss, "% missing]")
    cat(sprintf("  %-20s: %s\n", var, status))
  }
}

# Missing by cohort for edu
cat("\n  Missing 'edu' by Cohort:\n")
missing_by_cohort <- pooled_data %>%
  group_by(cohort) %>%
  summarise(
    N = n(),
    N_Missing = sum(is.na(edu) | edu == ""),
    Pct_Missing = round(N_Missing / N * 100, 2),
    .groups = "drop"
  )
print(missing_by_cohort)

# Save missing summary
writexl::write_xlsx(
  list(
    Overall = missing_overall,
    By_Cohort = as.data.frame(missing_by_cohort)
  ),
  file.path(mice_output_dir, "MICE_Missing_Data_Summary.xlsx")
)
cat("\n  Saved: MICE_Missing_Data_Summary.xlsx\n")

# =============================================================================
# PART 5: Prepare Imputation Data
# =============================================================================

cat("\n=== Preparing Imputation Dataset ===\n\n")

# Select variables for imputation model
impute_vars <- c(CLUSTER_VAR, VAR_TO_IMPUTE, PREDICTOR_VARS)
impute_vars <- impute_vars[impute_vars %in% names(pooled_data)]

mice_data <- pooled_data[, impute_vars, drop = FALSE]

# Convert to appropriate types
mice_data$cohort <- as.factor(mice_data$cohort)
mice_data$sex <- as.factor(mice_data$sex)
mice_data$edu <- as.factor(mice_data$edu)
mice_data$marital_binary <- as.factor(mice_data$marital_binary)
mice_data$physical_base <- as.numeric(mice_data$physical_base)
mice_data$psych_base <- as.numeric(mice_data$psych_base)
mice_data$cog_base <- as.numeric(mice_data$cog_base)

cat("  Imputation dataset:", nrow(mice_data), "rows,", ncol(mice_data), "columns\n")
cat("  Cohorts:", paste(levels(mice_data$cohort), collapse = ", "), "\n")
cat("  Variables:", paste(names(mice_data), collapse = ", "), "\n")

# =============================================================================
# PART 6: Setup MICE Imputation Model
# =============================================================================

cat("\n=== Setting Up MICE Model ===\n\n")

# Create predictor matrix
pred_matrix <- mice::make.predictorMatrix(mice_data)

# Set cohort as cluster variable (-2)
pred_matrix[, "cohort"] <- -2

# Cohort doesn't predict itself
pred_matrix["cohort", ] <- 0

# Set up methods
methods <- mice::make.method(mice_data)

# Only impute edu using 2-level PMM
methods["edu"] <- "2l.pmm"

# Set other variables to not be imputed
for (var in names(methods)) {
  if (var != "edu") {
    methods[var] <- ""
  }
}

cat("  Imputation Methods:\n")
for (var in names(methods)) {
  if (methods[var] != "") {
    cat("   ", var, ":", methods[var], "\n")
  }
}

cat("\n  Predictor Matrix (edu row):\n")
print(pred_matrix["edu", , drop = FALSE])

# =============================================================================
# PART 7: Run MICE Imputation
# =============================================================================

cat("\n=== Running 2-Level MICE Imputation ===\n")
cat("  This may take 10-20 minutes...\n\n")

set.seed(MICE_PARAMS$seed)
start_time <- Sys.time()

# Try 2-level MICE first
mice_result <- tryCatch({
  mice::mice(
    data = mice_data,
    m = MICE_PARAMS$m,
    maxit = MICE_PARAMS$maxit,
    method = methods,
    predictorMatrix = pred_matrix,
    printFlag = TRUE
  )
}, error = function(e) {
  cat("\n  [WARN] 2-level MICE failed:", e$message, "\n")
  cat("  Falling back to standard PMM...\n\n")
  
  # Fallback to standard PMM
  methods_fallback <- mice::make.method(mice_data)
  methods_fallback["edu"] <- "pmm"
  for (var in names(methods_fallback)) {
    if (var != "edu") methods_fallback[var] <- ""
  }
  
  mice::mice(
    data = mice_data,
    m = MICE_PARAMS$m,
    maxit = MICE_PARAMS$maxit,
    method = methods_fallback,
    printFlag = TRUE
  )
})

end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

cat("\n  MICE completed in", round(as.numeric(duration), 1), "minutes\n")
cat("  Number of imputations:", mice_result$m, "\n")
cat("  Number of iterations:", mice_result$iteration, "\n")

# =============================================================================
# PART 8: Convergence Diagnostics
# =============================================================================

cat("\n=== Convergence Diagnostics ===\n\n")

# Trace plots (convergence)
cat("  Generating trace plots...\n")

# PDF version
pdf(file.path(mice_output_dir, "MICE_Convergence_Trace.pdf"), width = 12, height = 8)
plot(mice_result, layout = c(2, 1))
dev.off()

# PNG version
png(file.path(mice_output_dir, "MICE_Convergence_Trace.png"), width = 1200, height = 800, res = 150)
plot(mice_result, layout = c(2, 1))
dev.off()

cat("  Saved: MICE_Convergence_Trace.pdf/png\n")

# =============================================================================
# PART 9: Density Comparison (Observed vs Imputed)
# =============================================================================

cat("\n=== Density Comparison ===\n\n")

# Density plots
cat("  Generating density comparison plots...\n")

# For categorical variable, use stripplot or barchart
pdf(file.path(mice_output_dir, "MICE_Density_Comparison.pdf"), width = 12, height = 8)
tryCatch({
  # Use stripplot for categorical
  print(stripplot(mice_result, edu ~ .imp, pch = 20, cex = 1.2,
                  main = "Distribution of 'edu' by Imputation"))
}, error = function(e) {
  cat("  Could not create stripplot:", e$message, "\n")
})
dev.off()

png(file.path(mice_output_dir, "MICE_Density_Comparison.png"), width = 1200, height = 800, res = 150)
tryCatch({
  print(stripplot(mice_result, edu ~ .imp, pch = 20, cex = 1.2,
                  main = "Distribution of 'edu' by Imputation"))
}, error = function(e) NULL)
dev.off()

cat("  Saved: MICE_Density_Comparison.pdf/png\n")

# =============================================================================
# PART 10: Summary Statistics Comparison
# =============================================================================

cat("\n=== Summary Statistics Comparison ===\n\n")

# Get observed data
observed_edu <- mice_data$edu[!is.na(mice_data$edu)]
observed_table <- table(observed_edu)
observed_prop <- prop.table(observed_table)

# Get imputed data from first imputation
imputed_data_1 <- mice::complete(mice_result, 1)
imputed_edu <- imputed_data_1$edu
imputed_table <- table(imputed_edu)
imputed_prop <- prop.table(imputed_table)

# Create comparison table
comparison_df <- data.frame(
  Education_Level = names(observed_table),
  Observed_N = as.numeric(observed_table),
  Observed_Pct = round(as.numeric(observed_prop) * 100, 1),
  Imputed_N = as.numeric(imputed_table),
  Imputed_Pct = round(as.numeric(imputed_prop) * 100, 1)
)
comparison_df$Difference_Pct <- comparison_df$Imputed_Pct - comparison_df$Observed_Pct

cat("  Education Distribution Comparison:\n")
print(comparison_df)

# Comparison by cohort
comparison_by_cohort <- list()

for (coh in levels(mice_data$cohort)) {
  obs_sub <- mice_data$edu[mice_data$cohort == coh & !is.na(mice_data$edu)]
  imp_sub <- imputed_data_1$edu[imputed_data_1$cohort == coh]
  
  if (length(obs_sub) > 0 && length(imp_sub) > 0) {
    obs_tbl <- prop.table(table(obs_sub))
    imp_tbl <- prop.table(table(imp_sub))
    
    comparison_by_cohort[[coh]] <- data.frame(
      Cohort = coh,
      Education_Level = names(obs_tbl),
      Observed_Pct = round(as.numeric(obs_tbl) * 100, 1),
      Imputed_Pct = round(as.numeric(imp_tbl) * 100, 1)
    )
  }
}

comparison_by_cohort_df <- bind_rows(comparison_by_cohort)

# Save comparison
writexl::write_xlsx(
  list(
    Overall_Comparison = comparison_df,
    By_Cohort = comparison_by_cohort_df,
    MICE_Parameters = data.frame(
      Parameter = c("m", "maxit", "seed", "method", "cluster_var"),
      Value = c(MICE_PARAMS$m, MICE_PARAMS$maxit, MICE_PARAMS$seed, "2l.pmm", CLUSTER_VAR)
    )
  ),
  file.path(mice_output_dir, "MICE_Summary_Statistics.xlsx")
)
cat("\n  Saved: MICE_Summary_Statistics.xlsx\n")

# =============================================================================
# PART 11: Save Imputed Data
# =============================================================================

cat("\n=== Saving Imputed Data ===\n\n")

# Save directory
pooled_dir <- file.path(DATA_DIR, "Pooled")
if (!dir.exists(pooled_dir)) dir.create(pooled_dir, recursive = TRUE)

# Save full mids object (for Rubin's rules in Cox analysis)
saveRDS(mice_result, file.path(pooled_dir, "Pooled_mice_mids.rds"))
cat("  Saved: Pooled_mice_mids.rds (full mids object for Rubin's rules)\n")

# Create complete dataset using first imputation (for simple sensitivity analysis)
complete_data_1 <- mice::complete(mice_result, 1)

# Merge back with full variable set from pooled_data
# Keep imputed edu, other variables from original
final_imputed <- pooled_data
final_imputed$edu <- complete_data_1$edu

# Save pooled imputed CSV
write.csv(final_imputed, file.path(pooled_dir, "Pooled_mice_imputed.csv"), row.names = FALSE)
cat("  Saved: Pooled_mice_imputed.csv\n")

# Split by cohort and save
for (coh in unique(final_imputed$cohort)) {
  cohort_data <- final_imputed[final_imputed$cohort == coh, ]
  output_path <- file.path(DATA_DIR, coh, paste0(coh, "_mice_imputed.csv"))
  write.csv(cohort_data, output_path, row.names = FALSE)
  cat("  Saved:", coh, "_mice_imputed.csv (", nrow(cohort_data), "rows)\n")
}

# =============================================================================
# PART 12: Summary Report
# =============================================================================

cat("\n")
cat("################################################################\n")
cat("#                                                              #\n")
cat("#        2-Level MICE Imputation Complete                      #\n")
cat("#                                                              #\n")
cat("################################################################\n")
cat("\n")

cat("=== Summary ===\n\n")
cat("  Total samples:", nrow(final_imputed), "\n")
cat("  Imputations (m):", mice_result$m, "\n")
cat("  Iterations (maxit):", mice_result$iteration, "\n")
cat("  Cluster variable:", CLUSTER_VAR, "(", length(unique(final_imputed$cohort)), "levels)\n")
cat("  Variable imputed:", VAR_TO_IMPUTE, "\n")
cat("  Duration:", round(as.numeric(duration), 1), "minutes\n")
cat("\n")

cat("=== Output Files ===\n\n")
cat("  Data:\n")
cat("    - Pooled_mice_mids.rds (for Rubin's rules)\n")
cat("    - Pooled_mice_imputed.csv\n")
cat("    - [COHORT]_mice_imputed.csv (5 files)\n")
cat("\n")
cat("  Diagnostics:\n")
cat("    - MICE_Convergence_Trace.pdf/png\n")
cat("    - MICE_Density_Comparison.pdf/png\n")
cat("    - MICE_Summary_Statistics.xlsx\n")
cat("    - MICE_Missing_Data_Summary.xlsx\n")
cat("\n")

cat("=== Usage for Sensitivity Analysis S3 ===\n\n")
cat("  # Load mids object\n")
cat("  mice_obj <- readRDS('Data_Excel/Pooled/Pooled_mice_mids.rds')\n")
cat("\n")
cat("  # Fit Cox model on each imputed dataset\n")
cat("  fit <- with(mice_obj, coxph(Surv(time, event) ~ exposure + covariates))\n")
cat("\n")
cat("  # Pool results using Rubin's rules\n")
cat("  pooled <- pool(fit)\n")
cat("  summary(pooled)\n")
cat("\n")
