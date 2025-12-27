###############################################################################
# 02_MICE_Imputation.R
# ============================================================================
# Project: Lifestyle Factors and PPC-MM Association Study
# Purpose: 2-Level Multiple Imputation by Chained Equations (MICE)
# ============================================================================
#
# METHODOLOGY:
#   - 2-Level MICE with country_name as Level 2 clustering variable
#   - Pooled data imputation, then split back to individual cohorts
#   - Uses mice + miceadds packages for multilevel imputation
#
# PARAMETERS:
#   - Level 2: country_name (20 countries)
#   - Level 1: individual
#   - m = 20 (number of imputed datasets)
#   - maxit = 30 (iterations)
#   - Method: 2l.pmm (multilevel predictive mean matching)
#
# VARIABLES TO IMPUTE:
#   - edu (1.44% missing overall)
#   - employment (1.67% missing overall, optional)
#
# EXCLUDED VARIABLES:
#   - rural (18.66% missing, ELSA 100% missing)
#
###############################################################################

# =============================================================================
# PART 1: Environment Setup
# =============================================================================

source("C:/Users/user/Desktop/Scientific Research/UKB Project/Jung Sun Lab/Project1_Lifestyle_PPM_Hexiao,Hongtao/Code/00_Functions_and_Setup.R")

# Additional packages for MICE
if (!require("miceadds")) install.packages("miceadds")
library(miceadds)

cat("\n")
cat("================================================================\n")
cat("   2-Level MICE Imputation                                      \n")
cat("   Level 2: country_name | Level 1: individual                  \n")
cat("================================================================\n\n")

# =============================================================================
# PART 2: Configuration
# =============================================================================

# MICE parameters
MICE_CONFIG <- list(
  m = 20,           # Number of imputed datasets
  maxit = 30,       # Maximum iterations
  seed = 12345,     # Random seed for reproducibility
  printFlag = TRUE  # Print progress
)

# Variables to impute
VARS_TO_IMPUTE <- c("edu")  # Primary: only edu
# VARS_TO_IMPUTE <- c("edu", "employment")  # Optional: add employment

# Level 2 clustering variable
CLUSTER_VAR <- "country_name"

# Predictor variables (complete, no missing)
PREDICTOR_VARS <- c(
  "age_baseline",    # Continuous

  "sex",             # Binary
  "marital_binary",  # Binary
  "cohort",          # Categorical
  "region",          # Categorical
  "physical_base",   # Binary
  "psych_base",      # Binary
  "cog_base"         # Binary
)

# Variables to exclude from imputation model
EXCLUDE_VARS <- c("rural")  # ELSA 100% missing

cat("MICE Configuration:\n")
cat("  m (imputed datasets):", MICE_CONFIG$m, "\n")
cat("  maxit (iterations):", MICE_CONFIG$maxit, "\n")
cat("  seed:", MICE_CONFIG$seed, "\n")
cat("  Level 2 cluster:", CLUSTER_VAR, "\n")
cat("  Variables to impute:", paste(VARS_TO_IMPUTE, collapse = ", "), "\n")
cat("  Excluded variables:", paste(EXCLUDE_VARS, collapse = ", "), "\n")
cat("\n")

# =============================================================================
# PART 3: Load and Pool MICE Data
# =============================================================================

cat("--- Loading MICE data from all cohorts ---\n\n")

# Define MICE data paths
mice_paths <- list(
  CHARLS = file.path(DATA_DIR, "CHARLS", "CHARLS_mice_data.csv"),
  ELSA   = file.path(DATA_DIR, "ELSA", "ELSA_mice_data.csv"),
  HRS    = file.path(DATA_DIR, "HRS", "HRS_mice_data.csv"),
  MHAS   = file.path(DATA_DIR, "MHAS", "MHAS_mice_data.csv"),
  SHARE  = file.path(DATA_DIR, "SHARE", "SHARE_mice_data.csv")
)

# Load each cohort's MICE data
cohort_mice_list <- list()

for (coh in names(mice_paths)) {
  file_path <- mice_paths[[coh]]
  
  if (file.exists(file_path)) {
    df <- read.csv(file_path, stringsAsFactors = FALSE)
    cat("  Loaded", coh, ":", nrow(df), "rows,", ncol(df), "columns\n")
    cohort_mice_list[[coh]] <- df
  } else {
    cat("  WARNING: File not found for", coh, "\n")
  }
}

# Find common columns across all cohorts
common_cols <- Reduce(intersect, lapply(cohort_mice_list, names))
cat("\n  Common columns across all cohorts:", length(common_cols), "\n")

# Subset to common columns and combine
for (coh in names(cohort_mice_list)) {
  cohort_mice_list[[coh]] <- cohort_mice_list[[coh]][, common_cols]
}

pooled_mice <- bind_rows(cohort_mice_list)
cat("\n  Pooled MICE data:", nrow(pooled_mice), "rows\n")

# =============================================================================
# PART 4: Missing Data Summary
# =============================================================================

cat("\n--- Missing Data Summary ---\n\n")

# Calculate missing rates for key variables
missing_summary <- data.frame(
  Variable = character(),
  N_Total = integer(),
  N_Missing = integer(),
  Pct_Missing = numeric(),
  stringsAsFactors = FALSE
)

key_vars <- c("age_baseline", "sex", "edu", "marital_binary", "employment", 
              "rural", "cohort", "country_name", "region",
              "physical_base", "psych_base", "cog_base")

for (var in key_vars) {
  if (var %in% names(pooled_mice)) {
    n_total <- nrow(pooled_mice)
    n_missing <- sum(is.na(pooled_mice[[var]]) | pooled_mice[[var]] == "")
    pct_missing <- round(n_missing / n_total * 100, 2)
    
    missing_summary <- rbind(missing_summary, data.frame(
      Variable = var,
      N_Total = n_total,
      N_Missing = n_missing,
      Pct_Missing = pct_missing,
      stringsAsFactors = FALSE
    ))
    
    cat(sprintf("  %-20s: %d missing (%.2f%%)\n", var, n_missing, pct_missing))
  }
}

# Missing by cohort
cat("\n--- Missing by Cohort ---\n\n")

missing_by_cohort <- data.frame()

for (coh in names(cohort_mice_list)) {
  df <- cohort_mice_list[[coh]]
  n <- nrow(df)
  
  for (var in VARS_TO_IMPUTE) {
    if (var %in% names(df)) {
      n_miss <- sum(is.na(df[[var]]) | df[[var]] == "")
      pct <- round(n_miss / n * 100, 2)
      
      missing_by_cohort <- rbind(missing_by_cohort, data.frame(
        Cohort = coh,
        Variable = var,
        N = n,
        N_Missing = n_miss,
        Pct_Missing = pct
      ))
    }
  }
}

print(missing_by_cohort)

# Save missing data summary
missing_dir <- file.path(OUTPUT_DIR, "01_Descriptive")
if (!dir.exists(missing_dir)) dir.create(missing_dir, recursive = TRUE)

writexl::write_xlsx(
  list(
    Overall = missing_summary,
    By_Cohort = missing_by_cohort
  ),
  file.path(missing_dir, "Missing_Data_Summary.xlsx")
)

cat("\n  Saved: Missing_Data_Summary.xlsx\n")

# =============================================================================
# PART 5: Prepare Data for 2-Level MICE
# =============================================================================

cat("\n--- Preparing Data for 2-Level MICE ---\n\n")

# Select variables for imputation model
impute_model_vars <- c(
  CLUSTER_VAR,       # Level 2 cluster
  VARS_TO_IMPUTE,    # Variables to impute
  PREDICTOR_VARS,    # Complete predictors
  "cohort"           # Keep for later splitting
)

# Remove variables not in dataset
impute_model_vars <- impute_model_vars[impute_model_vars %in% names(pooled_mice)]

# Create imputation dataset
mice_data <- pooled_mice[, impute_model_vars]

# Convert categorical variables to factors
factor_vars <- c("sex", "edu", "marital_binary", "cohort", "region", 
                 "physical_base", "psych_base", "cog_base", CLUSTER_VAR)

for (var in factor_vars) {
  if (var %in% names(mice_data)) {
    mice_data[[var]] <- as.factor(mice_data[[var]])
  }
}

cat("  Imputation dataset: ", nrow(mice_data), "rows,", ncol(mice_data), "columns\n")
cat("  Level 2 clusters (countries):", length(unique(mice_data[[CLUSTER_VAR]])), "\n")
cat("\n")

# =============================================================================
# PART 6: Run 2-Level MICE Imputation
# =============================================================================

cat("--- Running 2-Level MICE Imputation ---\n")
cat("  This may take several minutes...\n\n")

# Set up imputation methods
# -2 indicates the cluster variable
# 2l.pmm for multilevel predictive mean matching

# Create predictor matrix
pred_matrix <- mice::make.predictorMatrix(mice_data)

# Set cluster variable
pred_matrix[, CLUSTER_VAR] <- -2  # -2 indicates cluster variable

# Don't use cluster variable to predict itself
pred_matrix[CLUSTER_VAR, ] <- 0

# Set up methods
methods <- mice::make.method(mice_data)

# For variables to impute, use 2-level methods
for (var in VARS_TO_IMPUTE) {
  if (var %in% names(mice_data)) {
    # Check if binary or multi-category
    n_levels <- length(unique(na.omit(mice_data[[var]])))
    if (n_levels == 2) {
      methods[var] <- "2l.bin"  # Binary
    } else {
      methods[var] <- "2l.pmm"  # Multi-category: use pmm
    }
  }
}

# For complete variables, set method to empty
for (var in PREDICTOR_VARS) {
  if (var %in% names(methods)) {
    methods[var] <- ""
  }
}
methods[CLUSTER_VAR] <- ""
methods["cohort"] <- ""

cat("Imputation methods:\n")
print(methods[methods != ""])
cat("\n")

# Run MICE
set.seed(MICE_CONFIG$seed)

tryCatch({
  mice_result <- mice::mice(
    data = mice_data,
    m = MICE_CONFIG$m,
    maxit = MICE_CONFIG$maxit,
    method = methods,
    predictorMatrix = pred_matrix,
    printFlag = MICE_CONFIG$printFlag
  )
  
  cat("\n  MICE imputation completed successfully!\n")
  cat("  Number of imputed datasets:", mice_result$m, "\n")
  
}, error = function(e) {
  cat("\n  ERROR in MICE imputation:", e$message, "\n")
  cat("  Falling back to standard MICE (single-level)...\n")
  
  # Fallback to standard MICE
  methods_fallback <- mice::make.method(mice_data)
  for (var in VARS_TO_IMPUTE) {
    if (var %in% names(mice_data)) {
      methods_fallback[var] <- "pmm"
    }
  }
  for (var in c(PREDICTOR_VARS, CLUSTER_VAR, "cohort")) {
    if (var %in% names(methods_fallback)) {
      methods_fallback[var] <- ""
    }
  }
  
  mice_result <<- mice::mice(
    data = mice_data,
    m = MICE_CONFIG$m,
    maxit = MICE_CONFIG$maxit,
    method = methods_fallback,
    printFlag = MICE_CONFIG$printFlag
  )
  
  cat("  Standard MICE completed.\n")
})

# =============================================================================
# PART 7: Diagnostics
# =============================================================================

cat("\n--- Imputation Diagnostics ---\n\n")

# Convergence plots
diag_dir <- file.path(OUTPUT_DIR, "01_Descriptive")

pdf(file.path(diag_dir, "MICE_Convergence_Plot.pdf"), width = 10, height = 6)
plot(mice_result)
dev.off()
cat("  Saved: MICE_Convergence_Plot.pdf\n")

# Density plots for imputed variables
pdf(file.path(diag_dir, "MICE_Density_Plot.pdf"), width = 10, height = 6)
for (var in VARS_TO_IMPUTE) {
  if (var %in% names(mice_data)) {
    tryCatch({
      densityplot(mice_result, ~ get(var))
    }, error = function(e) {
      cat("  Could not create density plot for", var, "\n")
    })
  }
}
dev.off()
cat("  Saved: MICE_Density_Plot.pdf\n")

# =============================================================================
# PART 8: Create Complete Datasets and Split by Cohort
# =============================================================================

cat("\n--- Creating Imputed Datasets ---\n\n")

# Get the first complete dataset for main analysis
# (For full analysis, use with() and pool() functions)
imputed_complete <- mice::complete(mice_result, action = 1)

# Merge back with original pooled data (to get all variables)
# Keep only the imputed variables from mice_result, rest from original

# Get row indices
original_vars <- setdiff(names(pooled_mice), VARS_TO_IMPUTE)
original_vars <- original_vars[original_vars %in% names(pooled_mice)]

# Create final imputed dataset
pooled_imputed <- pooled_mice[, original_vars]
for (var in VARS_TO_IMPUTE) {
  if (var %in% names(imputed_complete)) {
    pooled_imputed[[var]] <- imputed_complete[[var]]
  }
}

cat("  Pooled imputed data:", nrow(pooled_imputed), "rows,", ncol(pooled_imputed), "columns\n")

# Verify no missing in imputed variables
cat("\n  Verification - Missing after imputation:\n")
for (var in VARS_TO_IMPUTE) {
  if (var %in% names(pooled_imputed)) {
    n_miss <- sum(is.na(pooled_imputed[[var]]) | pooled_imputed[[var]] == "")
    cat(sprintf("    %-15s: %d missing\n", var, n_miss))
  }
}

# =============================================================================
# PART 9: Split Back to Individual Cohorts
# =============================================================================

cat("\n--- Splitting Imputed Data by Cohort ---\n\n")

# Save pooled imputed data
pooled_dir <- file.path(DATA_DIR, "Pooled")
if (!dir.exists(pooled_dir)) dir.create(pooled_dir, recursive = TRUE)

# Save as RDS (preserves R object with all m imputations)
saveRDS(mice_result, file.path(pooled_dir, "Pooled_mice_imputed.rds"))
cat("  Saved: Pooled_mice_imputed.rds (full mice object with m=20 imputations)\n")

# Save first imputation as CSV
write.csv(pooled_imputed, file.path(pooled_dir, "Pooled_mice_imputed.csv"), row.names = FALSE)
cat("  Saved: Pooled_mice_imputed.csv (first imputation)\n")

# Split by cohort and save
cohorts <- unique(pooled_imputed$cohort)

for (coh in cohorts) {
  cohort_data <- pooled_imputed[pooled_imputed$cohort == coh, ]
  
  output_path <- file.path(DATA_DIR, coh, paste0(coh, "_mice_imputed.csv"))
  write.csv(cohort_data, output_path, row.names = FALSE)
  
  cat("  Saved:", coh, "-", nrow(cohort_data), "rows\n")
}

# =============================================================================
# PART 10: Summary Report
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   2-Level MICE Imputation Complete                             \n")
cat("================================================================\n\n")

cat("Summary:\n")
cat("  - Total samples imputed:", nrow(pooled_imputed), "\n")
cat("  - Number of imputations (m):", MICE_CONFIG$m, "\n")
cat("  - Level 2 clusters:", length(unique(pooled_imputed[[CLUSTER_VAR]])), "countries\n")
cat("  - Variables imputed:", paste(VARS_TO_IMPUTE, collapse = ", "), "\n")
cat("\n")

cat("Output files:\n")
cat("  - Pooled_mice_imputed.rds (full mice object)\n")
cat("  - Pooled_mice_imputed.csv (first imputation)\n")
for (coh in cohorts) {
  cat("  -", paste0(coh, "_mice_imputed.csv\n"))
}
cat("  - Missing_Data_Summary.xlsx\n")
cat("  - MICE_Convergence_Plot.pdf\n")
cat("  - MICE_Density_Plot.pdf\n")

cat("\n")
cat("Note: For pooled analysis with Rubin's rules, use:\n")
cat("  mice_obj <- readRDS('Pooled_mice_imputed.rds')\n")
cat("  fit <- with(mice_obj, coxph(...))\n")
cat("  pooled_result <- pool(fit)\n")
cat("\n")

