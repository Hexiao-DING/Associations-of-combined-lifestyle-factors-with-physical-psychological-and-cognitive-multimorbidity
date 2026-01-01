###############################################################################
# 04_Pooled_Cox_Analysis.R
# ============================================================================
# Project: Lifestyle Factors and PPC-MM Association Study
# Purpose: Comprehensive Cox proportional hazards regression (Lancet standard)
# ============================================================================
#
# ANALYSIS STRUCTURE:
#
# === PRIMARY ANALYSIS (Main Data) ===
#   A. Individual Lifestyle Factors (mutually adjusted)
#      - 4 factors: drink, smoke, PA, social
#      - Outcomes: Overall PPC-MM (Primary) + 4 subtypes (Secondary)
#
#   B. Cumulative Lifestyle Effect (4-level: 0/1/2/3+)
#      - P-trend test
#      - Dose-response plots
#
# === SENSITIVITY ANALYSES ===
#   S1: 5-level categories (0/1/2/3/4)
#   S2: Heavy drinking definition (replace unhealthy_drink with heavy_drink)
#   S3: MICE imputed data
#   S4: Drop1st (exclude first follow-up wave)
#
###############################################################################

# =============================================================================
# PART 1: Environment Setup
# =============================================================================

source("C:/Users/user/Desktop/Scientific Research/UKB Project/Jung Sun Lab/Project1_Lifestyle_PPM_Hexiao,Hongtao/Code/00_Functions_and_Setup.R")

cat("\n")
cat("================================================================\n")
cat("   Pooled Cox Regression Analysis (Comprehensive)               \n")
cat("   Lancet Standard Analysis Framework                           \n")
cat("================================================================\n\n")

# =============================================================================
# PART 2: Define Outcome Information
# =============================================================================

# Outcome definitions (Overall is PRIMARY, subtypes are SECONDARY)
outcome_info <- list(
  # PRIMARY OUTCOME
  Overall = list(
    event = "event_ppcmm", 
    time = "time_ppcmm_months",
    label = "Overall PPC-MM",
    type = "Primary"
  ),
  # SECONDARY OUTCOMES (4 subtypes)
  P1P2 = list(
    event = "event_mm_phys_psych", 
    time = "time_mm_phys_psych",
    label = "Physical-Psychological (P1P2)",
    type = "Secondary"
  ),
  P1C = list(
    event = "event_mm_phys_cog", 
    time = "time_mm_phys_cog",
    label = "Physical-Cognitive (P1C)",
    type = "Secondary"
  ),
  P2C = list(
    event = "event_mm_psych_cog", 
    time = "time_mm_psych_cog",
    label = "Psychological-Cognitive (P2C)",
    type = "Secondary"
  ),
  P1P2C = list(
    event = "event_mm_all_three", 
    time = "time_mm_all_three",
    label = "All Three (P1P2C)",
    type = "Secondary"
  )
)

all_outcomes <- names(outcome_info)

# =============================================================================
# PART 3: Data Loading Functions
# =============================================================================

#' Load and standardize data from all cohorts
#' @param data_paths List of file paths
#' @param data_type Type of data (main, sensitivity, mice)
#' @return Pooled data frame
load_and_pool_data <- function(data_paths, data_type = "main") {
  
  cat("--- Loading", data_type, "data from all cohorts ---\n")
  
  pooled_list <- list()
  
  for (coh in names(data_paths)) {
    file_path <- data_paths[[coh]]
    
    if (is.null(file_path) || !file.exists(file_path)) {
      cat("  File not found for", coh, ", skipping\n")
      next
    }
    
    df <- read.csv(file_path, stringsAsFactors = FALSE)
    bw <- BASELINE_WAVES[[coh]]
    
    cat("  Loaded", coh, "- N =", nrow(df), "\n")
    
    # Standardize variable names
    df_std <- df %>%
      mutate(
        cohort = coh,
        country = COHORT_COUNTRIES[[coh]]
      )
    
    # Rename wave-specific lifestyle variables
    score_var <- paste0("w", bw, "_unhealthy_score")
    drink_var <- paste0("w", bw, "_unhealthy_drink")
    smoke_var <- paste0("w", bw, "_unhealthy_smoke")
    pa_var <- paste0("w", bw, "_unhealthy_pa")
    soc_var <- paste0("w", bw, "_unhealthy_soc")
    cat_var <- paste0("w", bw, "_unhealthy_cat")
    heavy_drink_var <- paste0("w", bw, "_heavy_drink")
    score_heavy_var <- paste0("w", bw, "_unhealthy_score_heavy")
    
    # Debug: show which variables will be renamed
    cat("    Baseline wave:", bw, "- Looking for:", score_var, "\n")
    cat("    Variables found:", sum(c(score_var, drink_var, smoke_var, pa_var, soc_var) %in% names(df_std)), "/5\n")
    
    # Rename variables using base R (more robust than dplyr rename with !!)
    rename_var <- function(data, old_name, new_name) {
      if (old_name %in% names(data)) {
        names(data)[names(data) == old_name] <- new_name
        cat("      ✓ Renamed:", old_name, "->", new_name, "\n")
      }
      return(data)
    }
    
    df_std <- rename_var(df_std, score_var, "unhealthy_score")
    df_std <- rename_var(df_std, drink_var, "unhealthy_drink")
    df_std <- rename_var(df_std, smoke_var, "unhealthy_smoke")
    df_std <- rename_var(df_std, pa_var, "unhealthy_pa")
    df_std <- rename_var(df_std, soc_var, "unhealthy_soc")
    df_std <- rename_var(df_std, cat_var, "unhealthy_cat")
    df_std <- rename_var(df_std, heavy_drink_var, "heavy_drink")
    df_std <- rename_var(df_std, score_heavy_var, "unhealthy_score_heavy")
    
    # Select common variables for pooling
    # IMPORTANT: Include country_name and region from data
    common_vars <- c("cohort", "country", "country_name", "region",
                     "age_baseline", "sex", "edu", "marital_binary", "employment", "rural",
                     "unhealthy_score", "unhealthy_drink", "unhealthy_smoke",
                     "unhealthy_pa", "unhealthy_soc", "unhealthy_cat",
                     "heavy_drink", "unhealthy_score_heavy",
                     "event_ppcmm", "time_ppcmm_months",
                     "event_ppcmm_drop1st", "time_ppcmm_drop1st",
                     "event_mm_phys_psych", "event_mm_phys_cog",
                     "event_mm_psych_cog", "event_mm_all_three",
                     "time_mm_phys_psych", "time_mm_phys_cog",
                     "time_mm_psych_cog", "time_mm_all_three",
                     "physical_base", "psych_base", "cog_base")
    
    existing_vars <- common_vars[common_vars %in% names(df_std)]
    pooled_list[[coh]] <- df_std %>% select(all_of(existing_vars))
  }
  
  # Combine all cohorts
  pooled_data <- bind_rows(pooled_list)
  
  cat("\nPooled data summary:\n")
  cat("  Total sample size:", nrow(pooled_data), "\n")
  cat("  Cohorts:", paste(unique(pooled_data$cohort), collapse = ", "), "\n")
  
  return(pooled_data)
}

#' Prepare analysis variables
#' @param data Pooled data frame
#' @param use_heavy_drink Whether to use heavy_drink instead of unhealthy_drink
#' @return Prepared data frame
prepare_analysis_data <- function(data, use_heavy_drink = FALSE) {
  
  cat("\n--- Preparing analysis data ---\n")
  cat("Input rows:", nrow(data), "\n")
  
  df <- data %>%
    mutate(
      sex = relevel(factor(sex), ref = "Men"),
      edu = relevel(factor(edu), ref = "Primary or below"),
      cohort = factor(cohort),
      
      # Main analysis: 4-level categorical (0/1/2/3+)
      n_lifestyle_cat = factor(
        case_when(
          unhealthy_score == 0 ~ "0",
          unhealthy_score == 1 ~ "1",
          unhealthy_score == 2 ~ "2",
          unhealthy_score >= 3 ~ "3+"
        ), 
        levels = c("0", "1", "2", "3+")
      ),
      
      # 5-level score (0/1/2/3/4)
      n_lifestyle_5cat = factor(unhealthy_score, levels = 0:4),
      
      # Continuous for trend test
      n_lifestyle_continuous = unhealthy_score
    )
  
  # Make region a factor if it exists
  if ("region" %in% names(df)) {
    df$region <- factor(df$region)
  }
  
  # Heavy drink version if needed
  if (use_heavy_drink && "heavy_drink" %in% names(df) && "unhealthy_score_heavy" %in% names(df)) {
    df <- df %>%
      mutate(
        n_lifestyle_cat_heavy = factor(
          case_when(
            unhealthy_score_heavy == 0 ~ "0",
            unhealthy_score_heavy == 1 ~ "1",
            unhealthy_score_heavy == 2 ~ "2",
            unhealthy_score_heavy >= 3 ~ "3+"
          ), 
          levels = c("0", "1", "2", "3+")
        ),
        n_lifestyle_5cat_heavy = factor(unhealthy_score_heavy, levels = 0:4),
        n_lifestyle_continuous_heavy = unhealthy_score_heavy
      )
  }
  
  # Complete case filter - DON'T filter on region (it should exist in data)
  df <- df %>%
    filter(!is.na(age_baseline),
           !is.na(sex),
           !is.na(edu),
           !is.na(unhealthy_drink),
           !is.na(unhealthy_smoke),
           !is.na(unhealthy_pa),
           !is.na(unhealthy_soc))
  
  cat("After complete case filter:", nrow(df), "rows\n")
  
  # Show variable distributions
  cat("Sex distribution:\n")
  print(table(df$sex, useNA = "ifany"))
  cat("Education distribution:\n")
  print(table(df$edu, useNA = "ifany"))
  cat("Lifestyle category distribution:\n")
  print(table(df$n_lifestyle_cat, useNA = "ifany"))
  
  return(df)
}

# =============================================================================
# PART 4: Cox Model Functions
# =============================================================================

#' Run Cox model and extract results
#' @param data Analysis dataset
#' @param outcome_name Outcome name from outcome_info
#' @param exposure_var Exposure variable name
#' @param covariates Covariate names
#' @param model_label Model description
#' @param calc_trend Whether to calculate P for trend
#' @param use_strata Whether to use strata(cohort) for pooled analysis
#' @return Data frame with results
run_cox_model <- function(data, outcome_name, exposure_var, covariates, 
                          model_label, calc_trend = FALSE, use_strata = TRUE) {
  
  out_info <- outcome_info[[outcome_name]]
  event_var <- out_info$event
  time_var <- out_info$time
  
  # Check if variables exist
  if (!event_var %in% names(data) || !time_var %in% names(data)) {
    cat("    [SKIP] Outcome variables not found:", event_var, "or", time_var, "\n")
    return(NULL)
  }
  
  # Filter valid observations
  df <- data %>%
    filter(!is.na(.data[[event_var]]) & !is.na(.data[[time_var]]) & .data[[time_var]] > 0)
  
  n_total <- nrow(df)
  n_events <- sum(df[[event_var]], na.rm = TRUE)
  
  cat("    N =", n_total, ", Events =", n_events, "\n")
  
  if (n_events < 10) {
    cat("    [SKIP] Insufficient events (<10)\n")
    return(NULL)
  }
  
  # Check which covariates exist in data
  available_covs <- covariates[covariates %in% names(df)]
  if (length(available_covs) < length(covariates)) {
    missing_covs <- setdiff(covariates, available_covs)
    cat("    [WARN] Missing covariates:", paste(missing_covs, collapse = ", "), "\n")
  }
  
  # Build formula with stratification
  cov_str <- paste(available_covs, collapse = " + ")
  
  # Add strata(cohort) for pooled analysis
  if (use_strata && "cohort" %in% names(df)) {
    formula_str <- paste0("Surv(", time_var, ", ", event_var, ") ~ ", 
                          exposure_var, " + ", cov_str, " + strata(cohort)")
  } else {
    formula_str <- paste0("Surv(", time_var, ", ", event_var, ") ~ ", 
                          exposure_var, " + ", cov_str)
  }
  
  cat("    Formula:", formula_str, "\n")
  
  # Fit model
  fit <- tryCatch({
    coxph(as.formula(formula_str), data = df)
  }, error = function(e) {
    cat("    [ERROR] Cox model failed:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(fit)) return(NULL)
  
  cat("    [OK] Model fitted successfully\n")
  
  # Extract results
  s <- summary(fit)
  coef_mat <- s$coefficients
  conf_mat <- s$conf.int
  
  # Get exposure coefficients only
  exp_rows <- grepl(exposure_var, rownames(coef_mat))
  
  if (sum(exp_rows) == 0) return(NULL)
  
  results <- data.frame(
    Model = model_label,
    Outcome = outcome_name,
    Outcome_Label = out_info$label,
    Outcome_Type = out_info$type,
    Variable = rownames(coef_mat)[exp_rows],
    HR = round(coef_mat[exp_rows, "exp(coef)"], 3),
    Lower_CI = round(conf_mat[exp_rows, "lower .95"], 3),
    Upper_CI = round(conf_mat[exp_rows, "upper .95"], 3),
    P_value = format.pval(coef_mat[exp_rows, "Pr(>|z|)"], digits = 3),
    N = n_total,
    N_Events = n_events,
    stringsAsFactors = FALSE
  )
  
  # Calculate P for trend if requested
  if (calc_trend) {
    # Create continuous version for trend
    trend_var <- gsub("_cat$|_5cat$", "_continuous", exposure_var)
    if (grepl("heavy", exposure_var)) {
      trend_var <- "n_lifestyle_continuous_heavy"
    } else {
      trend_var <- "n_lifestyle_continuous"
    }
    
    if (trend_var %in% names(df)) {
      # Add stratification for trend test too
      if (use_strata && "cohort" %in% names(df)) {
        trend_formula <- paste0("Surv(", time_var, ", ", event_var, ") ~ ", 
                                trend_var, " + ", cov_str, " + strata(cohort)")
      } else {
        trend_formula <- paste0("Surv(", time_var, ", ", event_var, ") ~ ", 
                                trend_var, " + ", cov_str)
      }
      trend_fit <- tryCatch({
        coxph(as.formula(trend_formula), data = df)
      }, error = function(e) NULL)
      
      if (!is.null(trend_fit)) {
        trend_p <- summary(trend_fit)$coefficients[trend_var, "Pr(>|z|)"]
        results$P_trend <- format.pval(trend_p, digits = 3)
      }
    }
  }
  
  return(results)
}

#' Run individual factors analysis (mutually adjusted)
#' @param data Analysis dataset
#' @param outcome_name Outcome name
#' @param covariates Covariate names
#' @param use_heavy_drink Use heavy_drink instead of unhealthy_drink
#' @param use_strata Whether to use strata(cohort) for pooled analysis
#' @return Data frame with results
run_individual_factors <- function(data, outcome_name, covariates, 
                                   use_heavy_drink = FALSE, use_strata = TRUE) {
  
  out_info <- outcome_info[[outcome_name]]
  event_var <- out_info$event
  time_var <- out_info$time
  
  if (!event_var %in% names(data) || !time_var %in% names(data)) {
    cat("    [SKIP] Outcome variables not found\n")
    return(NULL)
  }
  
  df <- data %>%
    filter(!is.na(.data[[event_var]]) & !is.na(.data[[time_var]]) & .data[[time_var]] > 0)
  
  n_total <- nrow(df)
  n_events <- sum(df[[event_var]], na.rm = TRUE)
  
  cat("    N =", n_total, ", Events =", n_events, "\n")
  
  if (n_events < 10) {
    cat("    [SKIP] Insufficient events (<10)\n")
    return(NULL)
  }
  
  # Select drink variable
  drink_var <- ifelse(use_heavy_drink && "heavy_drink" %in% names(df), 
                      "heavy_drink", "unhealthy_drink")
  
  lifestyle_vars <- c(drink_var, "unhealthy_smoke", "unhealthy_pa", "unhealthy_soc")
  
  # Check which covariates exist in data
  available_covs <- covariates[covariates %in% names(df)]
  if (length(available_covs) < length(covariates)) {
    missing_covs <- setdiff(covariates, available_covs)
    cat("    [WARN] Missing covariates:", paste(missing_covs, collapse = ", "), "\n")
  }
  
  # Build formula with all 4 factors mutually adjusted + stratification
  if (use_strata && "cohort" %in% names(df)) {
    formula_str <- paste0("Surv(", time_var, ", ", event_var, ") ~ ",
                          paste(lifestyle_vars, collapse = " + "), " + ",
                          paste(available_covs, collapse = " + "), " + strata(cohort)")
  } else {
    formula_str <- paste0("Surv(", time_var, ", ", event_var, ") ~ ",
                          paste(lifestyle_vars, collapse = " + "), " + ",
                          paste(available_covs, collapse = " + "))
  }
  
  cat("    Formula:", formula_str, "\n")
  
  fit <- tryCatch({
    coxph(as.formula(formula_str), data = df)
  }, error = function(e) {
    cat("    [ERROR] Cox model failed:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(fit)) return(NULL)
  
  cat("    [OK] Model fitted successfully\n")
  
  s <- summary(fit)
  coef_mat <- s$coefficients
  conf_mat <- s$conf.int
  
  # Get lifestyle variable rows
  exp_rows <- rownames(coef_mat) %in% lifestyle_vars
  
  # Create nice labels
  var_labels <- c(
    "unhealthy_drink" = "Drinking (any)",
    "heavy_drink" = "Heavy Drinking",
    "unhealthy_smoke" = "Smoking",
    "unhealthy_pa" = "Physical Inactivity",
    "unhealthy_soc" = "Social Isolation"
  )
  
  results <- data.frame(
    Model = ifelse(use_heavy_drink, "Individual Factors (Heavy Drink)", "Individual Factors"),
    Outcome = outcome_name,
    Outcome_Label = out_info$label,
    Outcome_Type = out_info$type,
    Variable = rownames(coef_mat)[exp_rows],
    Variable_Label = var_labels[rownames(coef_mat)[exp_rows]],
    HR = round(coef_mat[exp_rows, "exp(coef)"], 3),
    Lower_CI = round(conf_mat[exp_rows, "lower .95"], 3),
    Upper_CI = round(conf_mat[exp_rows, "upper .95"], 3),
    P_value = format.pval(coef_mat[exp_rows, "Pr(>|z|)"], digits = 3),
    N = nrow(df),
    N_Events = n_events,
    stringsAsFactors = FALSE
  )
  
  return(results)
}

# =============================================================================
# PART 5: Load Main Data
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   LOADING DATA                                                  \n")
cat("================================================================\n")

# FIXED: Load standardized pooled data from 01_Pooled_Descriptive.R output
# This ensures variable names are consistent
pooled_rds <- file.path(OUTPUT_DIR, "Pooled_main_data.rds")
pooled_csv <- file.path(OUTPUT_DIR, "Pooled_main_data.csv")

if (file.exists(pooled_rds)) {
  cat("Loading standardized pooled data from RDS...\n")
  pooled_main <- readRDS(pooled_rds)
  cat("  [OK] Loaded from:", pooled_rds, "\n")
  cat("  Sample size:", nrow(pooled_main), "\n")
} else if (file.exists(pooled_csv)) {
  cat("Loading standardized pooled data from CSV...\n")
  pooled_main <- read.csv(pooled_csv, stringsAsFactors = FALSE)
  cat("  [OK] Loaded from:", pooled_csv, "\n")
  cat("  Sample size:", nrow(pooled_main), "\n")
} else {
  cat("  [WARNING] Standardized pooled data not found, loading from raw files...\n")
  pooled_main <- load_and_pool_data(DATA_PATHS, "main")
}

cat("  Columns:", length(names(pooled_main)), "\n")
cat("  Key variables: unhealthy_score =", "unhealthy_score" %in% names(pooled_main),
    ", event_ppcmm =", "event_ppcmm" %in% names(pooled_main), "\n")

# Check if region exists in data
cat("\n--- Checking key variables in pooled data ---\n")
cat("Variables in pooled data:", paste(head(names(pooled_main), 30), collapse = ", "), "...\n")
cat("region exists:", "region" %in% names(pooled_main), "\n")
cat("country_name exists:", "country_name" %in% names(pooled_main), "\n")

# If region doesn't exist, create it
if (!"region" %in% names(pooled_main) || all(is.na(pooled_main$region))) {
  cat("\n[INFO] Creating 'region' variable from cohort/country information...\n")
  
  # First, try using country_name if available
  if ("country_name" %in% names(pooled_main)) {
    pooled_main <- pooled_main %>%
      mutate(
        region = case_when(
          country_name %in% c("Austria", "Switzerland", "France", "Belgium", 
                              "Netherlands", "Germany", "Luxembourg") ~ "Western Europe",
          country_name %in% c("Estonia", "Sweden", "Denmark") ~ "Northern Europe",
          country_name %in% c("Spain", "Italy", "Portugal", "Slovenia", "Greece") ~ "Southern Europe",
          country_name %in% c("Czech Republic", "Poland", "Hungary") ~ "Eastern Europe",
          country_name == "England" ~ "Northern Europe",
          country_name == "China" ~ "Eastern Asia",
          country_name == "USA" ~ "Northern America",
          country_name == "Mexico" ~ "Central America",
          TRUE ~ "Other"
        )
      )
  } else {
    # Fallback: use cohort names
    pooled_main <- pooled_main %>%
      mutate(
        region = case_when(
          cohort == "CHARLS" ~ "Eastern Asia",
          cohort == "ELSA" ~ "Northern Europe",
          cohort == "HRS" ~ "Northern America",
          cohort == "MHAS" ~ "Central America",
          cohort == "SHARE" ~ "Europe",
          TRUE ~ "Other"
        )
      )
  }
  cat("  Region variable created.\n")
} else {
  cat("\n[INFO] Region variable already exists in data.\n")
}

# Check region distribution
cat("\nRegion distribution:\n")
print(table(pooled_main$region, useNA = "ifany"))

# Prepare analysis data
pooled_analysis <- prepare_analysis_data(pooled_main, use_heavy_drink = TRUE)

# Ensure region is a factor
if ("region" %in% names(pooled_analysis)) {
  pooled_analysis$region <- factor(pooled_analysis$region)
  cat("\nRegion levels:", paste(levels(pooled_analysis$region), collapse = ", "), "\n")
}

cat("\nComplete case sample size:", nrow(pooled_analysis), "\n")

# Debug: show variable summary
cat("\n--- Analysis Variables Summary ---\n")
cat("N:", nrow(pooled_analysis), "\n")
cat("Cohorts:", paste(unique(pooled_analysis$cohort), collapse = ", "), "\n")
cat("Events (Overall):", sum(pooled_analysis$event_ppcmm, na.rm = TRUE), "\n")
cat("Lifestyle score range:", min(pooled_analysis$unhealthy_score, na.rm = TRUE), "-", 
    max(pooled_analysis$unhealthy_score, na.rm = TRUE), "\n")

# Print event counts
cat("\nEvent counts by outcome:\n")
for (out_name in all_outcomes) {
  event_var <- outcome_info[[out_name]]$event
  if (event_var %in% names(pooled_analysis)) {
    n_events <- sum(pooled_analysis[[event_var]], na.rm = TRUE)
    cat("  ", out_name, "(", outcome_info[[out_name]]$type, "):", n_events, "events\n")
  }
}

# =============================================================================
# Initialize all result data frames (to prevent errors if some analyses fail)
# =============================================================================

cat("\n[INFO] Initializing result containers...\n")
# Primary analyses
individual_results_df <- data.frame()
cumulative_results_df <- data.frame()
# Sensitivity analyses
sens1_results_df <- data.frame()
sens2a_results_df <- data.frame()
sens2b_results_df <- data.frame()
sens3_individual_df <- data.frame()
sens3_cumulative_df <- data.frame()
sens4_results_df <- data.frame()
# Subgroup analyses (IMPORTANT: must initialize)
age_stratified_df <- data.frame()
sex_stratified_df <- data.frame()
interaction_df <- data.frame()
all_results <- data.frame()
cat("[INFO] All result containers initialized.\n")

# =============================================================================
# PART 6: PRIMARY ANALYSIS - Individual Lifestyle Factors
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   PRIMARY ANALYSIS A: Individual Lifestyle Factors              \n")
cat("   (4 factors mutually adjusted)                                 \n")
cat("================================================================\n")

# Use centralized covariate configuration from 00_Functions_and_Setup.R
# Covariates: age_baseline + sex + edu + region + strata(cohort)
# Check which covariates actually exist in the data
available_covariates <- COX_COVARIATES$pooled[COX_COVARIATES$pooled %in% names(pooled_analysis)]
strata_var <- COX_COVARIATES$strata_var

# If region is not available, remove it from covariates
if (!"region" %in% names(pooled_analysis)) {
  available_covariates <- available_covariates[available_covariates != "region"]
  cat("[WARNING] 'region' variable not found, using reduced covariate set.\n")
}

# Use available covariates
covariates <- available_covariates

cat("\nCox Model Specification:\n")
cat("  Covariates:", paste(covariates, collapse = " + "), "\n")
cat("  Stratification:", paste0("strata(", strata_var, ")"), "\n")
cat("  Available variables in data:", paste(names(pooled_analysis)[1:20], collapse = ", "), "...\n\n")

individual_results <- list()

for (out_name in all_outcomes) {
  cat("\n--- Outcome:", outcome_info[[out_name]]$label, 
      "(", outcome_info[[out_name]]$type, ") ---\n")
  
  result <- run_individual_factors(pooled_analysis, out_name, covariates, use_heavy_drink = FALSE)
  
  if (!is.null(result)) {
    individual_results[[out_name]] <- result
    print(result %>% select(Variable_Label, HR, Lower_CI, Upper_CI, P_value))
  }
}

individual_results_df <- bind_rows(individual_results)

# IMMEDIATE SAVE: Save primary individual results right away
cat("\n--- Immediate Save: Primary Individual Results ---\n")
tryCatch({
  if (nrow(individual_results_df) > 0) {
    write.csv(individual_results_df, file.path(OUTPUT_DIR, "Cox_Primary_Individual_Factors.csv"), row.names = FALSE)
    cat("  [OK] Saved Cox_Primary_Individual_Factors.csv (", nrow(individual_results_df), " rows)\n")
  } else {
    cat("  [WARN] No individual results to save\n")
  }
}, error = function(e) {
  cat("  [ERROR] Failed to save individual results:", e$message, "\n")
})

# =============================================================================
# PART 7: PRIMARY ANALYSIS - Cumulative Effect (4-level)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   PRIMARY ANALYSIS B: Cumulative Lifestyle Effect               \n")
cat("   (0/1/2/3+ categories, Reference = 0)                          \n")
cat("================================================================\n")

cumulative_results <- list()

for (out_name in all_outcomes) {
  cat("\n--- Outcome:", outcome_info[[out_name]]$label, 
      "(", outcome_info[[out_name]]$type, ") ---\n")
  
  result <- run_cox_model(
    pooled_analysis, 
    out_name, 
    "n_lifestyle_cat", 
    covariates,
    "Cumulative 4-level (0/1/2/3+)",
    calc_trend = TRUE
  )
  
  if (!is.null(result)) {
    cumulative_results[[out_name]] <- result
    print(result %>% select(Variable, HR, Lower_CI, Upper_CI, P_value, P_trend))
  }
}

cumulative_results_df <- bind_rows(cumulative_results)

# IMMEDIATE SAVE: Save cumulative results right away
cat("\n--- Immediate Save: Primary Cumulative Results ---\n")
tryCatch({
  if (nrow(cumulative_results_df) > 0) {
    write.csv(cumulative_results_df, file.path(OUTPUT_DIR, "Cox_Primary_Cumulative_4Level.csv"), row.names = FALSE)
    cat("  [OK] Saved Cox_Primary_Cumulative_4Level.csv (", nrow(cumulative_results_df), " rows)\n")
  } else {
    cat("  [WARN] No cumulative results to save\n")
  }
}, error = function(e) {
  cat("  [ERROR] Failed to save cumulative results:", e$message, "\n")
})

# =============================================================================
# PART 8: SENSITIVITY ANALYSIS S1 - 5-level Categories
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   SENSITIVITY ANALYSIS S1: 5-level Categories (0/1/2/3/4)       \n")
cat("================================================================\n")

sens1_results <- list()

for (out_name in all_outcomes) {
  cat("\n--- Outcome:", outcome_info[[out_name]]$label, "---\n")
  
  result <- run_cox_model(
    pooled_analysis, 
    out_name, 
    "n_lifestyle_5cat", 
    covariates,
    "S1: 5-level (0/1/2/3/4)",
    calc_trend = TRUE
  )
  
  if (!is.null(result)) {
    sens1_results[[out_name]] <- result
    print(result %>% select(Variable, HR, Lower_CI, Upper_CI, P_value))
  }
}

sens1_results_df <- bind_rows(sens1_results)

# =============================================================================
# PART 9: SENSITIVITY ANALYSIS S2 - Heavy Drinking Definition
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   SENSITIVITY ANALYSIS S2: Heavy Drinking Definition            \n")
cat("   (Replace unhealthy_drink with heavy_drink)                    \n")
cat("================================================================\n")

# S2a: Individual factors with heavy drink
sens2a_results <- list()

for (out_name in all_outcomes) {
  cat("\n--- Individual Factors (Heavy Drink) - Outcome:", 
      outcome_info[[out_name]]$label, "---\n")
  
  result <- run_individual_factors(pooled_analysis, out_name, covariates, use_heavy_drink = TRUE)
  
  if (!is.null(result)) {
    sens2a_results[[out_name]] <- result
    print(result %>% select(Variable_Label, HR, Lower_CI, Upper_CI, P_value))
  }
}

sens2a_results_df <- bind_rows(sens2a_results)

# S2b: Cumulative effect with heavy drink score
sens2b_results <- list()

if ("n_lifestyle_cat_heavy" %in% names(pooled_analysis)) {
  for (out_name in all_outcomes) {
    cat("\n--- Cumulative (Heavy Drink Score) - Outcome:", 
        outcome_info[[out_name]]$label, "---\n")
    
    result <- run_cox_model(
      pooled_analysis, 
      out_name, 
      "n_lifestyle_cat_heavy", 
      covariates,
      "S2: Heavy Drink - 4-level",
      calc_trend = TRUE
    )
    
    if (!is.null(result)) {
      sens2b_results[[out_name]] <- result
      print(result %>% select(Variable, HR, Lower_CI, Upper_CI, P_value))
    }
  }
}

sens2b_results_df <- bind_rows(sens2b_results)

# =============================================================================
# PART 10: SENSITIVITY ANALYSIS S3 - MICE Imputed Data (with Rubin's Rules)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   SENSITIVITY ANALYSIS S3: MICE Imputed Data                    \n")
cat("   Using Rubin's Rules for Multiple Imputation Pooling           \n")
cat("================================================================\n")

sens3_individual <- list()
sens3_cumulative <- list()

# Check for MICE mids object (preferred for Rubin's rules)
mice_mids_path <- file.path(DATA_DIR, "Pooled", "Pooled_mice_mids.rds")
mice_csv_path <- file.path(DATA_DIR, "Pooled", "Pooled_mice_imputed.csv")

cat("\nChecking MICE data availability:\n")
cat("  mids object:", ifelse(file.exists(mice_mids_path), "FOUND", "Not found"), "\n")
cat("  CSV file:", ifelse(file.exists(mice_csv_path), "FOUND", "Not found"), "\n")

if (file.exists(mice_csv_path)) {
  cat("\n  Loading MICE imputed data...\n")
  
  pooled_mice <- read.csv(mice_csv_path, stringsAsFactors = FALSE)
  cat("  Loaded:", nrow(pooled_mice), "rows\n")
  
  # Check if outcome variables exist
  if ("event_ppcmm" %in% names(pooled_mice) && "time_ppcmm_months" %in% names(pooled_mice)) {
    
    # Prepare analysis data
    pooled_mice_analysis <- prepare_analysis_data(pooled_mice, use_heavy_drink = FALSE)
    cat("  MICE analysis sample size:", nrow(pooled_mice_analysis), "\n")
    
    if (nrow(pooled_mice_analysis) > 100) {
      
      # Individual factors
      cat("\n--- S3: Individual Factors (MICE) ---\n")
      for (out_name in all_outcomes) {
        cat("  Processing:", outcome_info[[out_name]]$label, "\n")
        result <- run_individual_factors(pooled_mice_analysis, out_name, covariates, use_heavy_drink = FALSE)
        if (!is.null(result)) {
          result$Model <- "S3: MICE - Individual Factors"
          sens3_individual[[out_name]] <- result
          print(result %>% select(Variable_Label, HR, Lower_CI, Upper_CI, P_value))
        }
      }
      
      # Cumulative effect
      cat("\n--- S3: Cumulative Effect (MICE) ---\n")
      for (out_name in all_outcomes) {
        cat("  Processing:", outcome_info[[out_name]]$label, "\n")
        result <- run_cox_model(
          pooled_mice_analysis, 
          out_name, 
          "n_lifestyle_cat", 
          covariates,
          "S3: MICE - Cumulative 4-level",
          calc_trend = TRUE
        )
        if (!is.null(result)) {
          sens3_cumulative[[out_name]] <- result
          print(result %>% select(Variable, HR, Lower_CI, Upper_CI, P_value))
        }
      }
    } else {
      cat("  [WARN] Insufficient sample size after filtering\n")
    }
    
  } else {
    cat("  [WARN] Outcome variables not found in MICE data\n")
    cat("  Available columns:", paste(head(names(pooled_mice), 20), collapse = ", "), "...\n")
  }
  
} else {
  cat("\n  [INFO] No MICE data available.\n")
  cat("  Run 02_MICE_Imputation.R first to generate MICE imputed data.\n")
  cat("  Skipping S3 sensitivity analysis.\n")
}

sens3_individual_df <- bind_rows(sens3_individual)
sens3_cumulative_df <- bind_rows(sens3_cumulative)

cat("\n  S3 Results Summary:\n")
cat("    Individual factors:", nrow(sens3_individual_df), "rows\n")
cat("    Cumulative effect:", nrow(sens3_cumulative_df), "rows\n")

# =============================================================================
# PART 11: SENSITIVITY ANALYSIS S4 - Drop First Follow-up (using main data)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   SENSITIVITY ANALYSIS S4: Drop First Follow-up Wave            \n")
cat("   (Using event_ppcmm_drop1st outcome)                           \n")
cat("================================================================\n")

# Create modified outcome_info for drop1st
outcome_info_drop1st <- list(
  Overall_drop1st = list(
    event = "event_ppcmm_drop1st", 
    time = "time_ppcmm_drop1st",
    label = "Overall PPC-MM (Drop1st)",
    type = "Sensitivity"
  )
)

sens4_results <- list()

# Check if drop1st variables exist
if ("event_ppcmm_drop1st" %in% names(pooled_analysis)) {
  
  df_drop1st <- pooled_analysis %>%
    filter(!is.na(event_ppcmm_drop1st) & !is.na(time_ppcmm_drop1st) & time_ppcmm_drop1st > 0)
  
  n_events_drop1st <- sum(df_drop1st$event_ppcmm_drop1st, na.rm = TRUE)
  cat("Drop1st events:", n_events_drop1st, "\n")
  
  if (n_events_drop1st >= 10) {
    # Individual factors
    formula_str <- paste0("Surv(time_ppcmm_drop1st, event_ppcmm_drop1st) ~ ",
                          "unhealthy_drink + unhealthy_smoke + unhealthy_pa + unhealthy_soc + ",
                          paste(covariates, collapse = " + "))
    
    fit <- tryCatch({
      coxph(as.formula(formula_str), data = df_drop1st)
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      s <- summary(fit)
      lifestyle_vars <- c("unhealthy_drink", "unhealthy_smoke", "unhealthy_pa", "unhealthy_soc")
      exp_rows <- rownames(s$coefficients) %in% lifestyle_vars
      
      sens4_individual <- data.frame(
        Model = "S4: Drop1st - Individual Factors",
        Outcome = "Overall_drop1st",
        Outcome_Label = "Overall PPC-MM (Drop1st)",
        Variable = rownames(s$coefficients)[exp_rows],
        HR = round(s$coefficients[exp_rows, "exp(coef)"], 3),
        Lower_CI = round(s$conf.int[exp_rows, "lower .95"], 3),
        Upper_CI = round(s$conf.int[exp_rows, "upper .95"], 3),
        P_value = format.pval(s$coefficients[exp_rows, "Pr(>|z|)"], digits = 3),
        N = nrow(df_drop1st),
        N_Events = n_events_drop1st,
        stringsAsFactors = FALSE
      )
      sens4_results[["individual"]] <- sens4_individual
      cat("\n--- Individual Factors (Drop1st) ---\n")
      print(sens4_individual %>% select(Variable, HR, Lower_CI, Upper_CI, P_value))
    }
    
    # Cumulative effect
    formula_str2 <- paste0("Surv(time_ppcmm_drop1st, event_ppcmm_drop1st) ~ ",
                           "n_lifestyle_cat + ", paste(covariates, collapse = " + "))
    
    fit2 <- tryCatch({
      coxph(as.formula(formula_str2), data = df_drop1st)
    }, error = function(e) NULL)
    
    if (!is.null(fit2)) {
      s2 <- summary(fit2)
      exp_rows2 <- grepl("n_lifestyle_cat", rownames(s2$coefficients))
      
      # P for trend
      fit_trend <- coxph(as.formula(paste0("Surv(time_ppcmm_drop1st, event_ppcmm_drop1st) ~ ",
                                           "n_lifestyle_continuous + ", paste(covariates, collapse = " + "))),
                         data = df_drop1st)
      p_trend <- summary(fit_trend)$coefficients["n_lifestyle_continuous", "Pr(>|z|)"]
      
      sens4_cumulative <- data.frame(
        Model = "S4: Drop1st - Cumulative 4-level",
        Outcome = "Overall_drop1st",
        Outcome_Label = "Overall PPC-MM (Drop1st)",
        Variable = rownames(s2$coefficients)[exp_rows2],
        HR = round(s2$coefficients[exp_rows2, "exp(coef)"], 3),
        Lower_CI = round(s2$conf.int[exp_rows2, "lower .95"], 3),
        Upper_CI = round(s2$conf.int[exp_rows2, "upper .95"], 3),
        P_value = format.pval(s2$coefficients[exp_rows2, "Pr(>|z|)"], digits = 3),
        P_trend = format.pval(p_trend, digits = 3),
        N = nrow(df_drop1st),
        N_Events = n_events_drop1st,
        stringsAsFactors = FALSE
      )
      sens4_results[["cumulative"]] <- sens4_cumulative
      cat("\n--- Cumulative Effect (Drop1st) ---\n")
      print(sens4_cumulative %>% select(Variable, HR, Lower_CI, Upper_CI, P_value, P_trend))
    }
  }
} else {
  cat("Drop1st variables not available, skipping S4\n")
}

sens4_results_df <- bind_rows(sens4_results)

# =============================================================================
# PART 12: SUBGROUP ANALYSIS - Age Stratification
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   SUBGROUP: Age-Stratified Cox Regression                       \n")
cat("================================================================\n")

# Verify data availability
cat("\n>>> Checking pooled_analysis data...\n")
if (!exists("pooled_analysis") || is.null(pooled_analysis) || nrow(pooled_analysis) == 0) {
  cat("  [ERROR] pooled_analysis not available. Trying to recreate...\n")
  if (exists("pooled_main") && nrow(pooled_main) > 0) {
    pooled_analysis <- prepare_analysis_data(pooled_main, use_heavy_drink = TRUE)
    cat("  [OK] Recreated pooled_analysis with", nrow(pooled_analysis), "rows\n")
  } else {
    cat("  [FATAL] Cannot create pooled_analysis. Skipping age stratification.\n")
    age_stratified_df <- data.frame(Note = "Data not available")
    write.csv(age_stratified_df, file.path(OUTPUT_DIR, "Cox_Age_Stratified.csv"), row.names = FALSE)
  }
}

# Create age groups
cat(">>> Creating age groups...\n")
pooled_analysis <- pooled_analysis %>%
  mutate(
    age_group_4 = case_when(
      age_baseline >= 50 & age_baseline < 60 ~ "50-59",
      age_baseline >= 60 & age_baseline < 70 ~ "60-69",
      age_baseline >= 70 & age_baseline < 80 ~ "70-79",
      age_baseline >= 80 ~ "80+",
      TRUE ~ NA_character_
    )
  )

cat("\nAge group distribution:\n")
age_dist <- pooled_analysis %>%
  group_by(age_group_4) %>%
  summarise(N = n(), Events = sum(event_ppcmm, na.rm = TRUE), .groups = "drop")
print(age_dist)

# Run age-stratified analysis
age_groups <- c("50-59", "60-69", "70-79", "80+")
age_stratified_results <- list()

for (ag in age_groups) {
  cat("\n--- Age Group:", ag, "---\n")
  
  df_age <- pooled_analysis %>% filter(age_group_4 == ag)
  n_total <- nrow(df_age)
  n_events <- sum(df_age$event_ppcmm, na.rm = TRUE)
  
  cat("  N =", n_total, ", Events =", n_events, "\n")
  
  if (n_events < 10) {
    cat("  [SKIP] Insufficient events\n")
    next
  }
  
  # Local variable definitions
  lf_vars <- c("unhealthy_drink", "unhealthy_smoke", "unhealthy_pa", "unhealthy_soc")
  lf_labels <- c("Drinking", "Smoking", "Physical Inactivity", "Social Isolation")
  covs <- c("sex", "edu", "region")
  covs <- covs[covs %in% names(df_age)]
  
  # Build formula
  n_cohorts <- length(unique(df_age$cohort))
  if (n_cohorts > 1) {
    formula_str <- paste0("Surv(time_ppcmm_months, event_ppcmm) ~ ",
                          paste(lf_vars, collapse = " + "), " + ",
                          paste(covs, collapse = " + "), " + strata(cohort)")
  } else {
    formula_str <- paste0("Surv(time_ppcmm_months, event_ppcmm) ~ ",
                          paste(lf_vars, collapse = " + "), " + ",
                          paste(covs, collapse = " + "))
  }
  
  fit <- tryCatch({
    coxph(as.formula(formula_str), data = df_age)
  }, error = function(e) {
    cat("  [ERROR]", e$message, "\n")
    NULL
  })
  
  if (!is.null(fit)) {
    s <- summary(fit)
    for (i in seq_along(lf_vars)) {
      lv <- lf_vars[i]
      if (!lv %in% rownames(s$coefficients)) next
      
      age_stratified_results[[paste0(ag, "_", lv)]] <- data.frame(
        Analysis = "Age-Stratified",
        Age_Group = ag,
        Model = "Individual Factors",
        Outcome = "Overall PPC-MM",
        Variable = lv,
        Variable_Label = lf_labels[i],
        Prevalence = round(mean(df_age[[lv]], na.rm = TRUE), 4),
        HR = round(s$coefficients[lv, "exp(coef)"], 3),
        Lower_CI = round(s$conf.int[lv, "lower .95"], 3),
        Upper_CI = round(s$conf.int[lv, "upper .95"], 3),
        HR_CI = sprintf("%.2f (%.2f-%.2f)", 
                        s$coefficients[lv, "exp(coef)"],
                        s$conf.int[lv, "lower .95"],
                        s$conf.int[lv, "upper .95"]),
        P_value = s$coefficients[lv, "Pr(>|z|)"],
        P_value_fmt = format.pval(s$coefficients[lv, "Pr(>|z|)"], digits = 4),
        N = n_total,
        N_Events = n_events,
        stringsAsFactors = FALSE
      )
    }
    cat("  [OK] Model fitted\n")
  }
}

# Save age-stratified results
age_stratified_df <- if (length(age_stratified_results) > 0) {
  bind_rows(age_stratified_results)
} else {
  data.frame(Note = "No results generated")
}

cat("\n>>> Saving Age-Stratified Results:", nrow(age_stratified_df), "rows\n")
if (nrow(age_stratified_df) > 0 && !"Note" %in% names(age_stratified_df)) {
  print(age_stratified_df %>% select(Age_Group, Variable_Label, HR_CI, P_value_fmt))
}
write.csv(age_stratified_df, file.path(OUTPUT_DIR, "Cox_Age_Stratified.csv"), row.names = FALSE)
cat("[OK] Saved: Cox_Age_Stratified.csv\n")

# =============================================================================
# PART 13: SUBGROUP ANALYSIS - Sex Stratification
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   SUBGROUP: Sex-Stratified Cox Regression                       \n")
cat("================================================================\n")

cat("\nSex distribution:\n")
sex_dist <- pooled_analysis %>%
  group_by(sex) %>%
  summarise(N = n(), Events = sum(event_ppcmm, na.rm = TRUE), .groups = "drop")
print(sex_dist)

# Run sex-stratified analysis
sex_groups <- c("Men", "Women")
sex_stratified_results <- list()

for (sg in sex_groups) {
  cat("\n--- Sex:", sg, "---\n")
  
  df_sex <- pooled_analysis %>% filter(sex == sg)
  n_total <- nrow(df_sex)
  n_events <- sum(df_sex$event_ppcmm, na.rm = TRUE)
  
  cat("  N =", n_total, ", Events =", n_events, "\n")
  
  if (n_events < 10) {
    cat("  [SKIP] Insufficient events\n")
    next
  }
  
  # Local variable definitions
  lf_vars <- c("unhealthy_drink", "unhealthy_smoke", "unhealthy_pa", "unhealthy_soc")
  lf_labels <- c("Drinking", "Smoking", "Physical Inactivity", "Social Isolation")
  covs <- c("age_baseline", "edu", "region")
  covs <- covs[covs %in% names(df_sex)]
  
  # Build formula
  n_cohorts <- length(unique(df_sex$cohort))
  if (n_cohorts > 1) {
    formula_str <- paste0("Surv(time_ppcmm_months, event_ppcmm) ~ ",
                          paste(lf_vars, collapse = " + "), " + ",
                          paste(covs, collapse = " + "), " + strata(cohort)")
  } else {
    formula_str <- paste0("Surv(time_ppcmm_months, event_ppcmm) ~ ",
                          paste(lf_vars, collapse = " + "), " + ",
                          paste(covs, collapse = " + "))
  }
  
  fit <- tryCatch({
    coxph(as.formula(formula_str), data = df_sex)
  }, error = function(e) {
    cat("  [ERROR]", e$message, "\n")
    NULL
  })
  
  if (!is.null(fit)) {
    s <- summary(fit)
    for (i in seq_along(lf_vars)) {
      lv <- lf_vars[i]
      if (!lv %in% rownames(s$coefficients)) next
      
      sex_stratified_results[[paste0(sg, "_", lv)]] <- data.frame(
        Analysis = "Sex-Stratified",
        Sex = sg,
        Model = "Individual Factors",
        Outcome = "Overall PPC-MM",
        Variable = lv,
        Variable_Label = lf_labels[i],
        Prevalence = round(mean(df_sex[[lv]], na.rm = TRUE), 4),
        HR = round(s$coefficients[lv, "exp(coef)"], 3),
        Lower_CI = round(s$conf.int[lv, "lower .95"], 3),
        Upper_CI = round(s$conf.int[lv, "upper .95"], 3),
        HR_CI = sprintf("%.2f (%.2f-%.2f)", 
                        s$coefficients[lv, "exp(coef)"],
                        s$conf.int[lv, "lower .95"],
                        s$conf.int[lv, "upper .95"]),
        P_value = s$coefficients[lv, "Pr(>|z|)"],
        P_value_fmt = format.pval(s$coefficients[lv, "Pr(>|z|)"], digits = 4),
        N = n_total,
        N_Events = n_events,
        stringsAsFactors = FALSE
      )
    }
    cat("  [OK] Model fitted\n")
  }
}

# Save sex-stratified results
sex_stratified_df <- if (length(sex_stratified_results) > 0) {
  bind_rows(sex_stratified_results)
} else {
  data.frame(Note = "No results generated")
}

cat("\n>>> Saving Sex-Stratified Results:", nrow(sex_stratified_df), "rows\n")
if (nrow(sex_stratified_df) > 0 && !"Note" %in% names(sex_stratified_df)) {
  print(sex_stratified_df %>% select(Sex, Variable_Label, HR_CI, P_value_fmt))
}
write.csv(sex_stratified_df, file.path(OUTPUT_DIR, "Cox_Sex_Stratified.csv"), row.names = FALSE)
cat("[OK] Saved: Cox_Sex_Stratified.csv\n")

# =============================================================================
# PART 14: INTERACTION ANALYSIS
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   INTERACTION ANALYSIS: Lifestyle × Sex                         \n")
cat("================================================================\n")

# Local variable definitions
lf_vars <- c("unhealthy_drink", "unhealthy_smoke", "unhealthy_pa", "unhealthy_soc")
lf_labels <- c("Drinking", "Smoking", "Physical Inactivity", "Social Isolation")

interaction_results <- list()

cat("\n--- Testing Lifestyle × Sex Interactions ---\n")
for (i in seq_along(lf_vars)) {
  lv <- lf_vars[i]
  
  # Create interaction term
  pooled_analysis$int_term <- pooled_analysis[[lv]] * as.numeric(pooled_analysis$sex == "Women")
  
  formula_int <- paste0("Surv(time_ppcmm_months, event_ppcmm) ~ ",
                        lv, " + sex + int_term + age_baseline + edu + region + strata(cohort)")
  
  fit_int <- tryCatch({
    coxph(as.formula(formula_int), data = pooled_analysis)
  }, error = function(e) NULL)
  
  if (!is.null(fit_int) && "int_term" %in% rownames(summary(fit_int)$coefficients)) {
    s <- summary(fit_int)
    p_int <- s$coefficients["int_term", "Pr(>|z|)"]
    hr_int <- s$coefficients["int_term", "exp(coef)"]
    
    interaction_results[[paste0("sex_", lv)]] <- data.frame(
      Interaction = paste0(lf_labels[i], " × Sex"),
      Variable = lv,
      HR_Interaction = round(hr_int, 3),
      P_Interaction = p_int,
      P_fmt = format.pval(p_int, digits = 4),
      Significant = ifelse(p_int < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    )
    cat("  ", lf_labels[i], "× Sex: P =", format.pval(p_int, digits = 4), "\n")
  }
}

# Save interaction results
interaction_df <- if (length(interaction_results) > 0) {
  bind_rows(interaction_results)
} else {
  data.frame(Note = "No results generated")
}

cat("\n>>> Saving Interaction Results:", nrow(interaction_df), "rows\n")
if (nrow(interaction_df) > 0 && !"Note" %in% names(interaction_df)) {
  print(interaction_df)
}
write.csv(interaction_df, file.path(OUTPUT_DIR, "Cox_Interaction_Analysis.csv"), row.names = FALSE)
cat("[OK] Saved: Cox_Interaction_Analysis.csv\n")

# =============================================================================
# PART 15: Dose-Response Plots
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   DOSE-RESPONSE PLOTS                                           \n")
cat("================================================================\n")

# Create dose-response data
create_dose_response_data <- function(results_df, exposure_levels) {
  
  # Add reference (0 = 1.00)
  ref_row <- data.frame(
    n_lifestyle = 0,
    HR = 1.00,
    Lower_CI = 1.00,
    Upper_CI = 1.00
  )
  
  # Extract HRs for each level
  dr_data <- results_df %>%
    mutate(
      n_lifestyle = as.numeric(gsub("n_lifestyle_cat|n_lifestyle_5cat", "", Variable))
    ) %>%
    select(n_lifestyle, HR, Lower_CI, Upper_CI)
  
  # Combine with reference
  dr_data <- bind_rows(ref_row, dr_data) %>%
    arrange(n_lifestyle)
  
  return(dr_data)
}

# Plot function
plot_dose_response <- function(dr_data, outcome_label, p_trend = NULL, filename) {
  
  p <- ggplot(dr_data, aes(x = n_lifestyle, y = HR)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    geom_line(color = "#2E86AB", size = 1) +
    geom_point(size = 3, color = "#2E86AB") +
    geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), 
                  width = 0.1, color = "#2E86AB", size = 0.8) +
    scale_x_continuous(breaks = dr_data$n_lifestyle) +
    scale_y_log10(limits = c(0.8, max(dr_data$Upper_CI) * 1.1)) +
    labs(
      x = "Number of Unhealthy Lifestyle Factors",
      y = "Hazard Ratio (95% CI)",
      title = outcome_label,
      subtitle = if (!is.null(p_trend)) paste("P for trend =", p_trend) else NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(OUTPUT_DIR, filename), p, width = 6, height = 5, dpi = 300)
  cat("Saved:", filename, "\n")
  
  return(p)
}

# Generate dose-response plots for main analysis (4-level)
cat("\n--- Generating Dose-Response Plots (4-level) ---\n")

dr_plots_4level <- list()

for (out_name in all_outcomes) {
  result <- cumulative_results[[out_name]]
  if (!is.null(result) && nrow(result) > 0) {
    dr_data <- create_dose_response_data(result, c("1", "2", "3+"))
    dr_data$n_lifestyle[dr_data$n_lifestyle == 3] <- 3  # Keep as 3 for 3+
    
    p_trend <- result$P_trend[1]
    filename <- paste0("DoseResponse_4level_", out_name, ".png")
    
    dr_plots_4level[[out_name]] <- plot_dose_response(
      dr_data, 
      outcome_info[[out_name]]$label, 
      p_trend, 
      filename
    )
  }
}

# Generate dose-response plots for sensitivity (5-level)
cat("\n--- Generating Dose-Response Plots (5-level) ---\n")

dr_plots_5level <- list()

for (out_name in all_outcomes) {
  result <- sens1_results[[out_name]]
  if (!is.null(result) && nrow(result) > 0) {
    dr_data <- create_dose_response_data(result, c("1", "2", "3", "4"))
    
    p_trend <- result$P_trend[1]
    filename <- paste0("DoseResponse_5level_", out_name, ".png")
    
    dr_plots_5level[[out_name]] <- plot_dose_response(
      dr_data, 
      paste(outcome_info[[out_name]]$label, "(5-level)"), 
      p_trend, 
      filename
    )
  }
}

# =============================================================================
# PART 13: Combined Dose-Response Figure (Panel Plot)
# =============================================================================

cat("\n--- Creating Combined Dose-Response Panel ---\n")

# Create combined data for all outcomes
combined_dr_data <- data.frame()

for (out_name in all_outcomes) {
  result <- cumulative_results[[out_name]]
  if (!is.null(result) && nrow(result) > 0) {
    dr_data <- create_dose_response_data(result, c("1", "2", "3+"))
    dr_data$Outcome <- outcome_info[[out_name]]$label
    dr_data$Outcome_Type <- outcome_info[[out_name]]$type
    combined_dr_data <- bind_rows(combined_dr_data, dr_data)
  }
}

if (nrow(combined_dr_data) > 0) {
  # Reorder outcomes
  combined_dr_data$Outcome <- factor(combined_dr_data$Outcome, 
                                     levels = c("Overall PPC-MM", 
                                                "Physical-Psychological (P1P2)",
                                                "Physical-Cognitive (P1C)",
                                                "Psychological-Cognitive (P2C)",
                                                "All Three (P1P2C)"))
  
  p_combined <- ggplot(combined_dr_data, aes(x = n_lifestyle, y = HR, color = Outcome)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    geom_line(size = 1) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.1, size = 0.6) +
    scale_x_continuous(breaks = 0:3, labels = c("0", "1", "2", "3+")) +
    scale_y_log10() +
    scale_color_manual(values = c(
      "Overall PPC-MM" = "#E63946",
      "Physical-Psychological (P1P2)" = "#457B9D",
      "Physical-Cognitive (P1C)" = "#2A9D8F",
      "Psychological-Cognitive (P2C)" = "#E9C46A",
      "All Three (P1P2C)" = "#9B2335"
    )) +
    labs(
      x = "Number of Unhealthy Lifestyle Factors",
      y = "Hazard Ratio (95% CI)",
      title = "Dose-Response Relationship: Lifestyle Factors and PPC-MM",
      color = "Outcome"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(nrow = 2))
  
  ggsave(file.path(OUTPUT_DIR, "DoseResponse_Combined_AllOutcomes.png"), 
         p_combined, width = 10, height = 7, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, "DoseResponse_Combined_AllOutcomes.pdf"), 
         p_combined, width = 10, height = 7)
  cat("Saved: DoseResponse_Combined_AllOutcomes.png/pdf\n")
}

# =============================================================================
# PART 14: Combine and Save All Results
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   SAVING RESULTS                                                \n")
cat("================================================================\n")

# Combine all results (including subgroup analyses)
all_results <- bind_rows(
  individual_results_df %>% mutate(Analysis = "Primary: Individual Factors"),
  cumulative_results_df %>% mutate(Analysis = "Primary: Cumulative (4-level)"),
  sens1_results_df %>% mutate(Analysis = "S1: 5-level Categories"),
  sens2a_results_df %>% mutate(Analysis = "S2a: Heavy Drink - Individual"),
  sens2b_results_df %>% mutate(Analysis = "S2b: Heavy Drink - Cumulative"),
  sens3_individual_df %>% mutate(Analysis = "S3: MICE - Individual"),
  sens3_cumulative_df %>% mutate(Analysis = "S3: MICE - Cumulative"),
  sens4_results_df %>% mutate(Analysis = "S4: Drop1st"),
  # Include subgroup analyses in all_results
  if(exists("age_stratified_df") && nrow(age_stratified_df) > 0) age_stratified_df else NULL,
  if(exists("sex_stratified_df") && nrow(sex_stratified_df) > 0) sex_stratified_df else NULL,
  if(exists("interaction_df") && nrow(interaction_df) > 0) interaction_df else NULL
)

# Create analysis description
analysis_description <- data.frame(
  Analysis = c(
    "Primary: Individual Factors",
    "Primary: Cumulative (4-level)",
    "S1: 5-level Categories",
    "S2a: Heavy Drink - Individual",
    "S2b: Heavy Drink - Cumulative",
    "S3: MICE - Individual",
    "S3: MICE - Cumulative",
    "S4: Drop1st"
  ),
  Description = c(
    "Individual lifestyle factors (drink, smoke, PA, social) mutually adjusted; Outcomes: Overall (Primary) + 4 subtypes (Secondary)",
    "Cumulative lifestyle categories (0/1/2/3+) with P for trend; Outcomes: Overall (Primary) + 4 subtypes (Secondary)",
    "Sensitivity: 5-level categories (0/1/2/3/4) for finer dose-response assessment",
    "Sensitivity: Replace any drinking with heavy drinking definition for individual factors",
    "Sensitivity: Replace any drinking with heavy drinking for cumulative score",
    "Sensitivity: MICE imputed data for individual factors",
    "Sensitivity: MICE imputed data for cumulative effect",
    "Sensitivity: Exclude first follow-up wave to address reverse causality"
  ),
  Reference = c(
    "Healthy behavior for each factor",
    "0 unhealthy lifestyles",
    "0 unhealthy lifestyles",
    "Healthy behavior for each factor",
    "0 unhealthy lifestyles (heavy drink definition)",
    "Healthy behavior for each factor",
    "0 unhealthy lifestyles",
    "0 unhealthy lifestyles"
  ),
  stringsAsFactors = FALSE
)

# Prepare Excel output (including subgroup and interaction analyses)
results_for_excel <- list(
  "Analysis_Description" = analysis_description,
  "Primary_Individual" = individual_results_df,
  "Primary_Cumulative" = cumulative_results_df,
  "S1_5Level" = sens1_results_df,
  "S2a_HeavyDrink_Individual" = sens2a_results_df,
  "S2b_HeavyDrink_Cumulative" = sens2b_results_df,
  "S3_MICE_Individual" = sens3_individual_df,
  "S3_MICE_Cumulative" = sens3_cumulative_df,
  "S4_Drop1st" = sens4_results_df,
  "Age_Stratified" = if(exists("age_stratified_df") && nrow(age_stratified_df) > 0) age_stratified_df else data.frame(Note = "No data"),
  "Sex_Stratified" = if(exists("sex_stratified_df") && nrow(sex_stratified_df) > 0) sex_stratified_df else data.frame(Note = "No data"),
  "Interaction_Analysis" = if(exists("interaction_df") && nrow(interaction_df) > 0) interaction_df else data.frame(Note = "No data"),
  "All_Results" = all_results
)

# Save Excel and main CSV with robust error handling
cat("\n--- Saving Excel and CSV Files ---\n")

# Save Excel with tryCatch
tryCatch({
  cat("  Saving Excel file...\n")
  filepath_xlsx <- file.path(OUTPUT_DIR, "Pooled_Cox_Results_Comprehensive.xlsx")
  writexl::write_xlsx(results_for_excel, filepath_xlsx)
  cat("  [OK] Saved: Pooled_Cox_Results_Comprehensive.xlsx\n")
}, error = function(e) {
  cat("  [ERROR] Excel save failed:", e$message, "\n")
})

# Save main CSV
tryCatch({
  cat("  Saving main CSV file...\n")
  if (nrow(all_results) > 0) {
    filepath_csv <- file.path(OUTPUT_DIR, "Pooled_Cox_All_Results.csv")
    write.csv(all_results, filepath_csv, row.names = FALSE)
    cat("  [OK] Saved: Pooled_Cox_All_Results.csv (", nrow(all_results), " rows)\n")
  } else {
    cat("  [WARN] all_results is empty, skipping main CSV\n")
  }
}, error = function(e) {
  cat("  [ERROR] Main CSV save failed:", e$message, "\n")
})

# ===== Save Individual CSV Files for Each Analysis =====
cat("\n--- Saving Individual CSV Files (with error handling) ---\n")

# Helper function for safe CSV saving
safe_save_csv <- function(df, filename, df_name) {
  tryCatch({
    if (!is.null(df) && nrow(df) > 0) {
      filepath <- file.path(OUTPUT_DIR, filename)
      write.csv(df, filepath, row.names = FALSE)
      cat("  [OK]", filename, "(", nrow(df), "rows)\n")
    } else {
      cat("  [SKIP]", filename, "- empty or NULL data\n")
    }
  }, error = function(e) {
    cat("  [ERROR]", filename, "-", e$message, "\n")
  })
}

# 1. Primary Individual Factors
safe_save_csv(individual_results_df, "Cox_Primary_Individual_Factors.csv", "Primary Individual")

# 2. Primary Cumulative Effect (4-level)
safe_save_csv(cumulative_results_df, "Cox_Primary_Cumulative_4Level.csv", "Primary Cumulative")

# 3. Sensitivity 1: 5-Level
safe_save_csv(sens1_results_df, "Cox_S1_5Level_Categories.csv", "S1 5-Level")

# 4. Sensitivity 2a: Heavy Drink Individual
safe_save_csv(sens2a_results_df, "Cox_S2a_HeavyDrink_Individual.csv", "S2a Heavy Individual")

# 5. Sensitivity 2b: Heavy Drink Cumulative
safe_save_csv(sens2b_results_df, "Cox_S2b_HeavyDrink_Cumulative.csv", "S2b Heavy Cumulative")

# 6. Sensitivity 3 MICE Individual
safe_save_csv(sens3_individual_df, "Cox_S3_MICE_Individual.csv", "S3 MICE Individual")

# 7. Sensitivity 3 MICE Cumulative
safe_save_csv(sens3_cumulative_df, "Cox_S3_MICE_Cumulative.csv", "S3 MICE Cumulative")

# 8. Sensitivity 4: Drop1st
safe_save_csv(sens4_results_df, "Cox_S4_Drop1st.csv", "S4 Drop1st")

# 9. Age-Stratified Analysis (IMPORTANT: must save)
if (exists("age_stratified_df")) {
  safe_save_csv(age_stratified_df, "Cox_Age_Stratified.csv", "Age Stratified")
} else {
  cat("  [WARN] age_stratified_df not found\n")
}

# 10. Sex-Stratified Analysis (IMPORTANT: must save)
if (exists("sex_stratified_df")) {
  safe_save_csv(sex_stratified_df, "Cox_Sex_Stratified.csv", "Sex Stratified")
} else {
  cat("  [WARN] sex_stratified_df not found\n")
}

# 11. Interaction Analysis (IMPORTANT: must save)
if (exists("interaction_df")) {
  safe_save_csv(interaction_df, "Cox_Interaction_Analysis.csv", "Interaction Analysis")
} else {
  cat("  [WARN] interaction_df not found\n")
}

# ===== Create Summary Tables =====
cat("\n--- Creating Summary Tables for Publication ---\n")

# Summary table: HR (95% CI) format for publication
create_hr_summary <- function(results_df, analysis_name) {
  if (is.null(results_df) || nrow(results_df) == 0) return(NULL)
  
  results_df %>%
    mutate(
      HR_CI = paste0(HR, " (", Lower_CI, "-", Upper_CI, ")"),
      Analysis_Name = analysis_name
    ) %>%
    select(Analysis_Name, Outcome, Variable, HR_CI, P_value, N, N_Events)
}

hr_summary_list <- list()

if (!is.null(individual_results_df) && nrow(individual_results_df) > 0) {
  hr_summary_list[["Individual"]] <- create_hr_summary(individual_results_df, "Primary: Individual Factors")
}
if (!is.null(cumulative_results_df) && nrow(cumulative_results_df) > 0) {
  hr_summary_list[["Cumulative"]] <- create_hr_summary(cumulative_results_df, "Primary: Cumulative (4-level)")
}
if (!is.null(sens1_results_df) && nrow(sens1_results_df) > 0) {
  hr_summary_list[["S1"]] <- create_hr_summary(sens1_results_df, "S1: 5-Level")
}
if (!is.null(sens2a_results_df) && nrow(sens2a_results_df) > 0) {
  hr_summary_list[["S2a"]] <- create_hr_summary(sens2a_results_df, "S2a: Heavy Drink Individual")
}
if (!is.null(sens2b_results_df) && nrow(sens2b_results_df) > 0) {
  hr_summary_list[["S2b"]] <- create_hr_summary(sens2b_results_df, "S2b: Heavy Drink Cumulative")
}
if (!is.null(sens3_individual_df) && nrow(sens3_individual_df) > 0) {
  hr_summary_list[["S3_Ind"]] <- create_hr_summary(sens3_individual_df, "S3: MICE Individual")
}
if (!is.null(sens3_cumulative_df) && nrow(sens3_cumulative_df) > 0) {
  hr_summary_list[["S3_Cum"]] <- create_hr_summary(sens3_cumulative_df, "S3: MICE Cumulative")
}
if (!is.null(sens4_results_df) && nrow(sens4_results_df) > 0) {
  hr_summary_list[["S4"]] <- create_hr_summary(sens4_results_df, "S4: Drop1st")
}

hr_summary_df <- bind_rows(hr_summary_list)

tryCatch({
  if (nrow(hr_summary_df) > 0) {
    filepath <- file.path(OUTPUT_DIR, "Cox_Summary_HR_CI.csv")
    write.csv(hr_summary_df, filepath, row.names = FALSE)
    cat("  [OK] Saved: Cox_Summary_HR_CI.csv (", nrow(hr_summary_df), " rows)\n")
  } else {
    cat("  [WARN] No summary data to save\n")
  }
}, error = function(e) {
  cat("  [ERROR] Failed to save summary:", e$message, "\n")
})

# Save publication-ready table (formatted for journal submission)
tryCatch({
  pub_table <- hr_summary_df %>%
    filter(!is.na(HR_CI)) %>%
    mutate(
      Result = paste0(HR_CI, ", P=", P_value),
      Sample_Info = paste0("N=", N, ", Events=", N_Events)
    ) %>%
    select(Analysis_Name, Outcome, Variable, Result, Sample_Info)
  
  if (nrow(pub_table) > 0) {
    filepath <- file.path(OUTPUT_DIR, "Cox_Publication_Table.csv")
    write.csv(pub_table, filepath, row.names = FALSE)
    cat("  [OK] Saved: Cox_Publication_Table.csv (", nrow(pub_table), " rows)\n")
  }
}, error = function(e) {
  cat("  [ERROR] Failed to save publication table:", e$message, "\n")
})

# ===== Save Model Objects =====
cox_results <- list(
  primary_individual = individual_results_df,
  primary_cumulative = cumulative_results_df,
  sens1_5level = sens1_results_df,
  sens2a_heavy_individual = sens2a_results_df,
  sens2b_heavy_cumulative = sens2b_results_df,
  sens3_mice_individual = sens3_individual_df,
  sens3_mice_cumulative = sens3_cumulative_df,
  sens4_drop1st = sens4_results_df,
  all_results = all_results,
  pooled_data = pooled_analysis
)

# Save RDS with tryCatch
tryCatch({
  saveRDS(cox_results, file.path(OUTPUT_DIR, "Pooled_cox_results.rds"))
  cat("  [OK] Saved: Pooled_cox_results.rds\n")
}, error = function(e) {
  cat("  [ERROR] Failed to save Cox RDS:", e$message, "\n")
})

# =============================================================================
# PART 15: Summary Output
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   ANALYSIS SUMMARY                                              \n")
cat("================================================================\n")

cat("\n--- PRIMARY ANALYSIS: Individual Factors ---\n")
if (nrow(individual_results_df) > 0) {
  summary_ind <- individual_results_df %>%
    group_by(Outcome, Outcome_Type) %>%
    summarise(
      N_Events = first(N_Events),
      Significant = sum(as.numeric(gsub("<", "", P_value)) < 0.05, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_ind)
}

cat("\n--- PRIMARY ANALYSIS: Cumulative Effect (4-level) ---\n")
if (nrow(cumulative_results_df) > 0) {
  summary_cum <- cumulative_results_df %>%
    group_by(Outcome, Outcome_Type) %>%
    summarise(
      N_Events = first(N_Events),
      HR_3plus = HR[Variable == "n_lifestyle_cat3+"],
      P_trend = first(P_trend),
      .groups = "drop"
    )
  print(summary_cum)
}

cat("\n--- SENSITIVITY ANALYSES ---\n")
cat("S1 (5-level): ", nrow(sens1_results_df), " results\n")
cat("S2a (Heavy Drink - Individual): ", nrow(sens2a_results_df), " results\n")
cat("S2b (Heavy Drink - Cumulative): ", nrow(sens2b_results_df), " results\n")
cat("S3 (MICE - Individual): ", nrow(sens3_individual_df), " results\n")
cat("S3 (MICE - Cumulative): ", nrow(sens3_cumulative_df), " results\n")
cat("S4 (Drop1st): ", nrow(sens4_results_df), " results\n")

# =============================================================================
# COMPLETE
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   Pooled Cox Regression Analysis Complete                       \n")
cat("================================================================\n")
cat("Sample size:", nrow(pooled_analysis), "\n")
cat("Number of cohorts:", length(unique(pooled_analysis$cohort)), "\n")
cat("Total result entries:", nrow(all_results), "\n")

# List all generated files
cat("\n--- Generated Output Files ---\n")
cat("Directory:", OUTPUT_DIR, "\n")

# Check which files were actually created
cox_files <- c(
  "Pooled_Cox_Results_Comprehensive.xlsx",
  "Pooled_Cox_All_Results.csv",
  "Pooled_cox_results.rds",
  "Cox_Primary_Individual_Factors.csv",
  "Cox_Primary_Cumulative_4Level.csv",
  "Cox_S1_5Level_Categories.csv",
  "Cox_S2a_HeavyDrink_Individual.csv",
  "Cox_S2b_HeavyDrink_Cumulative.csv",
  "Cox_S3_MICE_Individual.csv",
  "Cox_S3_MICE_Cumulative.csv",
  "Cox_S4_Drop1st.csv",
  "Cox_Summary_HR_CI.csv",
  "Cox_Publication_Table.csv"
)

for (f in cox_files) {
  filepath <- file.path(OUTPUT_DIR, f)
  if (file.exists(filepath)) {
    cat("  [✓]", f, "\n")
  } else {
    cat("  [✗]", f, "(NOT CREATED)\n")
  }
}

cat("\n================================================================\n")
cat("  ★★★ 04_Pooled_Cox_Analysis.R EXECUTION COMPLETE ★★★           \n")
cat("================================================================\n")



