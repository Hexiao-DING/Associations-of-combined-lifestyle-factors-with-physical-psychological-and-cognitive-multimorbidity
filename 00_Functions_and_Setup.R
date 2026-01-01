###############################################################################
# 00_Functions_and_Setup.R
# ============================================================================
# Project: Lifestyle Factors and PPC-MM (Physical-Psychological-Cognitive 
#          Multi-Morbidity) Association Study
# Purpose: Common functions and environment setup
# Author: [Your Name]
# Date: 2024
# ============================================================================
# This script contains:
#   1. Required R packages loading
#   2. Global path settings
#   3. Common helper functions
#   4. Variable recoding functions
#   5. Descriptive statistics functions
#   6. Cox regression standardized functions
###############################################################################

# =============================================================================
# PART 1: Environment Setup and Package Loading
# =============================================================================

# Note: rm(list = ls()) removed to avoid clearing workspace when sourced
# from master script (09_Run_All_Analysis.R)

# Load required R packages
required_packages <- c(
  # Data manipulation
  "tidyverse",    # Data wrangling (dplyr, tidyr, ggplot2, etc.)
  "data.table",   # Efficient large dataset handling
  
  # Survival analysis
  "survival",     # Cox proportional hazards model
  "survminer",    # Survival curve visualization
  
  # Multilevel models and ICC
  "lme4",         # Mixed effects models
  "performance",  # ICC calculation
  "rptR",         # Repeatability/ICC with CI and P-values
  
  # Parallel computing (for faster bootstrap)
  "parallel",     # Base R parallel computing
  "foreach",      # Parallel foreach loops
  "doParallel",   # Parallel backend for foreach
  
  # Correlation analysis
  "psych",        # Phi coefficient calculation
  "corrplot",     # Correlation matrix visualization
  
  # Table output
  "tableone",     # Baseline characteristics table
  "gtsummary",    # Professional statistical tables
  "flextable",    # Word/Excel export
  
  # Other utilities
  "haven",        # Read Stata/SPSS files
  "writexl",      # Excel output
  "broom"         # Model results tidying
)

# Install and load packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

cat("All required R packages loaded successfully.\n")

# =============================================================================
# PART 2: Global Path Settings
# =============================================================================

# Project root directory
PROJECT_ROOT <- "C:/Users/user/Desktop/Scientific Research/UKB Project/Jung Sun Lab/Project1_Lifestyle_PPM_Hexiao,Hongtao"

# Data directory
DATA_DIR <- file.path(PROJECT_ROOT, "Data_Excel")

# Code directory
CODE_DIR <- file.path(PROJECT_ROOT, "Code")

# Output directory (create if not exists)
OUTPUT_DIR <- file.path(PROJECT_ROOT, "Output")
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Data file paths for each cohort - UPDATED with new data structure
# Each cohort has: main_analysis, sensitivity_drop1st, mice_data, full_cleaned

# Cohort baseline waves (for variable naming)
BASELINE_WAVES <- list(
  CHARLS = 1,
  ELSA   = 1,   # FIXED: ELSA data uses w1 as baseline (not w7)
  HRS    = 10,  # HRS data uses w10 as baseline
  SHARE  = 4,
  MHAS   = 3
)

# Main analysis data paths
DATA_PATHS <- list(
  CHARLS = file.path(DATA_DIR, "CHARLS", "CHARLS_main_analysis.csv"),
  ELSA   = file.path(DATA_DIR, "ELSA", "ELSA_main_analysis.csv"),
  HRS    = file.path(DATA_DIR, "HRS", "HRS_main_analysis.csv"),
  SHARE  = file.path(DATA_DIR, "SHARE", "SHARE_main_analysis.csv"),
  MHAS   = file.path(DATA_DIR, "MHAS", "MHAS_main_analysis.csv")
)

# Sensitivity analysis data paths (drop first follow-up wave)
SENSITIVITY_PATHS <- list(
  CHARLS = file.path(DATA_DIR, "CHARLS", "CHARLS_sensitivity_drop1st.csv"),
  ELSA   = file.path(DATA_DIR, "ELSA", "ELSA_sensitivity_drop1st.csv"),
  HRS    = file.path(DATA_DIR, "HRS", "HRS_sensitivity_drop1st.csv"),
  SHARE  = file.path(DATA_DIR, "SHARE", "SHARE_sensitivity_drop1st.csv"),
  MHAS   = file.path(DATA_DIR, "MHAS", "MHAS_sensitivity_drop1st.csv")
)

# MICE imputation data paths
MICE_PATHS <- list(
  CHARLS = file.path(DATA_DIR, "CHARLS", "CHARLS_mice_data.csv"),
  ELSA   = file.path(DATA_DIR, "ELSA", "ELSA_mice_data.csv"),
  HRS    = file.path(DATA_DIR, "HRS", "HRS_mice_data.csv"),
  SHARE  = file.path(DATA_DIR, "SHARE", "SHARE_mice_data.csv"),
  MHAS   = file.path(DATA_DIR, "MHAS", "MHAS_mice_data.csv")
)

# Full cleaned data paths
FULL_PATHS <- list(
  CHARLS = file.path(DATA_DIR, "CHARLS", "CHARLS_full_cleaned.csv"),
  ELSA   = file.path(DATA_DIR, "ELSA", "ELSA_full_cleaned.csv"),
  HRS    = file.path(DATA_DIR, "HRS", "HRS_full_cleaned.csv"),
  SHARE  = file.path(DATA_DIR, "SHARE", "SHARE_full_cleaned.csv"),
  MHAS   = file.path(DATA_DIR, "MHAS", "MHAS_full_cleaned.csv")
)

# Cohort country information
COHORT_COUNTRIES <- list(
  CHARLS = "China",
  ELSA   = "England",
  HRS    = "USA",
  SHARE  = "Europe (multi-country)",
  MHAS   = "Mexico"
)

# MICE imputed data paths (after imputation)
MICE_IMPUTED_PATHS <- list(
  CHARLS = file.path(DATA_DIR, "CHARLS", "CHARLS_mice_imputed.csv"),
  ELSA   = file.path(DATA_DIR, "ELSA", "ELSA_mice_imputed.csv"),
  HRS    = file.path(DATA_DIR, "HRS", "HRS_mice_imputed.csv"),
  SHARE  = file.path(DATA_DIR, "SHARE", "SHARE_mice_imputed.csv"),
  MHAS   = file.path(DATA_DIR, "MHAS", "MHAS_mice_imputed.csv"),
  POOLED = file.path(DATA_DIR, "Pooled", "Pooled_mice_imputed.csv")
)

cat("Path settings completed.\n")
cat("Project root:", PROJECT_ROOT, "\n")
cat("Data directory:", DATA_DIR, "\n")

# =============================================================================
# PART 2C: Cox Model Covariate Configuration (CENTRAL DEFINITION)
# =============================================================================
# All Cox regression models should use these covariate settings for consistency

# Primary covariates for pooled analysis
# - age_baseline: continuous age at baseline
# - sex: binary (Men/Women)
# - edu: 3-level education (Primary/Secondary/Tertiary)
# - region: geographic region (7 categories)
# - strata(cohort): stratification by cohort (controls baseline hazard differences)

COX_COVARIATES <- list(
  # For pooled analysis (with cohort stratification)
  pooled = c("age_baseline", "sex", "edu", "region"),
  
  # Stratification variable
  strata_var = "cohort",
  
  # For cohort-specific analysis (meta-analysis)
  cohort_specific = c("age_baseline", "sex", "edu"),
  
  # For SHARE only (has multiple countries within region)
  share_specific = c("age_baseline", "sex", "edu", "region")
)

# Variables excluded from Cox model (due to high missing or collinearity)
COX_EXCLUDED_VARS <- c(
  "rural",         # ELSA 100% missing
  "country_name",  # Too many categories (20), use region instead
  "employment"     # Optional, high missing in some cohorts
)

# Reference categories (for interpretation)
COX_REFERENCE_CATEGORIES <- list(
  sex = "Men",
  edu = "Primary",
  region = "Eastern Asia",
  n_lifestyle_cat = "0"
)

cat("Cox model covariate configuration loaded.\n")
cat("  Pooled covariates:", paste(COX_COVARIATES$pooled, collapse = ", "), "\n")
cat("  Strata variable:", COX_COVARIATES$strata_var, "\n")

# =============================================================================
# PART 2D: MICE Imputation Configuration
# =============================================================================

MICE_CONFIG <- list(
  # Imputation parameters
  m = 20,              # Number of imputed datasets
  maxit = 30,          # Maximum iterations
  seed = 12345,        # Random seed
  
  # Level 2 clustering variable (for 2-level MICE)
  cluster_var = "country_name",
  
  # Variables to impute
  vars_to_impute = c("edu"),  # Primary
  # vars_to_impute = c("edu", "employment"),  # Optional: add employment
  
  # Predictor variables (complete, no missing)
  predictor_vars = c(
    "age_baseline", "sex", "marital_binary", "cohort", "region",
    "physical_base", "psych_base", "cog_base"
  ),
  
  # Excluded from imputation model
  excluded_vars = c("rural")
)

cat("MICE configuration loaded.\n")
cat("  Imputed datasets (m):", MICE_CONFIG$m, "\n")
cat("  Cluster variable:", MICE_CONFIG$cluster_var, "\n")

# =============================================================================
# PART 2E: Variable Name Helper Functions
# =============================================================================

#' Get lifestyle variable names for a specific cohort
#' @param cohort Cohort name (CHARLS, ELSA, HRS, SHARE, MHAS)
#' @return Named list of lifestyle variable names
get_lifestyle_vars <- function(cohort) {
  bw <- BASELINE_WAVES[[cohort]]
  list(
    drink = paste0("w", bw, "_unhealthy_drink"),
    smoke = paste0("w", bw, "_unhealthy_smoke"),
    pa    = paste0("w", bw, "_unhealthy_pa"),
    soc   = paste0("w", bw, "_unhealthy_soc"),
    score = paste0("w", bw, "_unhealthy_score"),
    cat   = paste0("w", bw, "_unhealthy_cat"),
    heavy_drink = paste0("w", bw, "_heavy_drink"),
    score_heavy = paste0("w", bw, "_unhealthy_score_heavy")
  )
}

#' Get outcome event variable names for a specific cohort
#' @param cohort Cohort name
#' @return Named list of outcome variable names
get_outcome_vars <- function(cohort) {
  list(
    # Main PPC-MM outcome
    event_ppcmm = "event_ppcmm",
    time_ppcmm  = "time_ppcmm_months",
    
    # Sensitivity analysis (drop first follow-up)
    event_drop1st = "event_ppcmm_drop1st",
    time_drop1st  = "time_ppcmm_drop1st",
    
    # PPC-MM subtypes
    event_phys_psych = "event_mm_phys_psych",
    time_phys_psych  = "time_mm_phys_psych",
    
    event_phys_cog = "event_mm_phys_cog",
    time_phys_cog  = "time_mm_phys_cog",
    
    event_psych_cog = "event_mm_psych_cog",
    time_psych_cog  = "time_mm_psych_cog",
    
    event_all_three = "event_mm_all_three",
    time_all_three  = "time_mm_all_three"
  )
}

#' Standardize variable names for pooled analysis
#' @param df Data frame
#' @param cohort Cohort name
#' @return Data frame with standardized variable names
standardize_vars <- function(df, cohort) {
  bw <- BASELINE_WAVES[[cohort]]
  
  # Helper function to rename variables using base R (compatible with all dplyr versions)
  rename_if_exists <- function(data, old_name, new_name) {
    if (old_name %in% names(data)) {
      names(data)[names(data) == old_name] <- new_name
    }
    return(data)
  }
  
  # Map cohort-specific lifestyle variables to standard names
  df <- rename_if_exists(df, paste0("w", bw, "_unhealthy_drink"), "unhealthy_drink")
  df <- rename_if_exists(df, paste0("w", bw, "_unhealthy_smoke"), "unhealthy_smoke")
  df <- rename_if_exists(df, paste0("w", bw, "_unhealthy_pa"), "unhealthy_pa")
  df <- rename_if_exists(df, paste0("w", bw, "_unhealthy_soc"), "unhealthy_soc")
  df <- rename_if_exists(df, paste0("w", bw, "_unhealthy_score"), "unhealthy_score")
  df <- rename_if_exists(df, paste0("w", bw, "_unhealthy_cat"), "unhealthy_cat")
  df <- rename_if_exists(df, paste0("w", bw, "_heavy_drink"), "heavy_drink")
  df <- rename_if_exists(df, paste0("w", bw, "_unhealthy_score_heavy"), "unhealthy_score_heavy")
  
  return(df)
}

#' Load cohort data with error handling
#' @param cohort Cohort name
#' @param type Data type: "main", "sensitivity", "mice", "full"
#' @return Data frame or NULL if file not found
load_cohort_data <- function(cohort, type = "main") {
  paths <- switch(type,
    "main" = DATA_PATHS,
    "sensitivity" = SENSITIVITY_PATHS,
    "mice" = MICE_PATHS,
    "full" = FULL_PATHS
  )
  
  file_path <- paths[[cohort]]
  
  if (is.null(file_path) || !file.exists(file_path)) {
    cat("Warning: File not found for", cohort, "-", type, "\n")
    return(NULL)
  }
  
  tryCatch({
    df <- read.csv(file_path, stringsAsFactors = FALSE)
    cat("Loaded", cohort, type, "data: N =", nrow(df), "\n")
    return(df)
  }, error = function(e) {
    cat("Error loading", cohort, ":", e$message, "\n")
    return(NULL)
  })
}

# =============================================================================
# PART 3: Age Grouping Functions
# =============================================================================

#' Convert age to 8-group classification (for cognitive standardization)
#' @param age Age vector
#' @return Age group code (1-8)
age_to_group_8 <- function(age) {
  case_when(
    age >= 45 & age <= 49 ~ 1L,
    age >= 50 & age <= 54 ~ 2L,
    age >= 55 & age <= 59 ~ 3L,
    age >= 60 & age <= 64 ~ 4L,
    age >= 65 & age <= 69 ~ 5L,
    age >= 70 & age <= 74 ~ 6L,
    age >= 75 & age <= 79 ~ 7L,
    age >= 80             ~ 8L,
    TRUE                  ~ NA_integer_
  )
}

#' Convert age to 4-group classification (for stratified analysis)
#' @param age Age vector
#' @return Age group label
age_to_group_4 <- function(age) {
  case_when(
    age >= 50 & age <= 59 ~ "50-59",
    age >= 60 & age <= 69 ~ "60-69",
    age >= 70 & age <= 79 ~ "70-79",
    age >= 80             ~ "80+",
    TRUE                  ~ NA_character_
  )
}

# =============================================================================
# PART 4: Lifestyle Variable Recoding Functions
# =============================================================================

#' Recode drinking - CHARLS version
#' @description Drinking >= once per month is unhealthy (1), otherwise healthy (0)
recode_drink_charls <- function(x) {
  case_when(
    is.na(x)       ~ NA_real_,
    x %in% c(0, 1) ~ 0,      # No drink or < once/month -> healthy
    x >= 2         ~ 1,      # >= once/month -> unhealthy
    TRUE           ~ NA_real_
  )
}

#' Recode drinking - ELSA version
#' @description Any drinking in past 7 days is unhealthy (1)
recode_drink_elsa <- function(x) {
  case_when(
    is.na(x) ~ NA_real_,
    x == 0   ~ 0,            # No drinking -> healthy
    x > 0    ~ 1,            # Any drinking -> unhealthy
    TRUE     ~ NA_real_
  )
}

#' Recode drinking - HRS version
#' @description Drinking days per week > 0 is unhealthy
recode_drink_hrs <- function(days_per_week) {
  case_when(
    is.na(days_per_week) ~ NA_real_,
    days_per_week == 0   ~ 0,
    days_per_week > 0    ~ 1,
    TRUE                 ~ NA_real_
  )
}

#' Recode drinking - SHARE version
#' @description Any drinking in past 3 months is unhealthy
recode_drink_share <- function(x) {
  ifelse(is.na(x), NA_real_, ifelse(x == 1, 1, 0))
}

#' Recode drinking - MHAS version
#' @description Drinking >= 1 day per week is unhealthy
recode_drink_mhas <- function(x) {
  case_when(
    is.na(x) ~ NA_real_,
    x >= 1   ~ 1,
    x == 0   ~ 0,
    TRUE     ~ NA_real_
  )
}

#' Recode smoking - Universal version
#' @description Ever smoker is unhealthy (1), Never smoker is healthy (0)
recode_smoke_ever <- function(ever) {
  case_when(
    is.na(ever) ~ NA_real_,
    ever == 0   ~ 0,         # Never smoker -> healthy
    ever == 1   ~ 1,         # Ever/current smoker -> unhealthy
    TRUE        ~ NA_real_
  )
}

#' Recode physical activity - CHARLS version
#' @description No vigorous and moderate activity is unhealthy
recode_pa_charls <- function(vig, mod) {
  case_when(
    (!is.na(vig) & vig > 0) | (!is.na(mod) & mod > 0) ~ 0,  # Any activity -> healthy
    (!is.na(vig) & vig == 0) & (!is.na(mod) & mod == 0) ~ 1, # No activity -> unhealthy
    TRUE ~ NA_real_
  )
}

#' Recode physical activity - ELSA version
#' @description Frequency 2-3 (at least weekly) is healthy, 4-5 is unhealthy
recode_pa_elsa <- function(vig, mod) {
  case_when(
    (!is.na(vig) & vig %in% c(2, 3)) |
      (!is.na(mod) & mod %in% c(2, 3)) ~ 0,   # At least weekly -> healthy
    (!is.na(vig) & vig %in% c(4, 5)) &
      (!is.na(mod) & mod %in% c(4, 5)) ~ 1,   # Both inactive -> unhealthy
    TRUE ~ NA_real_
  )
}

#' Recode physical activity - HRS version
#' @description 1-3 (at least weekly) is healthy, 4-5 is unhealthy
recode_pa_hrs <- function(vig, mod) {
  case_when(
    (!is.na(vig) & vig %in% 1:3) |
      (!is.na(mod) & mod %in% 1:3) ~ 0,       # Active -> healthy
    (!is.na(vig) & vig %in% 4:5) &
      (!is.na(mod) & mod %in% 4:5) ~ 1,       # Inactive -> unhealthy
    TRUE ~ NA_real_
  )
}

#' Recode physical activity - SHARE version
#' @description 1-3 is healthy, 4 (Hardly ever) is unhealthy
recode_pa_share <- function(v, m) {
  case_when(
    (!is.na(v) & v < 4) | (!is.na(m) & m < 4) ~ 0,
    v == 4 & m == 4 ~ 1,
    TRUE ~ NA_real_
  )
}

#' Recode physical activity - MHAS version
#' @description vigact=1 (>=3 times/week) is healthy, 0 is unhealthy
recode_pa_mhas <- function(x) {
  case_when(
    is.na(x) ~ NA_real_,
    x == 1   ~ 0,            # Active -> healthy
    x == 0   ~ 1,            # Inactive -> unhealthy
    TRUE     ~ NA_real_
  )
}

#' Recode social participation - Universal version
#' @description Has participation (1) is healthy (0), no participation (0) is unhealthy (1)
recode_social <- function(x) {
  case_when(
    is.na(x) ~ NA_real_,
    x == 1   ~ 0,            # Has participation -> healthy
    x == 0   ~ 1,            # No participation -> unhealthy
    TRUE     ~ NA_real_
  )
}

# =============================================================================
# PART 5: Outcome Variable Processing Functions
# =============================================================================

#' Calculate physical disease status
#' @description Any of 7 chronic diseases = 1 is classified as having physical disease
#' @param df Data frame
#' @param disease_vars Disease variable names vector
#' @return 0/1/NA
calculate_physical <- function(df, disease_vars) {
  existing_vars <- disease_vars[disease_vars %in% names(df)]
  if (length(existing_vars) == 0) return(rep(NA_real_, nrow(df)))
  
  tmp <- df[, existing_vars, drop = FALSE]
  n_positive <- rowSums(tmp == 1, na.rm = TRUE)
  n_nonmiss <- rowSums(!is.na(tmp))
  
  ifelse(n_nonmiss == 0, NA_real_,
         ifelse(n_positive > 0, 1, 0))
}

#' Calculate depression status
#' @description Based on CESD score to determine depression (thresholds vary by scale)
#' @param cesd CESD score vector
#' @param threshold Threshold (CESD-10 uses 10, CESD-8 uses 4, MHAS uses 5)
#' @return 0/1/NA
calculate_depression <- function(cesd, threshold) {
  ifelse(is.na(cesd), NA_real_,
         ifelse(cesd >= threshold, 1, 0))
}

#' Calculate cognitive impairment status
#' @description Based on age-group standardized Z-scores, any domain < -1.5 is impairment
#' @param df Data frame with raw cognitive scores and age groups
#' @param mem_var Memory variable name
#' @param orient_var Orientation variable name (optional)
#' @param exec_var Executive function variable name
#' @param agegrp_var Age group variable name
#' @return 0/1/NA
calculate_cognition <- function(df, mem_var, exec_var, agegrp_var, 
                                orient_var = NULL) {
  
  # Prepare data
  vars <- c(agegrp_var, mem_var, exec_var)
  if (!is.null(orient_var) && orient_var %in% names(df)) {
    vars <- c(vars, orient_var)
    has_orient <- TRUE
  } else {
    has_orient <- FALSE
  }
  
  curr_data <- df[, vars[vars %in% names(df)], drop = FALSE]
  
  # Calculate statistics by age group
  stats <- curr_data %>%
    group_by(.data[[agegrp_var]]) %>%
    summarise(across(everything(), 
                     list(m = ~mean(., na.rm = TRUE), 
                          s = ~sd(., na.rm = TRUE))),
              .groups = "drop")
  
  # Merge statistics
  curr_data <- curr_data %>% left_join(stats, by = agegrp_var)
  
  # Calculate Z-scores
  z_mem <- (curr_data[[mem_var]] - curr_data[[paste0(mem_var, "_m")]]) / 
    curr_data[[paste0(mem_var, "_s")]]
  z_exec <- (curr_data[[exec_var]] - curr_data[[paste0(exec_var, "_m")]]) / 
    curr_data[[paste0(exec_var, "_s")]]
  
  if (has_orient) {
    z_orient <- (curr_data[[orient_var]] - curr_data[[paste0(orient_var, "_m")]]) / 
      curr_data[[paste0(orient_var, "_s")]]
  } else {
    z_orient <- NA
  }
  
  # Determine cognitive impairment
  imp_mem <- z_mem < -1.5
  imp_exec <- z_exec < -1.5
  imp_orient <- if (has_orient) z_orient < -1.5 else NA
  
  # Overall determination
  if (has_orient) {
    all_na <- is.na(imp_mem) & is.na(imp_exec) & is.na(imp_orient)
    any_imp <- (imp_mem %in% TRUE) | (imp_exec %in% TRUE) | (imp_orient %in% TRUE)
  } else {
    all_na <- is.na(imp_mem) & is.na(imp_exec)
    any_imp <- (imp_mem %in% TRUE) | (imp_exec %in% TRUE)
  }
  
  ifelse(all_na, NA_real_, ifelse(any_imp, 1, 0))
}

#' Calculate PPC-MM status
#' @description >= 2 of 3 domains impaired is PPC-MM positive
#' @param physical Physical disease status
#' @param psych Psychological/depression status
#' @param cog Cognitive impairment status
#' @return list(n_cond = number of impaired domains, ppc_mm = 0/1/NA)
calculate_ppcmm <- function(physical, psych, cog) {
  n_cond <- rowSums(cbind(physical, psych, cog), na.rm = TRUE)
  all_na <- is.na(physical) & is.na(psych) & is.na(cog)
  n_cond[all_na] <- NA
  
  ppc_mm <- ifelse(is.na(n_cond), NA_real_,
                   ifelse(n_cond >= 2, 1, 0))
  
  list(n_cond = n_cond, ppc_mm = ppc_mm)
}

# =============================================================================
# PART 6: Descriptive Statistics Functions
# =============================================================================

#' Generate baseline characteristics table
#' @param data Analysis dataset
#' @param strata Stratification variable (optional)
#' @param vars Variables to describe
#' @return tableone object
create_table_one <- function(data, vars, strata = NULL, 
                             cat_vars = NULL, nonnormal_vars = NULL) {
  CreateTableOne(
    vars = vars,
    strata = strata,
    data = data,
    factorVars = cat_vars,
    test = TRUE,
    smd = TRUE
  )
}

#' Calculate frequency and percentage for categorical variables
#' @param x Categorical variable
#' @param digits Decimal places
#' @return Formatted string
freq_pct <- function(x, digits = 1) {
  n <- sum(!is.na(x))
  tbl <- table(x, useNA = "no")
  paste0(tbl, " (", round(100 * tbl / n, digits), "%)")
}

#' Calculate mean and standard deviation for continuous variables
#' @param x Continuous variable
#' @param digits Decimal places
#' @return Formatted string
mean_sd <- function(x, digits = 1) {
  paste0(round(mean(x, na.rm = TRUE), digits), " +/- ",
         round(sd(x, na.rm = TRUE), digits))
}

#' Calculate median and interquartile range for continuous variables
#' @param x Continuous variable
#' @param digits Decimal places
#' @return Formatted string
median_iqr <- function(x, digits = 1) {
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  paste0(round(q[2], digits), " [", 
         round(q[1], digits), ", ", 
         round(q[3], digits), "]")
}

# =============================================================================
# PART 7: Cox Regression Standardized Functions
# =============================================================================

#' Extract Cox regression results as table
#' @param fit coxph model object
#' @param model_name Model name
#' @return Data frame
extract_cox_results <- function(fit, model_name = "Model") {
  if (is.null(fit)) return(NULL)
  
  s <- summary(fit)
  
  data.frame(
    Model = model_name,
    Variable = rownames(s$coefficients),
    HR = round(s$coefficients[, "exp(coef)"], 3),
    Lower_CI = round(s$conf.int[, "lower .95"], 3),
    Upper_CI = round(s$conf.int[, "upper .95"], 3),
    P_value = format.pval(s$coefficients[, "Pr(>|z|)"], digits = 3),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# PART 8: Phi Coefficient Calculation Functions
# =============================================================================

#' Calculate Phi coefficient between two binary variables
#' @param x Variable 1
#' @param y Variable 2
#' @return Phi coefficient value
calculate_phi <- function(x, y) {
  # Remove missing values
  complete_idx <- complete.cases(x, y)
  x <- x[complete_idx]
  y <- y[complete_idx]
  
  if (length(x) < 10) return(NA_real_)
  
  # Use cor function (equivalent to Phi for 0/1 variables)
  tryCatch({
    cor(x, y, method = "pearson")
  }, error = function(e) {
    NA_real_
  })
}

#' Calculate Phi coefficient matrix for lifestyle variables
#' @param data Data frame
#' @param lifestyle_vars Lifestyle variable names vector
#' @return Phi coefficient matrix
calculate_phi_matrix <- function(data, lifestyle_vars) {
  n_vars <- length(lifestyle_vars)
  phi_mat <- matrix(NA, nrow = n_vars, ncol = n_vars)
  rownames(phi_mat) <- colnames(phi_mat) <- names(lifestyle_vars)
  
  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      if (i == j) {
        phi_mat[i, j] <- 1
      } else {
        phi_mat[i, j] <- calculate_phi(
          data[[lifestyle_vars[i]]], 
          data[[lifestyle_vars[j]]]
        )
      }
    }
  }
  
  phi_mat
}

# =============================================================================
# PART 9: ICC Calculation Functions
# =============================================================================

#' Calculate Intraclass Correlation Coefficient (ICC)
#' @param outcome Outcome variable name
#' @param cluster Clustering variable name (e.g., cohort)
#' @param data Data frame
#' @return ICC value and confidence interval
calculate_icc <- function(outcome, cluster, data) {
  
  # Fit null model (random intercept model)
  form <- as.formula(paste0(outcome, " ~ 1 + (1|", cluster, ")"))
  
  tryCatch({
    # Use glmer for binary outcomes
    if (all(data[[outcome]] %in% c(0, 1, NA))) {
      fit <- glmer(form, data = data, family = binomial, 
                   control = glmerControl(optimizer = "bobyqa"))
    } else {
      fit <- lmer(form, data = data)
    }
    
    # Extract ICC
    icc_result <- performance::icc(fit)
    
    list(
      ICC_adjusted = icc_result$ICC_adjusted,
      ICC_conditional = icc_result$ICC_conditional
    )
    
  }, error = function(e) {
    cat("ICC calculation error:", e$message, "\n")
    list(ICC_adjusted = NA, ICC_conditional = NA)
  })
}

# =============================================================================
# PART 10: Output and Save Functions
# =============================================================================

#' Save results to Excel
#' @param results_list Results list
#' @param filename File name
#' @param output_dir Output directory
save_to_excel <- function(results_list, filename, output_dir = OUTPUT_DIR) {
  filepath <- file.path(output_dir, filename)
  writexl::write_xlsx(results_list, filepath)
  cat("Results saved to:", filepath, "\n")
}

#' Save results to CSV
#' @param df Data frame
#' @param filename File name
#' @param output_dir Output directory
save_to_csv <- function(df, filename, output_dir = OUTPUT_DIR) {
  filepath <- file.path(output_dir, filename)
  write.csv(df, filepath, row.names = FALSE)
  cat("Results saved to:", filepath, "\n")
}

# =============================================================================
# COMPLETE
# =============================================================================

cat("\n=======================================================\n")
cat("00_Functions_and_Setup.R loaded successfully\n")
cat("All common functions are now available\n")
cat("=======================================================\n")
