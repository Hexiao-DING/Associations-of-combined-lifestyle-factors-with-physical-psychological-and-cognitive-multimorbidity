###############################################################################
# 06_PAF_Analysis.R
# ============================================================================
# Project: Lifestyle Factors and PPC-MM Association Study
# Purpose: Calculate Population Attributable Fraction (PAF) - Lancet Standard
# ============================================================================
#
# ANALYSIS STRUCTURE:
#
# === PRIMARY ANALYSIS (Main Data) ===
#   A. Individual Lifestyle Factors PAF (mutually adjusted HR)
#      - Outcomes: Overall PPC-MM (Primary) + 4 subtypes (Secondary)
#
#   B. Cumulative Lifestyle PAF (4-level: 0/1/2/3+)
#      - Outcomes: Overall PPC-MM (Primary) + 4 subtypes (Secondary)
#
# === SENSITIVITY ANALYSES ===
#   S2: Heavy drinking definition
#   S3: MICE imputed data
#   S4: Drop1st (exclude first follow-up wave)
#
# PAF Formula (Miettinen):
#   PAF = P_case ? (HR - 1) / HR
#   Where P_case = proportion of cases with the exposure
#
###############################################################################

# =============================================================================
# PART 1: Environment Setup
# =============================================================================

source("C:/Users/user/Desktop/Scientific Research/UKB Project/Jung Sun Lab/Project1_Lifestyle_PPM_Hexiao,Hongtao/Code/00_Functions_and_Setup.R")

cat("\n")
cat("================================================================\n")
cat("   Population Attributable Fraction (PAF) Analysis              \n")
cat("   Lancet Standard Framework                                    \n")
cat("================================================================\n\n")

# =============================================================================
# PART 2: Outcome Definitions
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
  # SECONDARY OUTCOMES
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
# PART 3: PAF Calculation Functions
# =============================================================================

#' Calculate PAF for a single exposure level (Miettinen formula)
#' PAF = P_case ? (HR - 1) / HR
calc_paf_single <- function(p_case, hr) {
  paf <- ifelse(is.na(p_case) | is.na(hr) | hr <= 0, 
                NA_real_, 
                p_case * (hr - 1) / hr)
  return(paf)
}

#' Calculate PAF with 95% CI
#' Uses delta method approximation
calc_paf_with_ci <- function(p_case, hr, hr_lower, hr_upper) {
  paf <- calc_paf_single(p_case, hr)
  # CI uses inverse relationship: lower CI of PAF corresponds to upper CI of HR
  paf_lower <- calc_paf_single(p_case, hr_upper)
  paf_upper <- calc_paf_single(p_case, hr_lower)
  return(list(paf = paf, paf_lower = paf_lower, paf_upper = paf_upper))
}

#' Calculate PAF for individual lifestyle factors
calculate_individual_paf <- function(data, outcome_name, covariates, 
                                      use_heavy_drink = FALSE) {
  
  out_info <- outcome_info[[outcome_name]]
  event_var <- out_info$event
  time_var <- out_info$time
  
  if (!event_var %in% names(data) || !time_var %in% names(data)) {
    return(NULL)
  }
  
  # Filter valid data
  df <- data %>%
    filter(!is.na(.data[[event_var]]) & !is.na(.data[[time_var]]) & .data[[time_var]] > 0)
  
  n_events <- sum(df[[event_var]], na.rm = TRUE)
  if (n_events < 20) return(NULL)
  
  # Select drink variable
  drink_var <- ifelse(use_heavy_drink && "heavy_drink" %in% names(df), 
                      "heavy_drink", "unhealthy_drink")
  
  lifestyle_vars <- c(drink_var, "unhealthy_smoke", "unhealthy_pa", "unhealthy_soc")
  lifestyle_labels <- c(
    "unhealthy_drink" = "Drinking (any)",
    "heavy_drink" = "Heavy Drinking",
    "unhealthy_smoke" = "Smoking",
    "unhealthy_pa" = "Physical Inactivity",
    "unhealthy_soc" = "Social Isolation"
  )
  
  # Fit model with all 4 factors mutually adjusted
  formula_str <- paste0("Surv(", time_var, ", ", event_var, ") ~ ",
                        paste(lifestyle_vars, collapse = " + "), " + ",
                        paste(covariates, collapse = " + "))
  
  fit <- tryCatch({
    coxph(as.formula(formula_str), data = df)
  }, error = function(e) NULL)
  
  if (is.null(fit)) return(NULL)
  
  s <- summary(fit)
  cases_only <- df %>% filter(.data[[event_var]] == 1)
  
  results <- list()
  
  for (lv in lifestyle_vars) {
    if (!lv %in% rownames(s$coefficients)) next
    
    hr <- s$coefficients[lv, "exp(coef)"]
    hr_lower <- s$conf.int[lv, "lower .95"]
    hr_upper <- s$conf.int[lv, "upper .95"]
    p_case <- mean(cases_only[[lv]] == 1, na.rm = TRUE)
    
    paf_result <- calc_paf_with_ci(p_case, hr, hr_lower, hr_upper)
    
    results[[lv]] <- data.frame(
      Analysis = ifelse(use_heavy_drink, "S2: Heavy Drink - Individual PAF", "Primary: Individual PAF"),
      Outcome = outcome_name,
      Outcome_Label = out_info$label,
      Outcome_Type = out_info$type,
      Variable = lv,
      Variable_Label = lifestyle_labels[lv],
      HR = round(hr, 3),
      HR_Lower = round(hr_lower, 3),
      HR_Upper = round(hr_upper, 3),
      P_value = format.pval(s$coefficients[lv, "Pr(>|z|)"], digits = 3),
      P_case = round(p_case, 3),
      PAF = round(paf_result$paf, 4),
      PAF_Lower = round(paf_result$paf_lower, 4),
      PAF_Upper = round(paf_result$paf_upper, 4),
      PAF_pct = sprintf("%.1f%%", paf_result$paf * 100),
      N = nrow(df),
      N_Events = n_events,
      stringsAsFactors = FALSE
    )
  }
  
  return(bind_rows(results))
}

#' Calculate PAF for cumulative lifestyle categories
calculate_cumulative_paf <- function(data, outcome_name, exposure_var, 
                                      covariates, analysis_label) {
  
  out_info <- outcome_info[[outcome_name]]
  event_var <- out_info$event
  time_var <- out_info$time
  
  if (!event_var %in% names(data) || !time_var %in% names(data)) {
    return(NULL)
  }
  
  # Filter valid data
  df <- data %>%
    filter(!is.na(.data[[event_var]]) & !is.na(.data[[time_var]]) & .data[[time_var]] > 0)
  
  n_events <- sum(df[[event_var]], na.rm = TRUE)
  if (n_events < 20) return(NULL)
  
  # Fit Cox model
  formula_str <- paste0("Surv(", time_var, ", ", event_var, ") ~ ", 
                        exposure_var, " + ", paste(covariates, collapse = " + "))
  
  fit <- tryCatch({
    coxph(as.formula(formula_str), data = df)
  }, error = function(e) NULL)
  
  if (is.null(fit)) return(NULL)
  
  s <- summary(fit)
  coef_names <- rownames(s$coefficients)
  exp_coefs <- coef_names[grepl(exposure_var, coef_names)]
  
  if (length(exp_coefs) == 0) return(NULL)
  
  # Get case proportions
  cases_only <- df %>% filter(.data[[event_var]] == 1)
  n_cases <- nrow(cases_only)
  
  level_counts <- cases_only %>%
    group_by(.data[[exposure_var]]) %>%
    summarise(n_cases_level = n(), .groups = "drop") %>%
    mutate(p_case = n_cases_level / n_cases)
  
  results <- list()
  
  # Reference level (0)
  ref_p_case <- level_counts$p_case[level_counts[[exposure_var]] == "0"]
  if (length(ref_p_case) == 0) ref_p_case <- 0
  
  results[["ref"]] <- data.frame(
    Analysis = analysis_label,
    Outcome = outcome_name,
    Outcome_Label = out_info$label,
    Outcome_Type = out_info$type,
    Level = "0 (Reference)",
    HR = 1.00,
    HR_Lower = 1.00,
    HR_Upper = 1.00,
    P_value = "-",
    P_case = round(ref_p_case, 3),
    PAF = 0,
    PAF_Lower = 0,
    PAF_Upper = 0,
    PAF_pct = "0.0%",
    N = nrow(df),
    N_Events = n_events,
    stringsAsFactors = FALSE
  )
  
  # Exposure levels
  for (coef_name in exp_coefs) {
    level <- gsub(paste0("^", exposure_var), "", coef_name)
    
    hr <- s$coefficients[coef_name, "exp(coef)"]
    hr_lower <- s$conf.int[coef_name, "lower .95"]
    hr_upper <- s$conf.int[coef_name, "upper .95"]
    p_value <- s$coefficients[coef_name, "Pr(>|z|)"]
    
    # Get case proportion for this level
    level_data <- level_counts %>% filter(as.character(.data[[exposure_var]]) == level)
    p_case <- if (nrow(level_data) > 0) level_data$p_case else 0
    
    paf_result <- calc_paf_with_ci(p_case, hr, hr_lower, hr_upper)
    
    results[[coef_name]] <- data.frame(
      Analysis = analysis_label,
      Outcome = outcome_name,
      Outcome_Label = out_info$label,
      Outcome_Type = out_info$type,
      Level = level,
      HR = round(hr, 3),
      HR_Lower = round(hr_lower, 3),
      HR_Upper = round(hr_upper, 3),
      P_value = format.pval(p_value, digits = 3),
      P_case = round(p_case, 3),
      PAF = round(paf_result$paf, 4),
      PAF_Lower = round(paf_result$paf_lower, 4),
      PAF_Upper = round(paf_result$paf_upper, 4),
      PAF_pct = sprintf("%.1f%%", paf_result$paf * 100),
      N = nrow(df),
      N_Events = n_events,
      stringsAsFactors = FALSE
    )
  }
  
  results_df <- bind_rows(results)
  
  # Calculate combined PAF (sum of level-specific PAFs)
  combined_paf <- sum(results_df$PAF[results_df$Level != "0 (Reference)"], na.rm = TRUE)
  combined_paf <- min(combined_paf, 1.0)  # Cap at 100%
  
  combined_row <- data.frame(
    Analysis = analysis_label,
    Outcome = outcome_name,
    Outcome_Label = out_info$label,
    Outcome_Type = out_info$type,
    Level = "Combined (?? vs 0)",
    HR = NA,
    HR_Lower = NA,
    HR_Upper = NA,
    P_value = "-",
    P_case = round(1 - ref_p_case, 3),
    PAF = round(combined_paf, 4),
    PAF_Lower = NA,
    PAF_Upper = NA,
    PAF_pct = sprintf("%.1f%%", combined_paf * 100),
    N = nrow(df),
    N_Events = n_events,
    stringsAsFactors = FALSE
  )
  
  results_df <- bind_rows(results_df, combined_row)
  
  return(results_df)
}

# =============================================================================
# PART 4: Load and Prepare Data
# =============================================================================

cat("--- Loading and preparing data ---\n")

# Load and pool data (reuse function from 04)
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
    
    df_std <- df %>%
      mutate(cohort = coh, country = COHORT_COUNTRIES[[coh]])
    
    # Rename wave-specific variables
    var_mappings <- c(
      "unhealthy_score" = paste0("w", bw, "_unhealthy_score"),
      "unhealthy_drink" = paste0("w", bw, "_unhealthy_drink"),
      "unhealthy_smoke" = paste0("w", bw, "_unhealthy_smoke"),
      "unhealthy_pa" = paste0("w", bw, "_unhealthy_pa"),
      "unhealthy_soc" = paste0("w", bw, "_unhealthy_soc"),
      "heavy_drink" = paste0("w", bw, "_heavy_drink"),
      "unhealthy_score_heavy" = paste0("w", bw, "_unhealthy_score_heavy")
    )
    
    for (new_name in names(var_mappings)) {
      old_name <- var_mappings[new_name]
      if (old_name %in% names(df_std)) {
        df_std[[new_name]] <- df_std[[old_name]]
      }
    }
    
    # Keep existing country_name and region if they exist
    # (don't rename them, they should already be there)
    
    pooled_list[[coh]] <- df_std
  }
  
  pooled_data <- bind_rows(pooled_list)
  cat("\nPooled data summary: N =", nrow(pooled_data), "\n")
  
  return(pooled_data)
}

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

# Check key variables
cat("\n--- Checking key variables ---\n")
cat("region exists:", "region" %in% names(pooled_main), "\n")
cat("country_name exists:", "country_name" %in% names(pooled_main), "\n")

# Create region variable if not exists
if (!"region" %in% names(pooled_main) || all(is.na(pooled_main$region))) {
  cat("Creating 'region' variable from cohort/country information...\n")
  
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
} else {
  cat("[INFO] Region variable already exists in data.\n")
}

cat("\nRegion distribution:\n")
print(table(pooled_main$region, useNA = "ifany"))

# Prepare data
pooled_data <- pooled_main %>%
  filter(!is.na(age_baseline) & !is.na(sex) & !is.na(edu) &
         !is.na(unhealthy_drink) & !is.na(unhealthy_smoke) &
         !is.na(unhealthy_pa) & !is.na(unhealthy_soc)) %>%
  mutate(
    sex = factor(sex),
    edu = factor(edu),
    cohort = factor(cohort),
    n_lifestyle_cat = factor(
      case_when(
        unhealthy_score == 0 ~ "0",
        unhealthy_score == 1 ~ "1",
        unhealthy_score == 2 ~ "2",
        unhealthy_score >= 3 ~ "3+"
      ), 
      levels = c("0", "1", "2", "3+")
    ),
    n_lifestyle_5cat = factor(unhealthy_score, levels = 0:4)
  )

# Heavy drink version
if ("heavy_drink" %in% names(pooled_data) && "unhealthy_score_heavy" %in% names(pooled_data)) {
  pooled_data <- pooled_data %>%
    mutate(
      n_lifestyle_cat_heavy = factor(
        case_when(
          unhealthy_score_heavy == 0 ~ "0",
          unhealthy_score_heavy == 1 ~ "1",
          unhealthy_score_heavy == 2 ~ "2",
          unhealthy_score_heavy >= 3 ~ "3+"
        ), 
        levels = c("0", "1", "2", "3+")
      )
    )
}

cat("  Prepared data - N =", nrow(pooled_data), "\n")

# Print event counts
cat("\nEvent counts by outcome:\n")
for (out_name in all_outcomes) {
  event_var <- outcome_info[[out_name]]$event
  if (event_var %in% names(pooled_data)) {
    n_events <- sum(pooled_data[[event_var]], na.rm = TRUE)
    cat("  ", out_name, "(", outcome_info[[out_name]]$type, "):", n_events, "events\n")
  }
}

# Use centralized covariate configuration from 00_Functions_and_Setup.R
# Covariates: age_baseline + sex + edu + region
# Note: For PAF, we include cohort as a covariate (not strata) because 
# we need to calculate prevalence across the pooled population

# Check which covariates actually exist
available_covariates <- COX_COVARIATES$pooled[COX_COVARIATES$pooled %in% names(pooled_data)]
if (!"region" %in% names(pooled_data)) {
  available_covariates <- available_covariates[available_covariates != "region"]
  cat("[WARNING] 'region' variable not found, using reduced covariate set.\n")
}

# Ensure region is a factor if it exists
if ("region" %in% names(pooled_data)) {
  pooled_data$region <- factor(pooled_data$region)
}

covariates <- c(available_covariates, "cohort")

cat("\nPAF Analysis Covariates:", paste(covariates, collapse = " + "), "\n")

# =============================================================================
# Initialize all result data frames (to prevent errors if some analyses fail)
# =============================================================================

cat("\n[INFO] Initializing PAF result containers...\n")
# Primary analyses
individual_paf_df <- data.frame()
cumulative_paf_df <- data.frame()
# Sensitivity analyses
sens2a_df <- data.frame()
sens2b_df <- data.frame()
sens3_individual_df <- data.frame()
sens3_cumulative_df <- data.frame()
sens4_df <- data.frame()
# Subgroup analyses (IMPORTANT: must initialize)
age_stratified_individual_df <- data.frame()
age_stratified_cumulative_df <- data.frame()
sex_stratified_individual_df <- data.frame()
sex_stratified_cumulative_df <- data.frame()
all_paf_results <- data.frame()
cat("[INFO] All PAF result containers initialized.\n")

# =============================================================================
# PART 5: PRIMARY ANALYSIS - Individual Factors PAF
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   PRIMARY ANALYSIS A: Individual Lifestyle Factors PAF          \n")
cat("================================================================\n")

individual_paf_results <- list()

for (out_name in all_outcomes) {
  cat("\n--- Outcome:", outcome_info[[out_name]]$label, 
      "(", outcome_info[[out_name]]$type, ") ---\n")
  
  result <- calculate_individual_paf(pooled_data, out_name, covariates, use_heavy_drink = FALSE)
  
  if (!is.null(result)) {
    individual_paf_results[[out_name]] <- result
    print(result %>% select(Variable_Label, HR, P_case, PAF_pct))
  }
}

individual_paf_df <- bind_rows(individual_paf_results)

# IMMEDIATE SAVE: Save primary individual PAF results right away
cat("\n--- Immediate Save: Primary Individual PAF Results ---\n")
tryCatch({
  if (nrow(individual_paf_df) > 0) {
    write.csv(individual_paf_df, file.path(OUTPUT_DIR, "PAF_Primary_Individual.csv"), row.names = FALSE)
    cat("  [OK] Saved PAF_Primary_Individual.csv (", nrow(individual_paf_df), " rows)\n")
  } else {
    cat("  [WARN] No individual PAF results to save\n")
  }
}, error = function(e) {
  cat("  [ERROR] Failed to save individual PAF results:", e$message, "\n")
})

# =============================================================================
# PART 6: PRIMARY ANALYSIS - Cumulative Effect PAF (4-level)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   PRIMARY ANALYSIS B: Cumulative Lifestyle PAF (4-level)        \n")
cat("================================================================\n")

cumulative_paf_results <- list()

for (out_name in all_outcomes) {
  cat("\n--- Outcome:", outcome_info[[out_name]]$label, 
      "(", outcome_info[[out_name]]$type, ") ---\n")
  
  result <- calculate_cumulative_paf(
    pooled_data, out_name, "n_lifestyle_cat", covariates,
    "Primary: Cumulative PAF (4-level)"
  )
  
  if (!is.null(result)) {
    cumulative_paf_results[[out_name]] <- result
    combined <- result %>% filter(Level == "Combined (?? vs 0)")
    cat("  Combined PAF:", combined$PAF_pct, "\n")
  }
}

cumulative_paf_df <- bind_rows(cumulative_paf_results)

# IMMEDIATE SAVE: Save primary cumulative PAF results right away
cat("\n--- Immediate Save: Primary Cumulative PAF Results ---\n")
tryCatch({
  if (nrow(cumulative_paf_df) > 0) {
    write.csv(cumulative_paf_df, file.path(OUTPUT_DIR, "PAF_Primary_Cumulative.csv"), row.names = FALSE)
    cat("  [OK] Saved PAF_Primary_Cumulative.csv (", nrow(cumulative_paf_df), " rows)\n")
  } else {
    cat("  [WARN] No cumulative PAF results to save\n")
  }
}, error = function(e) {
  cat("  [ERROR] Failed to save cumulative PAF results:", e$message, "\n")
})

# =============================================================================
# PART 7: SENSITIVITY S2 - Heavy Drinking Definition
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   SENSITIVITY S2: Heavy Drinking Definition PAF                 \n")
cat("================================================================\n")

# S2a: Individual factors with heavy drink
sens2a_results <- list()

for (out_name in all_outcomes) {
  cat("\n--- Individual PAF (Heavy Drink) - Outcome:", 
      outcome_info[[out_name]]$label, "---\n")
  
  result <- calculate_individual_paf(pooled_data, out_name, covariates, use_heavy_drink = TRUE)
  
  if (!is.null(result)) {
    sens2a_results[[out_name]] <- result
  }
}

sens2a_df <- bind_rows(sens2a_results)

# S2b: Cumulative with heavy drink score
sens2b_results <- list()

if ("n_lifestyle_cat_heavy" %in% names(pooled_data)) {
  for (out_name in all_outcomes) {
    cat("\n--- Cumulative PAF (Heavy Drink) - Outcome:", 
        outcome_info[[out_name]]$label, "---\n")
    
    result <- calculate_cumulative_paf(
      pooled_data, out_name, "n_lifestyle_cat_heavy", covariates,
      "S2: Heavy Drink - Cumulative PAF"
    )
    
    if (!is.null(result)) {
      sens2b_results[[out_name]] <- result
      combined <- result %>% filter(Level == "Combined (?? vs 0)")
      cat("  Combined PAF:", combined$PAF_pct, "\n")
    }
  }
}

sens2b_df <- bind_rows(sens2b_results)

# =============================================================================
# PART 8: SENSITIVITY S3 - MICE Data
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   SENSITIVITY S3: MICE Imputed Data PAF                         \n")
cat("================================================================\n")

sens3_individual <- list()
sens3_cumulative <- list()

# Check for MICE imputed data
mice_csv_path <- file.path(DATA_DIR, "Pooled", "Pooled_mice_imputed.csv")

cat("\nChecking MICE data availability:\n")
cat("  CSV file:", ifelse(file.exists(mice_csv_path), "FOUND", "Not found"), "\n")

if (file.exists(mice_csv_path)) {
  pooled_mice <- read.csv(mice_csv_path, stringsAsFactors = FALSE)
  cat("  Loaded:", nrow(pooled_mice), "rows\n")
  
  if ("event_ppcmm" %in% names(pooled_mice)) {
    pooled_mice <- pooled_mice %>%
      filter(!is.na(age_baseline) & !is.na(sex) & !is.na(edu) &
             !is.na(unhealthy_drink) & !is.na(unhealthy_smoke) &
             !is.na(unhealthy_pa) & !is.na(unhealthy_soc)) %>%
      mutate(
        sex = factor(sex),
        edu = factor(edu),
        cohort = factor(cohort),
        n_lifestyle_cat = factor(
          case_when(
            unhealthy_score == 0 ~ "0",
            unhealthy_score == 1 ~ "1",
            unhealthy_score == 2 ~ "2",
            unhealthy_score >= 3 ~ "3+"
          ), 
          levels = c("0", "1", "2", "3+")
        )
      )
    
    cat("  MICE analysis sample size:", nrow(pooled_mice), "\n")
    
    if (nrow(pooled_mice) > 100) {
      # Individual factors
      cat("\n--- S3: Individual Factors PAF (MICE) ---\n")
      for (out_name in all_outcomes) {
        result <- calculate_individual_paf(pooled_mice, out_name, covariates, use_heavy_drink = FALSE)
        if (!is.null(result)) {
          result$Analysis <- "S3: MICE - Individual PAF"
          sens3_individual[[out_name]] <- result
        }
      }
      
      # Cumulative effect
      cat("\n--- S3: Cumulative Effect PAF (MICE) ---\n")
      for (out_name in all_outcomes) {
        result <- calculate_cumulative_paf(
          pooled_mice, out_name, "n_lifestyle_cat", covariates,
          "S3: MICE - Cumulative PAF"
        )
        if (!is.null(result)) {
          sens3_cumulative[[out_name]] <- result
        }
      }
    }
  } else {
    cat("  [WARN] Outcome variables not found in MICE data\n")
  }
} else {
  cat("\n  [INFO] No MICE data available.\n")
  cat("  Run 02_MICE_Imputation.R first to generate MICE imputed data.\n")
  cat("  Skipping S3 sensitivity analysis.\n")
}

sens3_individual_df <- bind_rows(sens3_individual)
sens3_cumulative_df <- bind_rows(sens3_cumulative)

cat("\n  S3 PAF Results:\n")
cat("    Individual:", nrow(sens3_individual_df), "rows\n")
cat("    Cumulative:", nrow(sens3_cumulative_df), "rows\n")

# =============================================================================
# PART 9: SENSITIVITY S4 - Drop First Follow-up
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   SENSITIVITY S4: Drop First Follow-up PAF                      \n")
cat("================================================================\n")

sens4_results <- list()

if ("event_ppcmm_drop1st" %in% names(pooled_data)) {
  
  df_drop1st <- pooled_data %>%
    filter(!is.na(event_ppcmm_drop1st) & !is.na(time_ppcmm_drop1st) & time_ppcmm_drop1st > 0)
  
  n_events <- sum(df_drop1st$event_ppcmm_drop1st, na.rm = TRUE)
  cat("Drop1st events:", n_events, "\n")
  
  if (n_events >= 20) {
    # Individual factors
    lifestyle_vars <- c("unhealthy_drink", "unhealthy_smoke", "unhealthy_pa", "unhealthy_soc")
    lifestyle_labels <- c("Drinking", "Smoking", "Physical Inactivity", "Social Isolation")
    
    formula_str <- paste0("Surv(time_ppcmm_drop1st, event_ppcmm_drop1st) ~ ",
                          paste(lifestyle_vars, collapse = " + "), " + ",
                          paste(covariates, collapse = " + "))
    
    fit <- tryCatch({
      coxph(as.formula(formula_str), data = df_drop1st)
    }, error = function(e) NULL)
    
    if (!is.null(fit)) {
      s <- summary(fit)
      cases_only <- df_drop1st %>% filter(event_ppcmm_drop1st == 1)
      
      indiv_results <- list()
      
      for (i in seq_along(lifestyle_vars)) {
        lv <- lifestyle_vars[i]
        if (!lv %in% rownames(s$coefficients)) next
        
        hr <- s$coefficients[lv, "exp(coef)"]
        hr_lower <- s$conf.int[lv, "lower .95"]
        hr_upper <- s$conf.int[lv, "upper .95"]
        p_case <- mean(cases_only[[lv]] == 1, na.rm = TRUE)
        
        paf_result <- calc_paf_with_ci(p_case, hr, hr_lower, hr_upper)
        
        indiv_results[[lv]] <- data.frame(
          Analysis = "S4: Drop1st - Individual PAF",
          Outcome = "Overall_drop1st",
          Outcome_Label = "Overall PPC-MM (Drop1st)",
          Variable = lv,
          Variable_Label = lifestyle_labels[i],
          HR = round(hr, 3),
          P_case = round(p_case, 3),
          PAF = round(paf_result$paf, 4),
          PAF_pct = sprintf("%.1f%%", paf_result$paf * 100),
          N = nrow(df_drop1st),
          N_Events = n_events,
          stringsAsFactors = FALSE
        )
      }
      
      sens4_results[["individual"]] <- bind_rows(indiv_results)
      cat("\n--- Individual Factors PAF (Drop1st) ---\n")
      print(sens4_results[["individual"]] %>% select(Variable_Label, HR, PAF_pct))
    }
    
    # Cumulative effect
    formula_str2 <- paste0("Surv(time_ppcmm_drop1st, event_ppcmm_drop1st) ~ ",
                           "n_lifestyle_cat + ", paste(covariates, collapse = " + "))
    
    fit2 <- tryCatch({
      coxph(as.formula(formula_str2), data = df_drop1st)
    }, error = function(e) NULL)
    
    if (!is.null(fit2)) {
      s2 <- summary(fit2)
      exp_coefs <- rownames(s2$coefficients)[grepl("n_lifestyle_cat", rownames(s2$coefficients))]
      cases_only <- df_drop1st %>% filter(event_ppcmm_drop1st == 1)
      
      level_counts <- cases_only %>%
        group_by(n_lifestyle_cat) %>%
        summarise(n_cases_level = n(), .groups = "drop") %>%
        mutate(p_case = n_cases_level / nrow(cases_only))
      
      cum_results <- list()
      
      for (coef_name in exp_coefs) {
        level <- gsub("n_lifestyle_cat", "", coef_name)
        hr <- s2$coefficients[coef_name, "exp(coef)"]
        
        level_data <- level_counts %>% filter(as.character(n_lifestyle_cat) == level)
        p_case <- if (nrow(level_data) > 0) level_data$p_case else 0
        
        paf <- calc_paf_single(p_case, hr)
        
        cum_results[[coef_name]] <- data.frame(
          Analysis = "S4: Drop1st - Cumulative PAF",
          Outcome = "Overall_drop1st",
          Level = level,
          HR = round(hr, 3),
          P_case = round(p_case, 3),
          PAF = round(paf, 4),
          PAF_pct = sprintf("%.1f%%", paf * 100),
          stringsAsFactors = FALSE
        )
      }
      
      sens4_results[["cumulative"]] <- bind_rows(cum_results)
      cat("\n--- Cumulative PAF (Drop1st) ---\n")
      print(sens4_results[["cumulative"]])
    }
  }
} else {
  cat("Drop1st variables not available, skipping S4\n")
}

sens4_df <- bind_rows(sens4_results)

# =============================================================================
# PART 10: SUBGROUP ANALYSIS - Age Stratification
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   SUBGROUP ANALYSIS: Age-Stratified PAF                         \n")
cat("================================================================\n")

# Verify data availability
cat("\n>>> Checking pooled_data...\n")
if (!exists("pooled_data") || is.null(pooled_data) || nrow(pooled_data) == 0) {
  cat("  [ERROR] pooled_data not available. Trying to recreate...\n")
  if (exists("pooled_main") && nrow(pooled_main) > 0) {
    pooled_data <- pooled_main %>%
      filter(!is.na(age_baseline) & !is.na(sex) & !is.na(edu) &
             !is.na(unhealthy_drink) & !is.na(unhealthy_smoke) &
             !is.na(unhealthy_pa) & !is.na(unhealthy_soc)) %>%
      mutate(
        sex = factor(sex),
        edu = factor(edu),
        cohort = factor(cohort),
        n_lifestyle_cat = factor(
          case_when(
            unhealthy_score == 0 ~ "0",
            unhealthy_score == 1 ~ "1",
            unhealthy_score == 2 ~ "2",
            unhealthy_score >= 3 ~ "3+"
          ), levels = c("0", "1", "2", "3+")
        )
      )
    cat("  [OK] Recreated pooled_data with", nrow(pooled_data), "rows\n")
  }
}

# Create age groups
cat(">>> Creating age groups...\n")
pooled_data <- pooled_data %>%
  mutate(
    age_group_4 = case_when(
      age_baseline >= 50 & age_baseline < 60 ~ "50-59",
      age_baseline >= 60 & age_baseline < 70 ~ "60-69",
      age_baseline >= 70 & age_baseline < 80 ~ "70-79",
      age_baseline >= 80 ~ "80+",
      TRUE ~ NA_character_
    )
  )

cat("\nAge group distribution with events:\n")
age_event_table <- pooled_data %>%
  group_by(age_group_4) %>%
  summarise(N = n(), Events = sum(event_ppcmm, na.rm=TRUE), .groups="drop")
print(age_event_table)

# Local variable definitions
lf_vars <- c("unhealthy_drink", "unhealthy_smoke", "unhealthy_pa", "unhealthy_soc")
lf_labels <- c("Drinking", "Smoking", "Physical Inactivity", "Social Isolation")

age_groups <- c("50-59", "60-69", "70-79", "80+")
age_stratified_paf <- list()

cat("\n--- Individual Factors PAF by Age Group ---\n")

for (ag in age_groups) {
  cat("\n  Age Group:", ag, "\n")
  
  df_age <- pooled_data %>% filter(age_group_4 == ag)
  n_total <- nrow(df_age)
  n_events <- sum(df_age$event_ppcmm, na.rm = TRUE)
  cat("    N =", n_total, ", Events =", n_events, "\n")
  
  if (n_events < 10) {
    cat("    [SKIP] Insufficient events\n")
    next
  }
  
  # Build formula - handle multiple cohorts
  covs <- c("sex", "edu", "region")
  covs <- covs[covs %in% names(df_age)]
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
    cat("    [ERROR] Model failed:", e$message, "\n")
    NULL
  })
  
  if (!is.null(fit)) {
    s <- summary(fit)
    cases_only <- df_age %>% filter(event_ppcmm == 1)
    
    for (i in seq_along(lf_vars)) {
      lv <- lf_vars[i]
      if (!lv %in% rownames(s$coefficients)) next
      
      hr <- s$coefficients[lv, "exp(coef)"]
      hr_lower <- s$conf.int[lv, "lower .95"]
      hr_upper <- s$conf.int[lv, "upper .95"]
      p_value <- s$coefficients[lv, "Pr(>|z|)"]
      p_case <- mean(cases_only[[lv]] == 1, na.rm = TRUE)
      
      paf_result <- calc_paf_with_ci(p_case, hr, hr_lower, hr_upper)
      
      age_stratified_paf[[paste0(ag, "_", lv)]] <- data.frame(
        Analysis = "Age-Stratified PAF",
        Age_Group = ag,
        Outcome = "Overall",
        Outcome_Label = "Overall PPC-MM",
        Outcome_Type = "Primary",
        Variable = lv,
        Variable_Label = lf_labels[i],
        HR = round(hr, 3),
        HR_Lower = round(hr_lower, 3),
        HR_Upper = round(hr_upper, 3),
        HR_CI = sprintf("%.2f (%.2f-%.2f)", hr, hr_lower, hr_upper),
        P_value = p_value,
        P_value_fmt = format.pval(p_value, digits = 4),
        Significance = case_when(
          p_value < 0.001 ~ "***",
          p_value < 0.01 ~ "**",
          p_value < 0.05 ~ "*",
          TRUE ~ ""
        ),
        P_case = round(p_case, 4),
        P_case_pct = sprintf("%.1f%%", p_case * 100),
        PAF = round(paf_result$paf, 4),
        PAF_Lower = round(paf_result$paf_lower, 4),
        PAF_Upper = round(paf_result$paf_upper, 4),
        PAF_pct = sprintf("%.1f%%", paf_result$paf * 100),
        PAF_CI = sprintf("%.1f%% (%.1f%%, %.1f%%)", 
                         paf_result$paf * 100, 
                         paf_result$paf_lower * 100, 
                         paf_result$paf_upper * 100),
        N = n_total,
        N_Events = n_events,
        stringsAsFactors = FALSE
      )
    }
    cat("    [OK] Individual factors PAF calculated\n")
  }
}

# Combine age-stratified individual PAF
if (length(age_stratified_paf) > 0) {
  age_stratified_individual_df <- bind_rows(age_stratified_paf)
} else {
  age_stratified_individual_df <- data.frame()
}

# Save ALWAYS
cat("\n=== Age-Stratified PAF Summary (Individual Factors) ===\n")
cat("Total rows:", nrow(age_stratified_individual_df), "\n")

tryCatch({
  paf_age_ind_file <- file.path(OUTPUT_DIR, "PAF_Age_Stratified_Individual.csv")
  if (nrow(age_stratified_individual_df) > 0) {
    summary_table <- age_stratified_individual_df %>%
      select(Age_Group, Variable_Label, HR_CI, PAF_CI, N_Events) %>%
      arrange(Age_Group, Variable_Label)
    print(summary_table, n = 50)
    write.csv(age_stratified_individual_df, paf_age_ind_file, row.names = FALSE)
  } else {
    placeholder <- data.frame(
      Note = "No age-stratified individual PAF results",
      Age_Groups = "50-59, 60-69, 70-79, 80+"
    )
    write.csv(placeholder, paf_age_ind_file, row.names = FALSE)
  }
  cat("\n  [OK] Saved: PAF_Age_Stratified_Individual.csv\n")
}, error = function(e) {
  cat("  [ERROR] Failed to save:", e$message, "\n")
})

# --- Cumulative PAF by Age Group ---
cat("\n--- Cumulative PAF by Age Group ---\n")

age_stratified_cumulative <- list()

for (ag in age_groups) {
  cat("\n  Age Group:", ag, "\n")
  
  df_age <- pooled_data %>%
    filter(age_group_4 == ag &
           !is.na(event_ppcmm) & !is.na(time_ppcmm_months) & time_ppcmm_months > 0)
  
  n_events <- sum(df_age$event_ppcmm, na.rm = TRUE)
  
  if (n_events < 20) {
    cat("    [SKIP] Insufficient events (<20)\n")
    next
  }
  
  # Cumulative effect
  age_cov <- covariates[covariates != "age_baseline"]
  formula_str <- paste0("Surv(time_ppcmm_months, event_ppcmm) ~ ",
                        "n_lifestyle_cat + ", paste(age_cov, collapse = " + "))
  
  fit <- tryCatch({
    coxph(as.formula(formula_str), data = df_age)
  }, error = function(e) NULL)
  
  if (!is.null(fit)) {
    s <- summary(fit)
    exp_coefs <- rownames(s$coefficients)[grepl("n_lifestyle_cat", rownames(s$coefficients))]
    cases_only <- df_age %>% filter(event_ppcmm == 1)
    
    level_counts <- cases_only %>%
      group_by(n_lifestyle_cat) %>%
      summarise(n_cases_level = n(), .groups = "drop") %>%
      mutate(p_case = n_cases_level / nrow(cases_only))
    
    combined_paf <- 0
    
    for (coef_name in exp_coefs) {
      level <- gsub("n_lifestyle_cat", "", coef_name)
      hr <- s$coefficients[coef_name, "exp(coef)"]
      hr_lower <- s$conf.int[coef_name, "lower .95"]
      hr_upper <- s$conf.int[coef_name, "upper .95"]
      p_value <- s$coefficients[coef_name, "Pr(>|z|)"]
      
      level_data <- level_counts %>% filter(as.character(n_lifestyle_cat) == level)
      p_case <- if (nrow(level_data) > 0) level_data$p_case else 0
      
      paf_result <- calc_paf_with_ci(p_case, hr, hr_lower, hr_upper)
      combined_paf <- combined_paf + paf_result$paf
      
      age_stratified_cumulative[[paste0(ag, "_", level)]] <- data.frame(
        Analysis = "Subgroup: Age-Stratified Cumulative PAF",
        Age_Group = ag,
        Outcome = "Overall",
        Outcome_Label = "Overall PPC-MM",
        Level = level,
        HR = round(hr, 3),
        HR_Lower = round(hr_lower, 3),
        HR_Upper = round(hr_upper, 3),
        HR_CI = sprintf("%.2f (%.2f-%.2f)", hr, hr_lower, hr_upper),
        P_value = p_value,
        P_value_fmt = format.pval(p_value, digits = 4),
        Significance = case_when(
          p_value < 0.001 ~ "***",
          p_value < 0.01 ~ "**",
          p_value < 0.05 ~ "*",
          TRUE ~ ""
        ),
        P_case = round(p_case, 4),
        P_case_pct = sprintf("%.1f%%", p_case * 100),
        PAF = round(paf_result$paf, 4),
        PAF_pct = sprintf("%.1f%%", paf_result$paf * 100),
        N = nrow(df_age),
        N_Events = n_events,
        stringsAsFactors = FALSE
      )
    }
    
    # Add combined PAF row
    age_stratified_cumulative[[paste0(ag, "_combined")]] <- data.frame(
      Analysis = "Subgroup: Age-Stratified Cumulative PAF",
      Age_Group = ag,
      Outcome = "Overall",
      Outcome_Label = "Overall PPC-MM",
      Level = "Combined (≥1 vs 0)",
      HR = NA,
      HR_Lower = NA,
      HR_Upper = NA,
      HR_CI = "-",
      P_value = NA,
      P_value_fmt = "-",
      Significance = "",
      P_case = NA,
      P_case_pct = "-",
      PAF = round(combined_paf, 4),
      PAF_pct = sprintf("%.1f%%", combined_paf * 100),
      N = nrow(df_age),
      N_Events = n_events,
      stringsAsFactors = FALSE
    )
    
    cat("    [OK] Combined PAF:", sprintf("%.1f%%", combined_paf * 100), "\n")
  }
}

# Combine age-stratified cumulative PAF
if (length(age_stratified_cumulative) > 0) {
  age_stratified_cumulative_df <- bind_rows(age_stratified_cumulative)
} else {
  age_stratified_cumulative_df <- data.frame()
}

# Save ALWAYS
cat("\n=== Age-Stratified Cumulative PAF Summary ===\n")
cat("Total rows:", nrow(age_stratified_cumulative_df), "\n")

tryCatch({
  paf_age_cum_file <- file.path(OUTPUT_DIR, "PAF_Age_Stratified_Cumulative.csv")
  if (nrow(age_stratified_cumulative_df) > 0) {
    summary_cum <- age_stratified_cumulative_df %>%
      select(Age_Group, Level, HR_CI, PAF_pct, N_Events) %>%
      arrange(Age_Group, Level)
    print(summary_cum, n = 50)
    write.csv(age_stratified_cumulative_df, paf_age_cum_file, row.names = FALSE)
  } else {
    placeholder <- data.frame(
      Note = "No age-stratified cumulative PAF results",
      Age_Groups = "50-59, 60-69, 70-79, 80+"
    )
    write.csv(placeholder, paf_age_cum_file, row.names = FALSE)
  }
  cat("\n  [OK] Saved: PAF_Age_Stratified_Cumulative.csv\n")
}, error = function(e) {
  cat("  [ERROR] Failed to save:", e$message, "\n")
})

# Combine age-stratified results
age_stratified_all_df <- bind_rows(
  age_stratified_individual_df,
  age_stratified_cumulative_df
)

# =============================================================================
# PART 10B: SUBGROUP ANALYSIS - Sex-Stratified PAF
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   SUBGROUP ANALYSIS: Sex-Stratified PAF                         \n")
cat("================================================================\n")

cat("\nSex distribution with events:\n")
sex_event_table <- pooled_data %>%
  group_by(sex) %>%
  summarise(N = n(), Events = sum(event_ppcmm, na.rm=TRUE), .groups="drop")
print(sex_event_table)

sex_groups <- c("Men", "Women")
sex_stratified_paf <- list()

cat("\n--- Individual Factors PAF by Sex ---\n")

for (sg in sex_groups) {
  cat("\n  Sex:", sg, "\n")
  
  df_sex <- pooled_data %>% filter(sex == sg)
  n_total <- nrow(df_sex)
  n_events <- sum(df_sex$event_ppcmm, na.rm = TRUE)
  cat("    N =", n_total, ", Events =", n_events, "\n")
  
  if (n_events < 10) {
    cat("    [SKIP] Insufficient events\n")
    next
  }
  
  # Build formula with local variables
  covs <- c("age_baseline", "edu", "region")
  covs <- covs[covs %in% names(df_sex)]
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
    cat("    [ERROR] Model failed:", e$message, "\n")
    NULL
  })
  
  if (!is.null(fit)) {
    s <- summary(fit)
    cases_only <- df_sex %>% filter(event_ppcmm == 1)
    
    for (i in seq_along(lf_vars)) {
      lv <- lf_vars[i]
      if (!lv %in% rownames(s$coefficients)) next
      
      hr <- s$coefficients[lv, "exp(coef)"]
      hr_lower <- s$conf.int[lv, "lower .95"]
      hr_upper <- s$conf.int[lv, "upper .95"]
      p_val <- s$coefficients[lv, "Pr(>|z|)"]
      p_case <- mean(cases_only[[lv]] == 1, na.rm = TRUE)
      paf_result <- calc_paf_with_ci(p_case, hr, hr_lower, hr_upper)
      
      sex_stratified_paf[[paste0(sg, "_", lv)]] <- data.frame(
        Analysis = "Sex-Stratified PAF",
        Sex = sg,
        Outcome = "Overall",
        Outcome_Label = "Overall PPC-MM",
        Outcome_Type = "Primary",
        Variable = lv,
        Variable_Label = lf_labels[i],
        HR = round(hr, 3),
        HR_Lower = round(hr_lower, 3),
        HR_Upper = round(hr_upper, 3),
        HR_CI = sprintf("%.2f (%.2f-%.2f)", hr, hr_lower, hr_upper),
        P_value = p_val,
        P_value_fmt = format.pval(p_val, digits = 4),
        Significance = case_when(
          p_val < 0.001 ~ "***",
          p_val < 0.01 ~ "**",
          p_val < 0.05 ~ "*",
          TRUE ~ ""
        ),
        P_case = round(p_case, 4),
        P_case_pct = sprintf("%.1f%%", p_case * 100),
        PAF = round(paf_result$paf, 4),
        PAF_Lower = round(paf_result$paf_lower, 4),
        PAF_Upper = round(paf_result$paf_upper, 4),
        PAF_pct = sprintf("%.1f%%", paf_result$paf * 100),
        PAF_CI = sprintf("%.1f%% (%.1f%%, %.1f%%)", 
                         paf_result$paf * 100, 
                         paf_result$paf_lower * 100, 
                         paf_result$paf_upper * 100),
        N = n_total,
        N_Events = n_events,
        stringsAsFactors = FALSE
      )
    }
    cat("    [OK] Individual PAF calculated\n")
  }
}

# Combine sex-stratified individual PAF
if (length(sex_stratified_paf) > 0) {
  sex_stratified_individual_df <- bind_rows(sex_stratified_paf)
} else {
  sex_stratified_individual_df <- data.frame()
}

# Save ALWAYS
cat("\n=== Sex-Stratified PAF Summary (Individual Factors) ===\n")
cat("Total rows:", nrow(sex_stratified_individual_df), "\n")

tryCatch({
  paf_sex_ind_file <- file.path(OUTPUT_DIR, "PAF_Sex_Stratified_Individual.csv")
  if (nrow(sex_stratified_individual_df) > 0) {
    summary_table <- sex_stratified_individual_df %>%
      select(Sex, Variable_Label, HR_CI, PAF_CI, N_Events) %>%
      arrange(Sex, Variable_Label)
    print(summary_table, n = 50)
    write.csv(sex_stratified_individual_df, paf_sex_ind_file, row.names = FALSE)
  } else {
    placeholder <- data.frame(
      Note = "No sex-stratified individual PAF results",
      Sex_Groups = "Men, Women"
    )
    write.csv(placeholder, paf_sex_ind_file, row.names = FALSE)
  }
  cat("\n  [OK] Saved: PAF_Sex_Stratified_Individual.csv\n")
}, error = function(e) {
  cat("  [ERROR] Failed to save:", e$message, "\n")
})

# --- Cumulative Effect by Sex ---
cat("\n--- Cumulative Effect PAF by Sex ---\n")

sex_stratified_cumulative <- list()

for (sg in sex_groups) {
  cat("\n  Sex:", sg, "\n")
  
  df_sex <- pooled_data %>%
    filter(sex == sg &
           !is.na(event_ppcmm) & !is.na(time_ppcmm_months) & time_ppcmm_months > 0)
  
  n_total <- nrow(df_sex)
  n_events <- sum(df_sex$event_ppcmm, na.rm = TRUE)
  
  if (n_events < 20) {
    cat("    [SKIP] Insufficient events (<20)\n")
    next
  }
  
  sex_cov <- covariates[covariates != "sex"]
  
  formula_cum <- paste0("Surv(time_ppcmm_months, event_ppcmm) ~ ",
                        "n_lifestyle_cat + ", paste(sex_cov, collapse = " + "), " + strata(cohort)")
  
  fit_cum <- tryCatch({
    coxph(as.formula(formula_cum), data = df_sex)
  }, error = function(e) NULL)
  
  if (!is.null(fit_cum)) {
    s <- summary(fit_cum)
    exp_rows <- grepl("n_lifestyle_cat", rownames(s$coefficients))
    cases_only <- df_sex %>% filter(event_ppcmm == 1)
    
    # Get level distribution in cases
    level_counts <- cases_only %>%
      group_by(n_lifestyle_cat) %>%
      summarise(n_cases_level = n(), .groups = "drop") %>%
      mutate(p_case = n_cases_level / nrow(cases_only))
    
    combined_paf <- 0
    
    for (coef_name in rownames(s$coefficients)[exp_rows]) {
      level <- gsub("n_lifestyle_cat", "", coef_name)
      hr <- s$coefficients[coef_name, "exp(coef)"]
      hr_lower <- s$conf.int[coef_name, "lower .95"]
      hr_upper <- s$conf.int[coef_name, "upper .95"]
      p_val <- s$coefficients[coef_name, "Pr(>|z|)"]
      
      level_data <- level_counts %>% filter(as.character(n_lifestyle_cat) == level)
      p_case <- if (nrow(level_data) > 0) level_data$p_case else 0
      
      paf_result <- calc_paf_with_ci(p_case, hr, hr_lower, hr_upper)
      combined_paf <- combined_paf + paf_result$paf
      
      sex_stratified_cumulative[[paste0(sg, "_", level)]] <- data.frame(
        Analysis = "Subgroup: Sex-Stratified Cumulative PAF",
        Sex = sg,
        Outcome = "Overall",
        Outcome_Label = "Overall PPC-MM",
        Level = level,
        HR = round(hr, 3),
        HR_Lower = round(hr_lower, 3),
        HR_Upper = round(hr_upper, 3),
        HR_CI = sprintf("%.2f (%.2f-%.2f)", hr, hr_lower, hr_upper),
        P_value = p_val,
        P_value_fmt = format.pval(p_val, digits = 4),
        Significance = case_when(
          p_val < 0.001 ~ "***",
          p_val < 0.01 ~ "**",
          p_val < 0.05 ~ "*",
          TRUE ~ ""
        ),
        P_case = round(p_case, 4),
        P_case_pct = sprintf("%.1f%%", p_case * 100),
        PAF = round(paf_result$paf, 4),
        PAF_pct = sprintf("%.1f%%", paf_result$paf * 100),
        N = n_total,
        N_Events = n_events,
        stringsAsFactors = FALSE
      )
    }
    
    # Add combined PAF row
    sex_stratified_cumulative[[paste0(sg, "_combined")]] <- data.frame(
      Analysis = "Subgroup: Sex-Stratified Cumulative PAF",
      Sex = sg,
      Outcome = "Overall",
      Outcome_Label = "Overall PPC-MM",
      Level = "Combined (≥1 vs 0)",
      HR = NA,
      HR_Lower = NA,
      HR_Upper = NA,
      HR_CI = "-",
      P_value = NA,
      P_value_fmt = "-",
      Significance = "",
      P_case = NA,
      P_case_pct = "-",
      PAF = round(combined_paf, 4),
      PAF_pct = sprintf("%.1f%%", combined_paf * 100),
      N = n_total,
      N_Events = n_events,
      stringsAsFactors = FALSE
    )
    
    cat("    [OK] Combined PAF:", sprintf("%.1f%%", combined_paf * 100), "\n")
  }
}

# Combine sex-stratified cumulative PAF
if (length(sex_stratified_cumulative) > 0) {
  sex_stratified_cumulative_df <- bind_rows(sex_stratified_cumulative)
} else {
  sex_stratified_cumulative_df <- data.frame()
}

# Save ALWAYS
cat("\n=== Sex-Stratified Cumulative PAF Summary ===\n")
cat("Total rows:", nrow(sex_stratified_cumulative_df), "\n")

tryCatch({
  paf_sex_cum_file <- file.path(OUTPUT_DIR, "PAF_Sex_Stratified_Cumulative.csv")
  if (nrow(sex_stratified_cumulative_df) > 0) {
    summary_cum <- sex_stratified_cumulative_df %>%
      select(Sex, Level, HR_CI, PAF_pct, N_Events) %>%
      arrange(Sex, Level)
    print(summary_cum, n = 50)
    write.csv(sex_stratified_cumulative_df, paf_sex_cum_file, row.names = FALSE)
  } else {
    placeholder <- data.frame(
      Note = "No sex-stratified cumulative PAF results",
      Sex_Groups = "Men, Women"
    )
    write.csv(placeholder, paf_sex_cum_file, row.names = FALSE)
  }
  cat("\n  [OK] Saved: PAF_Sex_Stratified_Cumulative.csv\n")
}, error = function(e) {
  cat("  [ERROR] Failed to save:", e$message, "\n")
})

# Combine sex-stratified results
sex_stratified_all_df <- bind_rows(
  sex_stratified_individual_df,
  sex_stratified_cumulative_df
)

# =============================================================================
# PART 11: Save Results
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   SAVING PAF RESULTS                                            \n")
cat("================================================================\n")

# Combine all results (including age and sex stratified)
all_paf_results <- bind_rows(
  individual_paf_df %>% mutate(Category = "Individual Factors"),
  cumulative_paf_df %>% mutate(Category = "Cumulative Effect"),
  sens2a_df %>% mutate(Category = "S2: Heavy Drink Individual"),
  sens2b_df %>% mutate(Category = "S2: Heavy Drink Cumulative"),
  sens3_individual_df %>% mutate(Category = "S3: MICE Individual"),
  sens3_cumulative_df %>% mutate(Category = "S3: MICE Cumulative"),
  sens4_df %>% mutate(Category = "S4: Drop1st"),
  age_stratified_individual_df %>% mutate(Category = "Subgroup: Age-Stratified Individual"),
  age_stratified_cumulative_df %>% mutate(Category = "Subgroup: Age-Stratified Cumulative"),
  sex_stratified_individual_df %>% mutate(Category = "Subgroup: Sex-Stratified Individual"),
  sex_stratified_cumulative_df %>% mutate(Category = "Subgroup: Sex-Stratified Cumulative")
)

cat("Total PAF results:", nrow(all_paf_results), "rows\n")

# Analysis description
paf_description <- data.frame(
  Analysis = c(
    "Primary: Individual PAF",
    "Primary: Cumulative PAF (4-level)",
    "S2: Heavy Drink Individual",
    "S2: Heavy Drink Cumulative",
    "S3: MICE Individual",
    "S3: MICE Cumulative",
    "S4: Drop1st",
    "Subgroup: Age-Stratified Individual PAF",
    "Subgroup: Age-Stratified Cumulative PAF",
    "Subgroup: Sex-Stratified Individual PAF",
    "Subgroup: Sex-Stratified Cumulative PAF"
  ),
  Description = c(
    "PAF for each lifestyle factor (drink, smoke, PA, social) from mutually adjusted HR",
    "PAF for cumulative categories (0/1/2/3+); Combined PAF = sum of level-specific PAFs",
    "Sensitivity: PAF using heavy drinking instead of any drinking",
    "Sensitivity: Cumulative PAF using heavy drinking score",
    "Sensitivity: PAF using MICE imputed data",
    "Sensitivity: Cumulative PAF using MICE data",
    "Sensitivity: PAF excluding first follow-up wave (reverse causality)",
    "Subgroup: PAF for individual factors stratified by age groups (50-59, 60-69, 70-79, 80+)",
    "Subgroup: Cumulative PAF stratified by age groups (50-59, 60-69, 70-79, 80+)",
    "Subgroup: PAF for individual factors stratified by sex (Men, Women)",
    "Subgroup: Cumulative PAF stratified by sex (Men, Women)"
  ),
  Interpretation = c(
    "Proportion of cases attributable to each specific lifestyle factor",
    "Combined PAF = proportion of cases preventable if all achieved 0 unhealthy lifestyles",
    "Tests robustness when using stricter drinking definition",
    "Tests robustness of cumulative PAF with stricter drinking definition",
    "Tests impact of missing data on PAF estimates",
    "Tests impact of missing data on cumulative PAF",
    "Tests robustness excluding early events that may represent reverse causation",
    "Examines if lifestyle-attributable risk varies by age",
    "Examines if population-level impact varies by age",
    "Examines if lifestyle-attributable risk varies by sex",
    "Examines if population-level impact varies by sex"
  ),
  stringsAsFactors = FALSE
)

# Save to Excel
paf_excel <- list(
  "Analysis_Description" = paf_description,
  "Primary_Individual" = individual_paf_df,
  "Primary_Cumulative" = cumulative_paf_df,
  "S2_HeavyDrink_Individual" = sens2a_df,
  "S2_HeavyDrink_Cumulative" = sens2b_df,
  "S3_MICE_Individual" = sens3_individual_df,
  "S3_MICE_Cumulative" = sens3_cumulative_df,
  "S4_Drop1st" = sens4_df,
  "Age_Stratified_Individual" = age_stratified_individual_df,
  "Age_Stratified_Cumulative" = age_stratified_cumulative_df,
  "Sex_Stratified_Individual" = sex_stratified_individual_df,
  "Sex_Stratified_Cumulative" = sex_stratified_cumulative_df,
  "All_Results" = all_paf_results
)

# ===== Save Excel and CSV with robust error handling =====
cat("\n--- Saving PAF Excel and CSV Files ---\n")

# Save Excel with tryCatch
tryCatch({
  cat("  Saving Excel file...\n")
  filepath_xlsx <- file.path(OUTPUT_DIR, "PAF_Analysis_Comprehensive.xlsx")
  writexl::write_xlsx(paf_excel, filepath_xlsx)
  cat("  [OK] Saved: PAF_Analysis_Comprehensive.xlsx\n")
}, error = function(e) {
  cat("  [ERROR] Excel save failed:", e$message, "\n")
})

# Save main CSV
tryCatch({
  cat("  Saving main CSV file...\n")
  if (nrow(all_paf_results) > 0) {
    filepath_csv <- file.path(OUTPUT_DIR, "PAF_Analysis_All_Results.csv")
    write.csv(all_paf_results, filepath_csv, row.names = FALSE)
    cat("  [OK] Saved: PAF_Analysis_All_Results.csv (", nrow(all_paf_results), " rows)\n")
  } else {
    cat("  [WARN] all_paf_results is empty, skipping main CSV\n")
  }
}, error = function(e) {
  cat("  [ERROR] Main CSV save failed:", e$message, "\n")
})

# ===== Save Individual CSV Files =====
cat("\n--- Saving Individual PAF CSV Files (with error handling) ---\n")

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

# 1. Primary Individual Factors PAF
safe_save_csv(individual_paf_df, "PAF_Primary_Individual.csv", "Primary Individual PAF")

# 2. Primary Cumulative PAF
safe_save_csv(cumulative_paf_df, "PAF_Primary_Cumulative.csv", "Primary Cumulative PAF")

# 3. Sensitivity S2 Heavy Drink
safe_save_csv(sens2a_df, "PAF_S2_HeavyDrink_Individual.csv", "S2a Heavy Individual PAF")
safe_save_csv(sens2b_df, "PAF_S2_HeavyDrink_Cumulative.csv", "S2b Heavy Cumulative PAF")

# 4. Sensitivity S3 MICE
safe_save_csv(sens3_individual_df, "PAF_S3_MICE_Individual.csv", "S3 MICE Individual PAF")
safe_save_csv(sens3_cumulative_df, "PAF_S3_MICE_Cumulative.csv", "S3 MICE Cumulative PAF")

# 5. Sensitivity S4 Drop1st
safe_save_csv(sens4_df, "PAF_S4_Drop1st.csv", "S4 Drop1st PAF")

# 6. Age-Stratified PAF (IMPORTANT: must save)
if (exists("age_stratified_individual_df")) {
  safe_save_csv(age_stratified_individual_df, "PAF_Age_Stratified_Individual.csv", "Age Stratified Individual PAF")
} else {
  cat("  [WARN] age_stratified_individual_df not found\n")
}
if (exists("age_stratified_cumulative_df")) {
  safe_save_csv(age_stratified_cumulative_df, "PAF_Age_Stratified_Cumulative.csv", "Age Stratified Cumulative PAF")
} else {
  cat("  [WARN] age_stratified_cumulative_df not found\n")
}

# 7. Sex-Stratified PAF (IMPORTANT: must save)
if (exists("sex_stratified_individual_df")) {
  safe_save_csv(sex_stratified_individual_df, "PAF_Sex_Stratified_Individual.csv", "Sex Stratified Individual PAF")
} else {
  cat("  [WARN] sex_stratified_individual_df not found\n")
}
if (exists("sex_stratified_cumulative_df")) {
  safe_save_csv(sex_stratified_cumulative_df, "PAF_Sex_Stratified_Cumulative.csv", "Sex Stratified Cumulative PAF")
} else {
  cat("  [WARN] sex_stratified_cumulative_df not found\n")
}

# ===== Create Summary Table for Publication =====
paf_summary <- data.frame()

# Individual PAF summary
if (nrow(individual_paf_df) > 0) {
  indiv_summ <- individual_paf_df %>%
    filter(Outcome == "Overall") %>%
    select(Variable_Label, HR, HR_Lower, HR_Upper, PAF, PAF_pct, N_Events) %>%
    mutate(Analysis = "Primary: Individual PAF")
  paf_summary <- bind_rows(paf_summary, indiv_summ)
}

# Cumulative PAF summary
if (nrow(cumulative_paf_df) > 0) {
  cum_summ <- cumulative_paf_df %>%
    filter(Outcome == "Overall") %>%
    select(Level, HR, HR_Lower, HR_Upper, PAF, PAF_pct, N_Events) %>%
    rename(Variable_Label = Level) %>%
    mutate(Analysis = "Primary: Cumulative PAF")
  paf_summary <- bind_rows(paf_summary, cum_summ)
}

tryCatch({
  if (nrow(paf_summary) > 0) {
    filepath <- file.path(OUTPUT_DIR, "PAF_Summary_Overall.csv")
    write.csv(paf_summary, filepath, row.names = FALSE)
    cat("  [OK] Saved: PAF_Summary_Overall.csv (", nrow(paf_summary), " rows)\n")
  }
}, error = function(e) {
  cat("  [ERROR] Failed to save PAF summary:", e$message, "\n")
})

# Save RDS with tryCatch
tryCatch({
  saveRDS(list(
    primary_individual = individual_paf_df,
    primary_cumulative = cumulative_paf_df,
    sens2a = sens2a_df,
    sens2b = sens2b_df,
    sens3_individual = sens3_individual_df,
    sens3_cumulative = sens3_cumulative_df,
    sens4 = sens4_df,
    all = all_paf_results
  ), file.path(OUTPUT_DIR, "PAF_results.rds"))
  cat("  [OK] Saved: PAF_results.rds\n")
}, error = function(e) {
  cat("  [ERROR] Failed to save PAF RDS:", e$message, "\n")
})

# =============================================================================
# PART 11: Summary Output
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   PAF ANALYSIS SUMMARY                                          \n")
cat("================================================================\n")

# Individual factors summary
cat("\n--- PRIMARY: Individual Factors PAF ---\n")
if (nrow(individual_paf_df) > 0) {
  indiv_summary <- individual_paf_df %>%
    group_by(Outcome, Outcome_Type) %>%
    summarise(
      Max_PAF_Factor = Variable_Label[which.max(PAF)],
      Max_PAF = max(PAF_pct),
      .groups = "drop"
    )
  print(indiv_summary)
}

# Cumulative summary
cat("\n--- PRIMARY: Cumulative PAF (Combined ?? vs 0) ---\n")
if (nrow(cumulative_paf_df) > 0) {
  cum_summary <- cumulative_paf_df %>%
    filter(Level == "Combined (?? vs 0)") %>%
    select(Outcome, Outcome_Type, PAF_pct, N_Events)
  print(cum_summary)
}

# Sensitivity summaries
cat("\n--- SENSITIVITY ANALYSES ---\n")
cat("S2 (Heavy Drink - Individual):", nrow(sens2a_df), "results\n")
cat("S2 (Heavy Drink - Cumulative):", nrow(sens2b_df), "results\n")
cat("S3 (MICE - Individual):", nrow(sens3_individual_df), "results\n")
cat("S3 (MICE - Cumulative):", nrow(sens3_cumulative_df), "results\n")
cat("S4 (Drop1st):", nrow(sens4_df), "results\n")

# =============================================================================
# COMPLETE
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   PAF Analysis Complete                                         \n")
cat("================================================================\n")
cat("Sample size:", nrow(pooled_data), "\n")
cat("Total PAF results:", nrow(all_paf_results), "\n")

# List all generated files
cat("\n--- Generated Output Files ---\n")
cat("Directory:", OUTPUT_DIR, "\n")

# Check which files were actually created
paf_files <- c(
  "PAF_Analysis_Comprehensive.xlsx",
  "PAF_Analysis_All_Results.csv",
  "PAF_results.rds",
  "PAF_Primary_Individual.csv",
  "PAF_Primary_Cumulative.csv",
  "PAF_S2_HeavyDrink_Individual.csv",
  "PAF_S2_HeavyDrink_Cumulative.csv",
  "PAF_S3_MICE_Individual.csv",
  "PAF_S3_MICE_Cumulative.csv",
  "PAF_S4_Drop1st.csv",
  "PAF_Age_Stratified_Individual.csv",
  "PAF_Age_Stratified_Cumulative.csv",
  "PAF_Summary_Overall.csv"
)

for (f in paf_files) {
  filepath <- file.path(OUTPUT_DIR, f)
  if (file.exists(filepath)) {
    cat("  [✓]", f, "\n")
  } else {
    cat("  [✗]", f, "(NOT CREATED)\n")
  }
}
