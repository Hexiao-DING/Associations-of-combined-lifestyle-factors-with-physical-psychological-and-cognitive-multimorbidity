###############################################################################
# 03_Pooled_Cox_Analysis.R
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
    
    # Rename if variables exist
    if (score_var %in% names(df_std)) df_std <- df_std %>% rename(unhealthy_score = !!score_var)
    if (drink_var %in% names(df_std)) df_std <- df_std %>% rename(unhealthy_drink = !!drink_var)
    if (smoke_var %in% names(df_std)) df_std <- df_std %>% rename(unhealthy_smoke = !!smoke_var)
    if (pa_var %in% names(df_std)) df_std <- df_std %>% rename(unhealthy_pa = !!pa_var)
    if (soc_var %in% names(df_std)) df_std <- df_std %>% rename(unhealthy_soc = !!soc_var)
    if (cat_var %in% names(df_std)) df_std <- df_std %>% rename(unhealthy_cat = !!cat_var)
    if (heavy_drink_var %in% names(df_std)) df_std <- df_std %>% rename(heavy_drink = !!heavy_drink_var)
    if (score_heavy_var %in% names(df_std)) df_std <- df_std %>% rename(unhealthy_score_heavy = !!score_heavy_var)
    
    # Select common variables for pooling
    common_vars <- c("cohort", "country", "age_baseline", "sex", "edu",
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
  
  # Complete case filter
  df <- df %>%
    filter(!is.na(age_baseline),
           !is.na(sex),
           !is.na(edu),
           !is.na(unhealthy_drink),
           !is.na(unhealthy_smoke),
           !is.na(unhealthy_pa),
           !is.na(unhealthy_soc))
  
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
#' @return Data frame with results
run_cox_model <- function(data, outcome_name, exposure_var, covariates, 
                          model_label, calc_trend = FALSE) {
  
  out_info <- outcome_info[[outcome_name]]
  event_var <- out_info$event
  time_var <- out_info$time
  
  # Check if variables exist
  if (!event_var %in% names(data) || !time_var %in% names(data)) {
    return(NULL)
  }
  
  # Filter valid observations
  df <- data %>%
    filter(!is.na(.data[[event_var]]) & !is.na(.data[[time_var]]) & .data[[time_var]] > 0)
  
  n_total <- nrow(df)
  n_events <- sum(df[[event_var]], na.rm = TRUE)
  
  if (n_events < 10) {
    return(NULL)
  }
  
  # Build formula
  cov_str <- paste(covariates, collapse = " + ")
  formula_str <- paste0("Surv(", time_var, ", ", event_var, ") ~ ", 
                        exposure_var, " + ", cov_str)
  
  # Fit model
  fit <- tryCatch({
    coxph(as.formula(formula_str), data = df)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(fit)) return(NULL)
  
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
      trend_formula <- paste0("Surv(", time_var, ", ", event_var, ") ~ ", 
                              trend_var, " + ", cov_str)
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
#' @return Data frame with results
run_individual_factors <- function(data, outcome_name, covariates, use_heavy_drink = FALSE) {
  
  out_info <- outcome_info[[outcome_name]]
  event_var <- out_info$event
  time_var <- out_info$time
  
  if (!event_var %in% names(data) || !time_var %in% names(data)) {
    return(NULL)
  }
  
  df <- data %>%
    filter(!is.na(.data[[event_var]]) & !is.na(.data[[time_var]]) & .data[[time_var]] > 0)
  
  n_events <- sum(df[[event_var]], na.rm = TRUE)
  
  if (n_events < 10) return(NULL)
  
  # Select drink variable
  drink_var <- ifelse(use_heavy_drink && "heavy_drink" %in% names(df), 
                      "heavy_drink", "unhealthy_drink")
  
  lifestyle_vars <- c(drink_var, "unhealthy_smoke", "unhealthy_pa", "unhealthy_soc")
  
  # Build formula with all 4 factors mutually adjusted
  formula_str <- paste0("Surv(", time_var, ", ", event_var, ") ~ ",
                        paste(lifestyle_vars, collapse = " + "), " + ",
                        paste(covariates, collapse = " + "))
  
  fit <- tryCatch({
    coxph(as.formula(formula_str), data = df)
  }, error = function(e) NULL)
  
  if (is.null(fit)) return(NULL)
  
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

# Load main data
pooled_main <- load_and_pool_data(DATA_PATHS, "main")
pooled_analysis <- prepare_analysis_data(pooled_main, use_heavy_drink = TRUE)

cat("\nComplete case sample size:", nrow(pooled_analysis), "\n")

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
# PART 6: PRIMARY ANALYSIS - Individual Lifestyle Factors
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   PRIMARY ANALYSIS A: Individual Lifestyle Factors              \n")
cat("   (4 factors mutually adjusted)                                 \n")
cat("================================================================\n")

covariates <- c("age_baseline", "sex", "edu", "cohort")
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
# PART 10: SENSITIVITY ANALYSIS S3 - MICE Data
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   SENSITIVITY ANALYSIS S3: MICE Imputed Data                    \n")
cat("================================================================\n")

# Check if MICE data exists
mice_files_exist <- sapply(MICE_PATHS, file.exists)
cat("MICE data available:", sum(mice_files_exist), "/", length(mice_files_exist), "cohorts\n")

sens3_individual <- list()
sens3_cumulative <- list()

if (sum(mice_files_exist) > 0) {
  pooled_mice <- load_and_pool_data(MICE_PATHS[mice_files_exist], "mice")
  
  if (nrow(pooled_mice) > 0) {
    pooled_mice_analysis <- prepare_analysis_data(pooled_mice, use_heavy_drink = FALSE)
    
    cat("MICE analysis sample size:", nrow(pooled_mice_analysis), "\n")
    
    # Individual factors
    for (out_name in all_outcomes) {
      result <- run_individual_factors(pooled_mice_analysis, out_name, covariates, use_heavy_drink = FALSE)
      if (!is.null(result)) {
        result$Model <- "S3: MICE - Individual Factors"
        sens3_individual[[out_name]] <- result
      }
    }
    
    # Cumulative effect
    for (out_name in all_outcomes) {
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
      }
    }
  }
} else {
  cat("No MICE data available, skipping S3\n")
}

sens3_individual_df <- bind_rows(sens3_individual)
sens3_cumulative_df <- bind_rows(sens3_cumulative)

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
# PART 12: Dose-Response Plots
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

# Combine all results
all_results <- bind_rows(
  individual_results_df %>% mutate(Analysis = "Primary: Individual Factors"),
  cumulative_results_df %>% mutate(Analysis = "Primary: Cumulative (4-level)"),
  sens1_results_df %>% mutate(Analysis = "S1: 5-level Categories"),
  sens2a_results_df %>% mutate(Analysis = "S2a: Heavy Drink - Individual"),
  sens2b_results_df %>% mutate(Analysis = "S2b: Heavy Drink - Cumulative"),
  sens3_individual_df %>% mutate(Analysis = "S3: MICE - Individual"),
  sens3_cumulative_df %>% mutate(Analysis = "S3: MICE - Cumulative"),
  sens4_results_df %>% mutate(Analysis = "S4: Drop1st")
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

# Prepare Excel output
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
  "All_Results" = all_results
)

# Save
save_to_excel(results_for_excel, "Pooled_Cox_Results_Comprehensive.xlsx")
save_to_csv(all_results, "Pooled_Cox_All_Results.csv")

# Save model objects
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

saveRDS(cox_results, file.path(OUTPUT_DIR, "Pooled_cox_results.rds"))

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
cat("\nOutputs saved to:", OUTPUT_DIR, "\n")
cat("  - Pooled_Cox_Results_Comprehensive.xlsx\n")
cat("  - Pooled_Cox_All_Results.csv\n")
cat("  - Pooled_cox_results.rds\n")
cat("  - DoseResponse_*.png/pdf\n")
