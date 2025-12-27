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

# Load and pool data (reuse function from 08)
load_and_pool_data <- function(data_paths, data_type = "main") {
  
  pooled_list <- list()
  
  for (coh in names(data_paths)) {
    file_path <- data_paths[[coh]]
    
    if (is.null(file_path) || !file.exists(file_path)) next
    
    df <- read.csv(file_path, stringsAsFactors = FALSE)
    bw <- BASELINE_WAVES[[coh]]
    
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
    
    pooled_list[[coh]] <- df_std
  }
  
  bind_rows(pooled_list)
}

# Load main data
pooled_main <- load_and_pool_data(DATA_PATHS, "main")

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
covariates <- c(COX_COVARIATES$pooled, "cohort")

cat("\nPAF Analysis Covariates:", paste(covariates, collapse = " + "), "\n")

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

mice_files_exist <- sapply(MICE_PATHS, file.exists)
cat("MICE data available:", sum(mice_files_exist), "/", length(mice_files_exist), "cohorts\n")

sens3_individual <- list()
sens3_cumulative <- list()

if (sum(mice_files_exist) > 0) {
  pooled_mice <- load_and_pool_data(MICE_PATHS[mice_files_exist], "mice")
  
  if (nrow(pooled_mice) > 0) {
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
    
    cat("MICE analysis sample size:", nrow(pooled_mice), "\n")
    
    # Individual factors
    for (out_name in all_outcomes) {
      result <- calculate_individual_paf(pooled_mice, out_name, covariates, use_heavy_drink = FALSE)
      if (!is.null(result)) {
        result$Analysis <- "S3: MICE - Individual PAF"
        sens3_individual[[out_name]] <- result
      }
    }
    
    # Cumulative effect
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
  cat("No MICE data available, skipping S3\n")
}

sens3_individual_df <- bind_rows(sens3_individual)
sens3_cumulative_df <- bind_rows(sens3_cumulative)

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
# PART 10: Save Results
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("   SAVING PAF RESULTS                                            \n")
cat("================================================================\n")

# Combine all results
all_paf_results <- bind_rows(
  individual_paf_df %>% mutate(Category = "Individual Factors"),
  cumulative_paf_df %>% mutate(Category = "Cumulative Effect"),
  sens2a_df %>% mutate(Category = "S2: Heavy Drink Individual"),
  sens2b_df %>% mutate(Category = "S2: Heavy Drink Cumulative"),
  sens3_individual_df %>% mutate(Category = "S3: MICE Individual"),
  sens3_cumulative_df %>% mutate(Category = "S3: MICE Cumulative"),
  sens4_df %>% mutate(Category = "S4: Drop1st")
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
    "S4: Drop1st"
  ),
  Description = c(
    "PAF for each lifestyle factor (drink, smoke, PA, social) from mutually adjusted HR",
    "PAF for cumulative categories (0/1/2/3+); Combined PAF = sum of level-specific PAFs",
    "Sensitivity: PAF using heavy drinking instead of any drinking",
    "Sensitivity: Cumulative PAF using heavy drinking score",
    "Sensitivity: PAF using MICE imputed data",
    "Sensitivity: Cumulative PAF using MICE data",
    "Sensitivity: PAF excluding first follow-up wave (reverse causality)"
  ),
  Interpretation = c(
    "Proportion of cases attributable to each specific lifestyle factor",
    "Combined PAF = proportion of cases preventable if all achieved 0 unhealthy lifestyles",
    "Tests robustness when using stricter drinking definition",
    "Tests robustness of cumulative PAF with stricter drinking definition",
    "Tests impact of missing data on PAF estimates",
    "Tests impact of missing data on cumulative PAF",
    "Tests robustness excluding early events that may represent reverse causation"
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
  "All_Results" = all_paf_results
)

save_to_excel(paf_excel, "PAF_Analysis_Comprehensive.xlsx")
save_to_csv(all_paf_results, "PAF_Analysis_All_Results.csv")

# Save RDS
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
cat("\nOutputs saved to:", OUTPUT_DIR, "\n")
cat("  - PAF_Analysis_Comprehensive.xlsx\n")
cat("  - PAF_Analysis_All_Results.csv\n")
cat("  - PAF_results.rds\n")
