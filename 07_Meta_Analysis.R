###############################################################################
# 07_Meta_Analysis.R
# ============================================================================
# Project: Lifestyle Factors and PPC-MM Association Study
# Purpose: Comprehensive Meta-analysis (Lancet Standard)
# ============================================================================
#
# ANALYSIS STRUCTURE:
#
# === PRIMARY ANALYSIS (Main Data) ===
#   A. Individual Lifestyle Factors (mutually adjusted) - Meta across cohorts
#   B. Cumulative Effect (4-level: 0/1/2/3+) - Meta across cohorts
#   Outcomes: Overall PPC-MM (Primary) + 4 subtypes (Secondary)
#
# === SENSITIVITY ANALYSES ===
#   S1: 5-level categories (0/1/2/3/4)
#   S2: Heavy drinking definition
#   S3: MICE imputed data
#   S4: Drop1st (exclude first follow-up wave)
#
# STATISTICAL METHODS:
#   - Fixed-effects model: Inverse variance method
#   - Random-effects model: DerSimonian-Laird method
#   - Heterogeneity: Q statistic, I², ?²
#   - Leave-one-out sensitivity analysis
#   - Funnel plot and Egger's test for publication bias
#
###############################################################################

# =============================================================================
# PART 1: ENVIRONMENT SETUP
# =============================================================================

cat("\n")
cat("###############################################################################\n")
cat("#                                                                             #\n")
cat("#          META-ANALYSIS: COMPREHENSIVE LANCET STANDARD                      #\n")
cat("#                                                                             #\n")
cat("###############################################################################\n")
cat("\n")

source("C:/Users/user/Desktop/Scientific Research/UKB Project/Jung Sun Lab/Project1_Lifestyle_PPM_Hexiao,Hongtao/Code/00_Functions_and_Setup.R")

# Install and load meta-analysis packages
required_meta_pkgs <- c("meta", "metafor", "grid", "gridExtra")

for (pkg in required_meta_pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

cat("[INFO] All required packages loaded.\n\n")

# =============================================================================
# PART 2: DEFINITIONS
# =============================================================================

# Outcome definitions (Overall is PRIMARY)
outcome_info <- list(
  Overall = list(
    event = "event_ppcmm", time = "time_ppcmm_months",
    label = "Overall PPC-MM", type = "Primary"
  ),
  P1P2 = list(
    event = "event_mm_phys_psych", time = "time_mm_phys_psych",
    label = "Physical-Psychological (P1P2)", type = "Secondary"
  ),
  P1C = list(
    event = "event_mm_phys_cog", time = "time_mm_phys_cog",
    label = "Physical-Cognitive (P1C)", type = "Secondary"
  ),
  P2C = list(
    event = "event_mm_psych_cog", time = "time_mm_psych_cog",
    label = "Psychological-Cognitive (P2C)", type = "Secondary"
  ),
  P1P2C = list(
    event = "event_mm_all_three", time = "time_mm_all_three",
    label = "All Three (P1P2C)", type = "Secondary"
  )
)

all_outcomes <- names(outcome_info)

# Cohort info
cohort_info <- list(
  CHARLS = list(wave = BASELINE_WAVES$CHARLS, country = COHORT_COUNTRIES$CHARLS),
  ELSA = list(wave = BASELINE_WAVES$ELSA, country = COHORT_COUNTRIES$ELSA),
  HRS = list(wave = BASELINE_WAVES$HRS, country = COHORT_COUNTRIES$HRS),
  SHARE = list(wave = BASELINE_WAVES$SHARE, country = COHORT_COUNTRIES$SHARE),
  MHAS = list(wave = BASELINE_WAVES$MHAS, country = COHORT_COUNTRIES$MHAS)
)

# Create output directories
fig_dir <- file.path(OUTPUT_DIR, "Figures", "Forest")
funnel_dir <- file.path(OUTPUT_DIR, "Figures", "Funnel")
sensitivity_dir <- file.path(OUTPUT_DIR, "Figures", "Sensitivity")

for (dir_path in c(fig_dir, funnel_dir, sensitivity_dir)) {
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
}

# =============================================================================
# PART 3: HELPER FUNCTIONS
# =============================================================================

close_all_devices <- function() {
  while (dev.cur() > 1) tryCatch(dev.off(), error = function(e) NULL)
}

save_plot_safe <- function(plot_expr, filepath_base, width = 10, height = 6, res = 300) {
  close_all_devices()
  results <- list(pdf = FALSE, png = FALSE)
  
  # PDF
  filepath_pdf <- paste0(filepath_base, ".pdf")
  tryCatch({
    pdf(file = filepath_pdf, width = width, height = height, useDingbats = FALSE)
    eval(plot_expr)
    dev.off()
    Sys.sleep(0.1)
    if (file.exists(filepath_pdf) && file.info(filepath_pdf)$size > 500) {
      results$pdf <- TRUE
    }
  }, error = function(e) close_all_devices())
  
  # PNG
  filepath_png <- paste0(filepath_base, ".png")
  tryCatch({
    png(filename = filepath_png, width = width, height = height, units = "in", res = res)
    eval(plot_expr)
    dev.off()
    Sys.sleep(0.1)
    if (file.exists(filepath_png) && file.info(filepath_png)$size > 1000) {
      results$png <- TRUE
    }
  }, error = function(e) close_all_devices())
  
  return(results)
}

extract_hr_from_cox <- function(fit) {
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  data.frame(
    variable = rownames(s$coefficients),
    log_hr = s$coefficients[, "coef"],
    se_log_hr = s$coefficients[, "se(coef)"],
    hr = s$coefficients[, "exp(coef)"],
    hr_lower = s$conf.int[, "lower .95"],
    hr_upper = s$conf.int[, "upper .95"],
    p_value = s$coefficients[, "Pr(>|z|)"],
    n = fit$n,
    n_events = fit$nevent,
    stringsAsFactors = FALSE
  )
}

run_cohort_cox <- function(df, outcome_name, exposure_var, covariates, 
                            event_override = NULL, time_override = NULL) {
  
  event_var <- ifelse(is.null(event_override), outcome_info[[outcome_name]]$event, event_override)
  time_var <- ifelse(is.null(time_override), outcome_info[[outcome_name]]$time, time_override)
  
  if (!event_var %in% names(df) || !time_var %in% names(df)) return(NULL)
  
  df_valid <- df %>%
    filter(!is.na(.data[[event_var]]) & !is.na(.data[[time_var]]) & .data[[time_var]] > 0)
  
  n_events <- sum(df_valid[[event_var]], na.rm = TRUE)
  if (n_events < 10) return(NULL)
  
  formula_str <- paste0("Surv(", time_var, ", ", event_var, ") ~ ", 
                        exposure_var, " + ", paste(covariates, collapse = " + "))
  
  tryCatch(coxph(as.formula(formula_str), data = df_valid), error = function(e) NULL)
}

run_meta_analysis <- function(meta_data, title = "") {
  if (nrow(meta_data) < 2) return(NULL)
  tryCatch({
    metagen(
      TE = meta_data$log_hr, seTE = meta_data$se_log_hr, studlab = meta_data$study,
      sm = "HR", method.tau = "DL", common = TRUE, random = TRUE, prediction = TRUE, title = title
    )
  }, error = function(e) NULL)
}

extract_ma_summary <- function(ma, analysis, outcome, exposure_level) {
  if (is.null(ma)) return(NULL)
  
  # Handle outcome type lookup safely
  outcome_base <- sub("_drop1st$", "", outcome)  # Remove _drop1st suffix if present
  outcome_type <- if (!is.null(outcome_info[[outcome_base]])) {
    outcome_info[[outcome_base]]$type
  } else {
    "Sensitivity"
  }
  
  data.frame(
    Analysis = analysis, Outcome = outcome, Exposure_Level = exposure_level,
    Outcome_Type = outcome_type,
    N_Studies = ma$k,
    Fixed_HR = round(exp(ma$TE.common), 3),
    Fixed_Lower = round(exp(ma$lower.common), 3),
    Fixed_Upper = round(exp(ma$upper.common), 3),
    Fixed_P = format.pval(ma$pval.common, digits = 3),
    Random_HR = round(exp(ma$TE.random), 3),
    Random_Lower = round(exp(ma$lower.random), 3),
    Random_Upper = round(exp(ma$upper.random), 3),
    Random_P = format.pval(ma$pval.random, digits = 3),
    Q = round(ma$Q, 2), I2 = round(ma$I2 * 100, 1), Tau2 = round(ma$tau2, 4),
    Recommended = ifelse(ma$I2 > 0.5, "Random", "Fixed"),
    stringsAsFactors = FALSE
  )
}

run_leave_one_out <- function(ma, meta_data) {
  if (is.null(ma) || ma$k < 3) return(NULL)
  loo_results <- list()
  for (i in seq_len(ma$k)) {
    ma_loo <- tryCatch({
      metagen(TE = meta_data$log_hr[-i], seTE = meta_data$se_log_hr[-i],
              studlab = meta_data$study[-i], sm = "HR", method.tau = "DL")
    }, error = function(e) NULL)
    if (!is.null(ma_loo)) {
      loo_results[[i]] <- data.frame(
        Excluded = meta_data$study[i],
        Random_HR = round(exp(ma_loo$TE.random), 3),
        Random_Lower = round(exp(ma_loo$lower.random), 3),
        Random_Upper = round(exp(ma_loo$upper.random), 3),
        I2 = round(ma_loo$I2 * 100, 1)
      )
    }
  }
  bind_rows(loo_results)
}

run_eggers_test <- function(ma) {
  if (is.null(ma) || ma$k < 3) return(NULL)
  egger <- tryCatch(metabias(ma, method.bias = "linreg", k.min = 3), error = function(e) NULL)
  if (!is.null(egger)) {
    data.frame(
      Intercept = round(egger$estimate[1], 3),
      P_value = format.pval(egger$p.value, digits = 3),
      Interpretation = ifelse(egger$p.value < 0.1, "Potential bias", "No evidence")
    )
  }
}

# =============================================================================
# PART 4: FIT COHORT-SPECIFIC COX MODELS
# =============================================================================

cat("=============================================================================\n")
cat("PART 4: FITTING COHORT-SPECIFIC COX MODELS\n")
cat("=============================================================================\n\n")

all_study_results <- list()

for (coh in names(cohort_info)) {
  
  cat("--------------------------------------------------------------------------------\n")
  cat("Processing:", coh, "\n")
  cat("--------------------------------------------------------------------------------\n")
  
  data_file <- DATA_PATHS[[coh]]
  if (is.null(data_file) || !file.exists(data_file)) {
    cat("  [SKIP] Data file not found\n\n")
    next
  }
  
  df <- read.csv(data_file, stringsAsFactors = FALSE)
  bw <- cohort_info[[coh]]$wave
  
  # Variable names
  score_var <- paste0("w", bw, "_unhealthy_score")
  heavy_score_var <- paste0("w", bw, "_unhealthy_score_heavy")
  heavy_drink_var <- paste0("w", bw, "_heavy_drink")
  
  if (!score_var %in% names(df)) {
    cat("  [SKIP] Required variable not found\n\n")
    next
  }
  
  # Prepare data
  df <- df %>%
    filter(!is.na(age_baseline) & !is.na(sex) & !is.na(edu)) %>%
    mutate(
      sex = factor(sex), edu = factor(edu),
      n_lifestyle_cat = factor(case_when(
        .data[[score_var]] == 0 ~ "0", .data[[score_var]] == 1 ~ "1",
        .data[[score_var]] == 2 ~ "2", .data[[score_var]] >= 3 ~ "3+"
      ), levels = c("0", "1", "2", "3+")),
      n_lifestyle_5cat = factor(.data[[score_var]], levels = 0:4)
    )
  
  # Heavy drink variables if available
  if (heavy_score_var %in% names(df)) {
    df <- df %>%
      mutate(n_lifestyle_cat_heavy = factor(case_when(
        .data[[heavy_score_var]] == 0 ~ "0", .data[[heavy_score_var]] == 1 ~ "1",
        .data[[heavy_score_var]] == 2 ~ "2", .data[[heavy_score_var]] >= 3 ~ "3+"
      ), levels = c("0", "1", "2", "3+")))
  }
  
  cat("  N =", nrow(df), "\n")
  
  # Use centralized covariate configuration
  # For cohort-specific analysis: age_baseline, sex, edu
  # For SHARE (multi-country): add region
  if (coh == "SHARE" && "region" %in% names(df)) {
    covariates <- COX_COVARIATES$share_specific
    cat("  Covariates (SHARE):", paste(covariates, collapse = " + "), "\n")
  } else {
    covariates <- COX_COVARIATES$cohort_specific
    cat("  Covariates:", paste(covariates, collapse = " + "), "\n")
  }
  
  # ---- Primary Analysis: All outcomes → 4-level ----
  cat("\n  [Primary] All Outcomes → 4-level:\n")
  
  for (out_name in all_outcomes) {
    fit <- run_cohort_cox(df, out_name, "n_lifestyle_cat", covariates)
    if (!is.null(fit)) {
      hr_results <- extract_hr_from_cox(fit)
      exp_rows <- grepl("n_lifestyle_cat", hr_results$variable)
      cat("    ", out_name, ": variables =", paste(hr_results$variable[exp_rows], collapse = ", "), "\n")
      if (sum(exp_rows) > 0) {
        hr_results <- hr_results[exp_rows, ]
        hr_results$study <- coh
        hr_results$outcome <- out_name
        hr_results$analysis <- "Primary_4Level"
        all_study_results[[paste0("Primary_", out_name, "_", coh)]] <- hr_results
        cat("      -> Added ", nrow(hr_results), " rows\n")
      }
    } else {
      cat("    ", out_name, ": Cox model failed or insufficient events\n")
    }
  }
  
  # ---- S1: 5-level ----
  cat("\n  [S1] 5-level:\n")
  for (out_name in all_outcomes) {
    fit <- run_cohort_cox(df, out_name, "n_lifestyle_5cat", covariates)
    if (!is.null(fit)) {
      hr_results <- extract_hr_from_cox(fit)
      exp_rows <- grepl("n_lifestyle_5cat", hr_results$variable)
      if (sum(exp_rows) > 0) {
        hr_results <- hr_results[exp_rows, ]
        hr_results$study <- coh
        hr_results$outcome <- out_name
        hr_results$analysis <- "S1_5Level"
        all_study_results[[paste0("S1_", out_name, "_", coh)]] <- hr_results
        cat("    ", out_name, ": OK\n")
      }
    }
  }
  
  # ---- S2: Heavy drink ----
  if ("n_lifestyle_cat_heavy" %in% names(df)) {
    cat("\n  [S2] Heavy Drink:\n")
    for (out_name in all_outcomes) {
      fit <- run_cohort_cox(df, out_name, "n_lifestyle_cat_heavy", covariates)
      if (!is.null(fit)) {
        hr_results <- extract_hr_from_cox(fit)
        exp_rows <- grepl("n_lifestyle_cat_heavy", hr_results$variable)
        if (sum(exp_rows) > 0) {
          hr_results <- hr_results[exp_rows, ]
          hr_results$study <- coh
          hr_results$outcome <- out_name
          hr_results$analysis <- "S2_HeavyDrink"
          all_study_results[[paste0("S2_", out_name, "_", coh)]] <- hr_results
          cat("    ", out_name, ": OK\n")
        }
      }
    }
  }
  
  # ---- S4: Drop1st ----
  if ("event_ppcmm_drop1st" %in% names(df)) {
    cat("\n  [S4] Drop1st:\n")
    fit <- run_cohort_cox(df, "Overall", "n_lifestyle_cat", covariates,
                          event_override = "event_ppcmm_drop1st",
                          time_override = "time_ppcmm_drop1st")
    if (!is.null(fit)) {
      hr_results <- extract_hr_from_cox(fit)
      exp_rows <- grepl("n_lifestyle_cat", hr_results$variable)
      if (sum(exp_rows) > 0) {
        hr_results <- hr_results[exp_rows, ]
        hr_results$study <- coh
        hr_results$outcome <- "Overall_drop1st"
        hr_results$analysis <- "S4_Drop1st"
        all_study_results[[paste0("S4_", coh)]] <- hr_results
        cat("    Overall (drop1st): OK\n")
      }
    }
  }
  
  cat("\n")
}

# ---- S3: MICE data ----
cat("--------------------------------------------------------------------------------\n")
cat("Processing: MICE Data\n")
cat("--------------------------------------------------------------------------------\n")

for (coh in names(MICE_PATHS)) {
  mice_file <- MICE_PATHS[[coh]]
  if (!file.exists(mice_file)) next
  
  df <- read.csv(mice_file, stringsAsFactors = FALSE)
  bw <- cohort_info[[coh]]$wave
  score_var <- paste0("w", bw, "_unhealthy_score")
  
  if (!score_var %in% names(df)) next
  
  df <- df %>%
    filter(!is.na(age_baseline) & !is.na(sex) & !is.na(edu)) %>%
    mutate(
      sex = factor(sex), edu = factor(edu),
      n_lifestyle_cat = factor(case_when(
        .data[[score_var]] == 0 ~ "0", .data[[score_var]] == 1 ~ "1",
        .data[[score_var]] == 2 ~ "2", .data[[score_var]] >= 3 ~ "3+"
      ), levels = c("0", "1", "2", "3+"))
    )
  
  # Use centralized covariate configuration
  if (coh == "SHARE" && "region" %in% names(df)) {
    covariates <- COX_COVARIATES$share_specific
  } else {
    covariates <- COX_COVARIATES$cohort_specific
  }
  
  for (out_name in all_outcomes) {
    fit <- run_cohort_cox(df, out_name, "n_lifestyle_cat", covariates)
    if (!is.null(fit)) {
      hr_results <- extract_hr_from_cox(fit)
      exp_rows <- grepl("n_lifestyle_cat", hr_results$variable)
      if (sum(exp_rows) > 0) {
        hr_results <- hr_results[exp_rows, ]
        hr_results$study <- coh
        hr_results$outcome <- out_name
        hr_results$analysis <- "S3_MICE"
        all_study_results[[paste0("S3_", out_name, "_", coh)]] <- hr_results
      }
    }
  }
  cat("  ", coh, ": Processed\n")
}

# Combine results
study_hr_data <- bind_rows(all_study_results)
cat("\n================================================================================\n")
cat("Total study-specific HR results:", nrow(study_hr_data), "\n")
cat("================================================================================\n\n")

# DEBUG: Show data structure
if (nrow(study_hr_data) > 0) {
  cat("[DEBUG] study_hr_data columns:", paste(names(study_hr_data), collapse = ", "), "\n")
  cat("[DEBUG] Unique analyses:", paste(unique(study_hr_data$analysis), collapse = ", "), "\n")
  cat("[DEBUG] Unique outcomes:", paste(unique(study_hr_data$outcome), collapse = ", "), "\n")
  cat("[DEBUG] Sample variables:", paste(head(unique(study_hr_data$variable), 10), collapse = ", "), "\n")
  cat("[DEBUG] Unique studies:", paste(unique(study_hr_data$study), collapse = ", "), "\n\n")
} else {
  cat("[ERROR] study_hr_data is EMPTY! No Cox models succeeded.\n\n")
}

save_to_csv(study_hr_data, "Study_Specific_HR_Results.csv")

# =============================================================================
# PART 5: STUDY CHARACTERISTICS
# =============================================================================

cat("=============================================================================\n")
cat("PART 5: STUDY CHARACTERISTICS\n")
cat("=============================================================================\n\n")

study_characteristics <- data.frame()

for (coh in names(cohort_info)) {
  data_file <- DATA_PATHS[[coh]]
  if (!file.exists(data_file)) next
  
  df <- read.csv(data_file, stringsAsFactors = FALSE)
  bw <- cohort_info[[coh]]$wave
  score_var <- paste0("w", bw, "_unhealthy_score")
  
  study_characteristics <- bind_rows(study_characteristics, data.frame(
    Study = coh,
    Country = cohort_info[[coh]]$country,
    N = nrow(df),
    Pct_Female = round(mean(df$sex == "Women", na.rm = TRUE) * 100, 1),
    Mean_Age = round(mean(df$age_baseline, na.rm = TRUE), 1),
    Mean_Followup = round(mean(df$time_ppcmm_months, na.rm = TRUE), 1),
    Event_Rate = round(mean(df$event_ppcmm, na.rm = TRUE) * 100, 1),
    Mean_Score = round(mean(df[[score_var]], na.rm = TRUE), 2)
  ))
}

cat("Study Characteristics:\n")
print(study_characteristics)

# =============================================================================
# PART 6: RUN META-ANALYSES
# =============================================================================

cat("\n=============================================================================\n")
cat("PART 6: RUNNING META-ANALYSES\n")
cat("=============================================================================\n\n")

ma_results <- list()
ma_summaries <- list()
loo_results <- list()
egger_results <- list()

# Function to run meta for a specific analysis/outcome/level
# FIXED: Renamed 'outcome' parameter to 'outcome_name' to avoid column name conflict
run_meta_for_level <- function(analysis_pattern, var_pattern, analysis_label, outcome_name, level) {
  
  # Use base R subsetting to avoid dplyr column name conflicts
  meta_data <- study_hr_data[
    study_hr_data$analysis == analysis_pattern & 
    study_hr_data$outcome == outcome_name & 
    grepl(var_pattern, study_hr_data$variable, fixed = TRUE), 
  ]
  
  cat("      [DEBUG] Filtering:", analysis_pattern, "/", outcome_name, "/", var_pattern, 
      "-> N =", nrow(meta_data), "\n")
  
  if (nrow(meta_data) < 2) {
    cat("      [SKIP] < 2 studies\n")
    return(NULL)
  }
  
  ma <- run_meta_analysis(meta_data, paste0(outcome_name, " - ", level, " vs 0"))
  if (is.null(ma)) {
    cat("      [SKIP] meta-analysis failed\n")
    return(NULL)
  }
  
  key <- paste0(analysis_label, "_", outcome_name, "_", gsub("\\+", "plus", level))
  ma_results[[key]] <<- ma
  ma_summaries[[key]] <<- extract_ma_summary(ma, analysis_label, outcome_name, level)
  
  loo <- run_leave_one_out(ma, meta_data)
  if (!is.null(loo)) {
    loo$Analysis <- analysis_label
    loo$Outcome <- outcome_name
    loo$Level <- level
    loo_results[[key]] <<- loo
  }
  
  egger <- run_eggers_test(ma)
  if (!is.null(egger)) {
    egger$Analysis <- analysis_label
    egger$Outcome <- outcome_name
    egger$Level <- level
    egger_results[[key]] <<- egger
  }
  
  cat("    ", outcome_name, "(", level, "): Fixed HR =", round(exp(ma$TE.common), 2),
      ", Random HR =", round(exp(ma$TE.random), 2), ", I² =", round(ma$I2 * 100, 1), "%\n")
  
  return(ma)
}

# ---- Primary Analysis ----
cat("--------------------------------------------------------------------------------\n")
cat("PRIMARY ANALYSIS: All Outcomes ? 4-level (0/1/2/3+)\n")
cat("--------------------------------------------------------------------------------\n")

for (out_name in all_outcomes) {
  cat("\n  Outcome:", outcome_info[[out_name]]$label, "(", outcome_info[[out_name]]$type, ")\n")
  for (level in c("1", "2", "3+")) {
    run_meta_for_level("Primary_4Level", paste0("n_lifestyle_cat", level), 
                       "Primary", out_name, level)
  }
}

# ---- S1: 5-level ----
cat("\n--------------------------------------------------------------------------------\n")
cat("S1: 5-level Categories (0/1/2/3/4)\n")
cat("--------------------------------------------------------------------------------\n")

for (out_name in all_outcomes) {
  cat("\n  Outcome:", outcome_info[[out_name]]$label, "\n")
  for (level in c("1", "2", "3", "4")) {
    run_meta_for_level("S1_5Level", paste0("n_lifestyle_5cat", level), 
                       "S1_5Level", out_name, level)
  }
}

# ---- S2: Heavy Drink ----
cat("\n--------------------------------------------------------------------------------\n")
cat("S2: Heavy Drinking Definition\n")
cat("--------------------------------------------------------------------------------\n")

for (out_name in all_outcomes) {
  cat("\n  Outcome:", outcome_info[[out_name]]$label, "\n")
  for (level in c("1", "2", "3+")) {
    run_meta_for_level("S2_HeavyDrink", paste0("n_lifestyle_cat_heavy", level), 
                       "S2_HeavyDrink", out_name, level)
  }
}

# ---- S3: MICE ----
cat("\n--------------------------------------------------------------------------------\n")
cat("S3: MICE Data\n")
cat("--------------------------------------------------------------------------------\n")

for (out_name in all_outcomes) {
  cat("\n  Outcome:", outcome_info[[out_name]]$label, "\n")
  for (level in c("1", "2", "3+")) {
    run_meta_for_level("S3_MICE", paste0("n_lifestyle_cat", level), 
                       "S3_MICE", out_name, level)
  }
}

# ---- S4: Drop1st ----
cat("\n--------------------------------------------------------------------------------\n")
cat("S4: Drop First Follow-up\n")
cat("--------------------------------------------------------------------------------\n")

for (level in c("1", "2", "3+")) {
  run_meta_for_level("S4_Drop1st", paste0("n_lifestyle_cat", level), 
                     "S4_Drop1st", "Overall_drop1st", level)
}

# Combine summaries
ma_summary_df <- bind_rows(ma_summaries)
loo_summary_df <- bind_rows(loo_results)
egger_summary_df <- bind_rows(egger_results)

cat("\n================================================================================\n")
cat("META-ANALYSIS RESULTS SUMMARY\n")
cat("================================================================================\n")
cat("Total meta-analyses completed:", length(ma_results), "\n")
cat("Meta-analysis summary rows:", nrow(ma_summary_df), "\n")
cat("Leave-one-out analyses:", nrow(loo_summary_df), "\n")
cat("Egger's tests:", nrow(egger_summary_df), "\n")

if (length(ma_results) == 0) {
  cat("\n[WARNING] No meta-analyses were completed!\n")
  cat("Possible reasons:\n")
  cat("  - study_hr_data may be empty\n")
  cat("  - Filter conditions may not match any data\n")
  cat("  - All cohorts may have failed Cox models\n\n")
} else {
  cat("\nSuccessful meta-analyses:\n")
  for (key in names(ma_results)) {
    cat("  -", key, "\n")
  }
}
cat("================================================================================\n\n")

# =============================================================================
# PART 7: GENERATE PLOTS
# =============================================================================

cat("=============================================================================\n")
cat("PART 7: GENERATING PLOTS\n")
cat("=============================================================================\n\n")

close_all_devices()

# Check if we have any meta-analysis results
cat("Number of meta-analyses to plot:", length(ma_results), "\n\n")

# Forest plots
cat("--- Forest Plots ---\n")
forest_count <- 0

for (key in names(ma_results)) {
  ma <- ma_results[[key]]
  if (is.null(ma) || ma$k < 2) {
    cat("  [SKIP]", key, "- insufficient studies\n")
    next
  }
  
  cat("  [PLOT]", key, "(k =", ma$k, ")\n")
  
  filepath_pdf <- file.path(fig_dir, paste0("Forest_", key, ".pdf"))
  filepath_png <- file.path(fig_dir, paste0("Forest_", key, ".png"))
  
  # PDF
  tryCatch({
    pdf(file = filepath_pdf, width = 12, height = 7, useDingbats = FALSE)
    forest(ma, sortvar = ma$TE, prediction = TRUE, print.tau2 = TRUE, print.I2 = TRUE,
           common = TRUE, random = TRUE, leftcols = "studlab", leftlabs = "Study",
           rightcols = c("effect", "ci"), rightlabs = c("HR", "95% CI"),
           col.diamond = "steelblue", col.diamond.common = "darkblue", fontsize = 10)
    dev.off()
    cat("    [OK] PDF saved\n")
    forest_count <- forest_count + 1
  }, error = function(e) {
    cat("    [ERROR] PDF:", e$message, "\n")
    try(dev.off(), silent = TRUE)
  })
  
  # PNG
  tryCatch({
    png(filename = filepath_png, width = 12, height = 7, units = "in", res = 300)
    forest(ma, sortvar = ma$TE, prediction = TRUE, print.tau2 = TRUE, print.I2 = TRUE,
           common = TRUE, random = TRUE, leftcols = "studlab", leftlabs = "Study",
           rightcols = c("effect", "ci"), rightlabs = c("HR", "95% CI"),
           col.diamond = "steelblue", col.diamond.common = "darkblue", fontsize = 10)
    dev.off()
    cat("    [OK] PNG saved\n")
  }, error = function(e) {
    cat("    [ERROR] PNG:", e$message, "\n")
    try(dev.off(), silent = TRUE)
  })
}

cat("\nForest plots generated:", forest_count, "\n")

# Funnel plots
cat("\n--- Funnel Plots ---\n")
funnel_count <- 0

for (key in names(ma_results)) {
  ma <- ma_results[[key]]
  if (is.null(ma) || ma$k < 3) next
  
  cat("  [PLOT]", key, "\n")
  
  filepath_pdf <- file.path(funnel_dir, paste0("Funnel_", key, ".pdf"))
  filepath_png <- file.path(funnel_dir, paste0("Funnel_", key, ".png"))
  
  # PDF
  tryCatch({
    pdf(file = filepath_pdf, width = 8, height = 7, useDingbats = FALSE)
    funnel(ma, xlab = "Hazard Ratio (log)", studlab = TRUE, col = "steelblue")
    title(main = paste0("Funnel Plot: ", key))
    dev.off()
    cat("    [OK] PDF\n")
    funnel_count <- funnel_count + 1
  }, error = function(e) {
    cat("    [ERROR] PDF:", e$message, "\n")
    try(dev.off(), silent = TRUE)
  })
  
  # PNG
  tryCatch({
    png(filename = filepath_png, width = 8, height = 7, units = "in", res = 300)
    funnel(ma, xlab = "Hazard Ratio (log)", studlab = TRUE, col = "steelblue")
    title(main = paste0("Funnel Plot: ", key))
    dev.off()
    cat("    [OK] PNG\n")
  }, error = function(e) {
    cat("    [ERROR] PNG:", e$message, "\n")
    try(dev.off(), silent = TRUE)
  })
}

cat("\nFunnel plots generated:", funnel_count, "\n")
cat("Plots generation completed.\n\n")

# =============================================================================
# PART 8: SAVE RESULTS
# =============================================================================

cat("=============================================================================\n")
cat("PART 8: SAVING RESULTS\n")
cat("=============================================================================\n\n")

# Analysis description
ma_description <- data.frame(
  Analysis = c("Primary", "S1_5Level", "S2_HeavyDrink", "S3_MICE", "S4_Drop1st"),
  Description = c(
    "Primary: All outcomes ? 4-level (0/1/2/3+)",
    "Sensitivity 1: 5-level categories (0/1/2/3/4)",
    "Sensitivity 2: Heavy drinking definition",
    "Sensitivity 3: MICE imputed data",
    "Sensitivity 4: Drop first follow-up wave"
  )
)

# Save Excel
results_for_excel <- list(
  "Analysis_Description" = ma_description,
  "Study_Characteristics" = study_characteristics,
  "MA_Summary_All" = ma_summary_df,
  "Primary_Analysis" = ma_summary_df %>% filter(Analysis == "Primary"),
  "S1_5Level" = ma_summary_df %>% filter(Analysis == "S1_5Level"),
  "S2_HeavyDrink" = ma_summary_df %>% filter(Analysis == "S2_HeavyDrink"),
  "S3_MICE" = ma_summary_df %>% filter(Analysis == "S3_MICE"),
  "S4_Drop1st" = ma_summary_df %>% filter(Analysis == "S4_Drop1st"),
  "Leave_One_Out" = loo_summary_df,
  "Eggers_Test" = egger_summary_df,
  "Study_HR" = study_hr_data
)

save_to_excel(results_for_excel, "Meta_Analysis_Comprehensive.xlsx")
save_to_csv(ma_summary_df, "Meta_Analysis_Summary.csv")

# ===== Save Individual CSV Files for Each Component =====
cat("\n--- Saving Individual CSV Files ---\n")

# 1. Study Characteristics
if (nrow(study_characteristics) > 0) {
  save_to_csv(study_characteristics, "Meta_Study_Characteristics.csv")
  cat("  Saved: Meta_Study_Characteristics.csv (", nrow(study_characteristics), " rows)\n")
}

# 2. Study-Specific HR Results
if (nrow(study_hr_data) > 0) {
  save_to_csv(study_hr_data, "Meta_Study_Specific_HR.csv")
  cat("  Saved: Meta_Study_Specific_HR.csv (", nrow(study_hr_data), " rows)\n")
}

# 3. Meta-Analysis Summary by Analysis Type
if (nrow(ma_summary_df) > 0) {
  # Primary Analysis Only
  primary_ma <- ma_summary_df %>% filter(Analysis == "Primary")
  if (nrow(primary_ma) > 0) {
    save_to_csv(primary_ma, "Meta_Primary_Summary.csv")
    cat("  Saved: Meta_Primary_Summary.csv (", nrow(primary_ma), " rows)\n")
  }
  
  # Sensitivity Analyses
  for (sens_analysis in c("S1_5Level", "S2_HeavyDrink", "S3_MICE", "S4_Drop1st")) {
    sens_ma <- ma_summary_df %>% filter(Analysis == sens_analysis)
    if (nrow(sens_ma) > 0) {
      filename <- paste0("Meta_", sens_analysis, "_Summary.csv")
      save_to_csv(sens_ma, filename)
      cat("  Saved:", filename, "(", nrow(sens_ma), "rows)\n")
    }
  }
}

# 4. Leave-One-Out Results
if (nrow(loo_summary_df) > 0) {
  save_to_csv(loo_summary_df, "Meta_LeaveOneOut.csv")
  cat("  Saved: Meta_LeaveOneOut.csv (", nrow(loo_summary_df), " rows)\n")
}

# 5. Egger's Test Results
if (nrow(egger_summary_df) > 0) {
  save_to_csv(egger_summary_df, "Meta_EggersTest.csv")
  cat("  Saved: Meta_EggersTest.csv (", nrow(egger_summary_df), " rows)\n")
}

# ===== Create Publication-Ready Summary Table =====
if (nrow(ma_summary_df) > 0) {
  # Create HR (95% CI) format
  pub_summary <- ma_summary_df %>%
    mutate(
      Fixed_HR_CI = paste0(Fixed_HR, " (", Fixed_Lower, "-", Fixed_Upper, ")"),
      Random_HR_CI = paste0(Random_HR, " (", Random_Lower, "-", Random_Upper, ")"),
      I2_pct = paste0(I2, "%")
    ) %>%
    select(Analysis, Outcome, Exposure_Level, N_Studies, 
           Fixed_HR_CI, Fixed_P, Random_HR_CI, Random_P, I2_pct, Recommended)
  
  save_to_csv(pub_summary, "Meta_Summary_Publication.csv")
  cat("  Saved: Meta_Summary_Publication.csv\n")
}

# Save RDS
saveRDS(list(
  ma_results = ma_results,
  ma_summary = ma_summary_df,
  study_hr = study_hr_data,
  study_char = study_characteristics,
  loo = loo_summary_df,
  egger = egger_summary_df
), file.path(OUTPUT_DIR, "Meta_analysis_objects.rds"))

# =============================================================================
# PART 9: SUMMARY
# =============================================================================

cat("\n=============================================================================\n")
cat("ANALYSIS SUMMARY\n")
cat("=============================================================================\n\n")

cat("--- PRIMARY ANALYSIS: 3+ vs 0 ---\n")
primary_3plus <- ma_summary_df %>% filter(Analysis == "Primary" & Exposure_Level == "3+")
if (nrow(primary_3plus) > 0) {
  for (i in seq_len(nrow(primary_3plus))) {
    row <- primary_3plus[i,]
    cat("  ", row$Outcome, "(", row$Outcome_Type, "):\n")
    cat("    Fixed: HR =", row$Fixed_HR, "(", row$Fixed_Lower, "-", row$Fixed_Upper, ")\n")
    cat("    Random: HR =", row$Random_HR, "(", row$Random_Lower, "-", row$Random_Upper, ")\n")
    cat("    I² =", row$I2, "%, Recommended:", row$Recommended, "\n\n")
  }
}

cat("--- SENSITIVITY ANALYSES ---\n")
for (analysis in c("S1_5Level", "S2_HeavyDrink", "S3_MICE", "S4_Drop1st")) {
  n_results <- sum(ma_summary_df$Analysis == analysis)
  cat("  ", analysis, ":", n_results, "results\n")
}

# =============================================================================
# COMPLETE
# =============================================================================

cat("\n")
cat("###############################################################################\n")
cat("#                META-ANALYSIS COMPLETED SUCCESSFULLY                        #\n")
cat("###############################################################################\n")
cat("\nTotal meta-analyses:", length(ma_results), "\n")
cat("Results saved to:", OUTPUT_DIR, "\n")
cat("  - Meta_Analysis_Comprehensive.xlsx\n")
cat("  - Meta_Analysis_Summary.csv\n")
cat("  - Forest plots:", fig_dir, "\n")
cat("  - Funnel plots:", funnel_dir, "\n")

close_all_devices()
