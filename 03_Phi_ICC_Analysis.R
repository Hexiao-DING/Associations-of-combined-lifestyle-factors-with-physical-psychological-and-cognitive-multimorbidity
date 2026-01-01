###############################################################################
# 03_Phi_ICC_Analysis.R
# ============================================================================
# Project: Lifestyle Factors and PPC-MM Association Study
# Purpose: Phi Coefficient and ICC Analysis with P-values and 95% CI
# ============================================================================
# This script includes:
#   1. Phi coefficient calculation with:
#      - Bootstrap 95% confidence intervals (parallel)
#      - Chi-square test P-values
#      - Bonferroni correction for multiple comparisons
#   2. ICC calculation using glmer + performance::icc with:
#      - Delta method 95% confidence intervals
#      - P-values from likelihood ratio test (LRT)
#      - Bonferroni correction
#   3. Publication-ready output tables
###############################################################################

# =============================================================================
# PART 1: Environment Setup
# =============================================================================

source("C:/Users/user/Desktop/Scientific Research/UKB Project/Jung Sun Lab/Project1_Lifestyle_PPM_Hexiao,Hongtao/Code/00_Functions_and_Setup.R")

# Install parallel computing packages for Phi bootstrap
if (!requireNamespace("parallel", quietly = TRUE)) {
  install.packages("parallel", dependencies = TRUE)
}
if (!requireNamespace("foreach", quietly = TRUE)) {
  install.packages("foreach", dependencies = TRUE)
}
if (!requireNamespace("doParallel", quietly = TRUE)) {
  install.packages("doParallel", dependencies = TRUE)
}

library(parallel)
library(foreach)
library(doParallel)

cat("\n")
cat("================================================================\n")
cat("     Phi Coefficient and ICC Analysis (with P-values & CI)      \n")
cat("         *** PARALLEL COMPUTING FOR PHI BOOTSTRAP ***           \n")
cat("================================================================\n\n")

# =============================================================================
# PART 1B: Setup Parallel Computing (for Phi Bootstrap)
# =============================================================================

# Detect number of CPU cores
n_cores <- detectCores()
# Use all cores minus 1 (leave one for system)
n_workers <- max(1, n_cores - 1)

cat("Parallel Computing Setup:\n")
cat("  - Total CPU cores detected:", n_cores, "\n")
cat("  - Workers to use:", n_workers, "\n")

# Register parallel backend
cl <- makeCluster(n_workers)
registerDoParallel(cl)

cat("  - Parallel cluster registered successfully!\n\n")

# =============================================================================
# PART 2: Load Pooled Data
# =============================================================================

cat("--- Loading pooled data ---\n")

pooled_data <- readRDS(file.path(OUTPUT_DIR, "Pooled_main_data.rds"))

cat("Pooled data loaded, sample size:", nrow(pooled_data), "\n")

# =============================================================================
# PART 3: Helper Functions for Phi with P-value and CI
# =============================================================================

#' Calculate Phi coefficient with Chi-square test P-value and PARALLEL Bootstrap CI
#' @param x Binary variable 1
#' @param y Binary variable 2
#' @param n_boot Number of bootstrap samples for CI
#' @param conf_level Confidence level (default 0.95)
#' @param use_parallel Use parallel computing (default TRUE)
#' @return List with phi, chi2, p_value, ci_lower, ci_upper
calculate_phi_complete <- function(x, y, n_boot = 1000, conf_level = 0.95, use_parallel = TRUE) {
  # Remove missing values
  complete_idx <- complete.cases(x, y)
  x <- x[complete_idx]
  y <- y[complete_idx]
  n <- length(x)
  
  if (n < 10) {
    return(list(
      phi = NA_real_,
      chi2 = NA_real_,
      p_value = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      n = n
    ))
  }
  
  # Calculate Phi coefficient
  phi_val <- tryCatch({
    cor(x, y, method = "pearson")
  }, error = function(e) NA_real_)
  
  # Calculate Chi-square statistic and P-value
  chi2_result <- tryCatch({
    tbl <- table(factor(x, levels = c(0, 1)), factor(y, levels = c(0, 1)))
    test <- chisq.test(tbl, correct = FALSE)
    list(chi2 = test$statistic, p_value = test$p.value)
  }, error = function(e) {
    list(chi2 = NA_real_, p_value = NA_real_)
  })
  
  # PARALLEL Bootstrap confidence interval
  boot_phis <- tryCatch({
    if (use_parallel && foreach::getDoParRegistered()) {
      # Parallel bootstrap using foreach
      foreach(i = 1:n_boot, .combine = c, .packages = c("stats")) %dopar% {
        idx <- sample(n, n, replace = TRUE)
        cor(x[idx], y[idx], method = "pearson")
      }
    } else {
      # Sequential fallback
      replicate(n_boot, {
        idx <- sample(n, n, replace = TRUE)
        cor(x[idx], y[idx], method = "pearson")
      })
    }
  }, error = function(e) {
    # Fallback to sequential if parallel fails
    replicate(n_boot, {
      idx <- sample(n, n, replace = TRUE)
      cor(x[idx], y[idx], method = "pearson")
    })
  })
  
  # Calculate CI from bootstrap distribution
  alpha <- 1 - conf_level
  ci <- quantile(boot_phis, c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  
  list(
    phi = phi_val,
    chi2 = as.numeric(chi2_result$chi2),
    p_value = chi2_result$p_value,
    ci_lower = ci[1],
    ci_upper = ci[2],
    n = n
  )
}

#' Format P-value with 4 decimal places
#' @param p P-value
#' @return Formatted string
format_pvalue <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.0001) return("<0.0001")
  sprintf("%.4f", p)
}

#' Apply Bonferroni correction
#' @param p_values Vector of P-values
#' @return Vector of corrected P-values
bonferroni_correct <- function(p_values) {
  n_tests <- sum(!is.na(p_values))
  pmin(p_values * n_tests, 1)
}

# =============================================================================
# PART 4: Phi Coefficient Analysis
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("         Phi Coefficient Analysis (Bootstrap CI + P-values)     \n")
cat("================================================================\n")

# Lifestyle variables
lifestyle_vars <- c("unhealthy_drink", "unhealthy_smoke", 
                    "unhealthy_pa", "unhealthy_soc")

lifestyle_labels <- c("Drinking", "Smoking", "Physical Inactivity", "Social Isolation")
names(lifestyle_labels) <- lifestyle_vars

# Number of pairwise comparisons for Bonferroni
n_pairs <- choose(length(lifestyle_vars), 2)  # 6 pairs
cat("Number of pairwise comparisons:", n_pairs, "\n")
cat("Bonferroni-adjusted significance threshold:", 0.05 / n_pairs, "\n\n")

# -----------------------------------------------------------------------------
# 4.1 Overall Phi Coefficient (Pooled Data)
# -----------------------------------------------------------------------------

cat("=== Calculating Phi Coefficients (Pooled Data) ===\n")
cat("Running bootstrap (1000 iterations per pair)...\n\n")

phi_results_overall <- data.frame(
  Variable_1 = character(),
  Variable_2 = character(),
  N = integer(),
  Phi = numeric(),
  CI_95_Lower = numeric(),
  CI_95_Upper = numeric(),
  Chi_Square = numeric(),
  P_value = numeric(),
  P_value_Bonf = numeric(),
  Significance = character(),
  stringsAsFactors = FALSE
)

for (i in 1:(length(lifestyle_vars) - 1)) {
  for (j in (i + 1):length(lifestyle_vars)) {
    cat("  Processing:", lifestyle_labels[i], "vs", lifestyle_labels[j], "...")
    
    result <- calculate_phi_complete(
      pooled_data[[lifestyle_vars[i]]],
      pooled_data[[lifestyle_vars[j]]],
      n_boot = 1000
    )
    
    phi_results_overall <- rbind(phi_results_overall, data.frame(
      Variable_1 = lifestyle_labels[i],
      Variable_2 = lifestyle_labels[j],
      N = result$n,
      Phi = round(result$phi, 4),
      CI_95_Lower = round(result$ci_lower, 4),
      CI_95_Upper = round(result$ci_upper, 4),
      Chi_Square = round(result$chi2, 2),
      P_value = result$p_value,
      P_value_Bonf = NA_real_,  # Will calculate after
      Significance = "",
      stringsAsFactors = FALSE
    ))
    
    cat(" Phi =", round(result$phi, 3), "\n")
  }
}

# Apply Bonferroni correction
phi_results_overall$P_value_Bonf <- bonferroni_correct(phi_results_overall$P_value)

# Add significance indicators
phi_results_overall$Significance <- sapply(phi_results_overall$P_value_Bonf, function(p) {
  if (is.na(p)) return("")
  if (p < 0.0001) return("****")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
})

# Format P-values for display
phi_results_overall$P_value_formatted <- sapply(phi_results_overall$P_value, format_pvalue)
phi_results_overall$P_value_Bonf_formatted <- sapply(phi_results_overall$P_value_Bonf, format_pvalue)

# Print results
cat("\n=== Phi Coefficient Results (Pooled Data) ===\n")
print(phi_results_overall[, c("Variable_1", "Variable_2", "N", "Phi", 
                               "CI_95_Lower", "CI_95_Upper", "Chi_Square",
                               "P_value_formatted", "P_value_Bonf_formatted", "Significance")])

# Create publication-ready table
phi_table_publication <- phi_results_overall %>%
  mutate(
    `95% CI` = paste0("(", sprintf("%.3f", CI_95_Lower), ", ", sprintf("%.3f", CI_95_Upper), ")"),
    `Phi (95% CI)` = paste0(sprintf("%.3f", Phi), " ", `95% CI`),
    `χ²` = sprintf("%.2f", Chi_Square),
    `P-value` = P_value_formatted,
    `P-value (Bonferroni)` = P_value_Bonf_formatted
  ) %>%
  select(Variable_1, Variable_2, N, `Phi (95% CI)`, `χ²`, `P-value`, `P-value (Bonferroni)`, Significance)

# -----------------------------------------------------------------------------
# 4.2 Phi Coefficients by Cohort
# -----------------------------------------------------------------------------

cat("\n\n=== Phi Coefficients by Cohort ===\n")

cohort_names <- unique(pooled_data$cohort)
phi_by_cohort_list <- list()

for (coh in cohort_names) {
  cat("\n--- Cohort:", coh, "---\n")
  
  df_cohort <- pooled_data %>% filter(cohort == coh)
  
  if (nrow(df_cohort) < 50) {
    cat("Sample size too small (N =", nrow(df_cohort), "), skipping\n")
    next
  }
  
  cat("Sample size: N =", nrow(df_cohort), "\n")
  
  phi_cohort_results <- data.frame(
    Cohort = character(),
    Variable_1 = character(),
    Variable_2 = character(),
    N = integer(),
    Phi = numeric(),
    CI_95_Lower = numeric(),
    CI_95_Upper = numeric(),
    Chi_Square = numeric(),
    P_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:(length(lifestyle_vars) - 1)) {
    for (j in (i + 1):length(lifestyle_vars)) {
      result <- calculate_phi_complete(
        df_cohort[[lifestyle_vars[i]]],
        df_cohort[[lifestyle_vars[j]]],
        n_boot = 500  # Fewer iterations for cohort-level
      )
      
      phi_cohort_results <- rbind(phi_cohort_results, data.frame(
        Cohort = coh,
        Variable_1 = lifestyle_labels[i],
        Variable_2 = lifestyle_labels[j],
        N = result$n,
        Phi = round(result$phi, 4),
        CI_95_Lower = round(result$ci_lower, 4),
        CI_95_Upper = round(result$ci_upper, 4),
        Chi_Square = round(result$chi2, 2),
        P_value = result$p_value,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Apply Bonferroni correction within cohort
  phi_cohort_results$P_value_Bonf <- bonferroni_correct(phi_cohort_results$P_value)
  phi_cohort_results$Significance <- sapply(phi_cohort_results$P_value_Bonf, function(p) {
    if (is.na(p)) return("")
    if (p < 0.0001) return("****")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    return("ns")
  })
  
  phi_by_cohort_list[[coh]] <- phi_cohort_results
  
  # Print summary
  cat("\nResults for", coh, ":\n")
  for (k in 1:nrow(phi_cohort_results)) {
    row <- phi_cohort_results[k, ]
    cat(sprintf("  %s vs %s: Phi = %.3f (%.3f, %.3f), P = %s %s\n",
                row$Variable_1, row$Variable_2,
                row$Phi, row$CI_95_Lower, row$CI_95_Upper,
                format_pvalue(row$P_value_Bonf), row$Significance))
  }
}

# Combine all cohort results
phi_all_cohorts <- bind_rows(phi_by_cohort_list)

# -----------------------------------------------------------------------------
# 4.3 Phi Coefficient Matrix (for heatmap)
# -----------------------------------------------------------------------------

# Create matrix format
phi_matrix_overall <- matrix(NA, 
                             nrow = length(lifestyle_vars), 
                             ncol = length(lifestyle_vars),
                             dimnames = list(lifestyle_labels, lifestyle_labels))

p_matrix_overall <- matrix(NA, 
                           nrow = length(lifestyle_vars), 
                           ncol = length(lifestyle_vars),
                           dimnames = list(lifestyle_labels, lifestyle_labels))

for (i in 1:length(lifestyle_vars)) {
  for (j in 1:length(lifestyle_vars)) {
    if (i == j) {
      phi_matrix_overall[i, j] <- 1.000
      p_matrix_overall[i, j] <- NA
    } else if (i < j) {
      idx <- which(phi_results_overall$Variable_1 == lifestyle_labels[i] & 
                     phi_results_overall$Variable_2 == lifestyle_labels[j])
      phi_matrix_overall[i, j] <- phi_results_overall$Phi[idx]
      phi_matrix_overall[j, i] <- phi_results_overall$Phi[idx]
      p_matrix_overall[i, j] <- phi_results_overall$P_value_Bonf[idx]
      p_matrix_overall[j, i] <- phi_results_overall$P_value_Bonf[idx]
    }
  }
}

cat("\n=== Phi Coefficient Matrix ===\n")
print(round(phi_matrix_overall, 3))

cat("\n=== P-value Matrix (Bonferroni-corrected) ===\n")
print(round(p_matrix_overall, 4))

# =============================================================================
# PART 5: ICC Analysis (using performance::icc with LRT P-value)
# =============================================================================

cat("\n\n")
cat("================================================================\n")
cat("         ICC Analysis (glmer + performance::icc)                \n")
cat("================================================================\n")

#' Calculate ICC with CI and P-value using glmer
#' @param var_name Variable name
#' @param cluster_var Clustering variable name
#' @param data Data frame
#' @return List with icc, ci_lower, ci_upper, p_value
calculate_icc_complete <- function(var_name, cluster_var, data) {
  
  # Prepare data
  df <- data[, c(var_name, cluster_var)]
  df <- df[complete.cases(df), ]
  df[[cluster_var]] <- as.factor(df[[cluster_var]])
  
  n_obs <- nrow(df)
  n_clusters <- length(unique(df[[cluster_var]]))
  
  if (n_obs < 50 || n_clusters < 2) {
    return(list(icc = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, 
                p_value = NA_real_, n_obs = n_obs, n_clusters = n_clusters))
  }
  
  # Fit random intercept model
  formula_full <- as.formula(paste0(var_name, " ~ 1 + (1|", cluster_var, ")"))
  formula_null <- as.formula(paste0(var_name, " ~ 1"))
  
  result <- tryCatch({
    # Fit models
    fit_full <- glmer(formula_full, data = df, family = binomial,
                      control = glmerControl(optimizer = "bobyqa", 
                                             optCtrl = list(maxfun = 100000)))
    fit_null <- glm(formula_null, data = df, family = binomial)
    
    # Extract ICC using performance package
    icc_obj <- performance::icc(fit_full)
    icc_val <- icc_obj$ICC_adjusted
    
    # Calculate CI using confint (profile likelihood)
    ci <- tryCatch({
      vc <- as.data.frame(VarCorr(fit_full))
      var_cluster <- vc$vcov[1]
      var_resid <- pi^2 / 3  # logistic distribution variance
      
      # Approximate CI using delta method
      se_icc <- sqrt((2 * (1 - icc_val)^2 * icc_val^2) / (n_clusters - 1))
      ci_lower <- max(0, icc_val - 1.96 * se_icc)
      ci_upper <- min(1, icc_val + 1.96 * se_icc)
      
      c(ci_lower, ci_upper)
    }, error = function(e) c(NA_real_, NA_real_))
    
    # LRT P-value
    lrt <- anova(fit_full, fit_null, test = "Chisq")
    p_val <- lrt$`Pr(>Chisq)`[2]
    
    list(
      icc = icc_val,
      ci_lower = ci[1],
      ci_upper = ci[2],
      p_value = p_val,
      n_obs = n_obs,
      n_clusters = n_clusters
    )
    
  }, error = function(e) {
    cat(" [glmer error:", e$message, "]")
    list(icc = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, 
         p_value = NA_real_, n_obs = n_obs, n_clusters = n_clusters)
  })
  
  return(result)
}

# -----------------------------------------------------------------------------
# 5.1 ICC for Lifestyle Variables
# -----------------------------------------------------------------------------

cat("\n=== ICC for Lifestyle Variables ===\n")
cat("Using glmer + performance::icc with LRT P-value...\n\n")

icc_lifestyle_results <- data.frame(
  Variable = character(),
  ICC = numeric(),
  CI_95_Lower = numeric(),
  CI_95_Upper = numeric(),
  P_value = numeric(),
  P_value_Bonf = numeric(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

n_lifestyle_tests <- length(lifestyle_vars)

for (lv in lifestyle_vars) {
  cat("  Processing:", lifestyle_labels[lv], "...")
  
  icc_result <- calculate_icc_complete(lv, "cohort", pooled_data)
  
  # Interpretation
  icc_val <- icc_result$icc
  interp <- dplyr::case_when(
    is.na(icc_val) ~ "Cannot calculate",
    icc_val < 0.05 ~ "Very low clustering",
    icc_val < 0.10 ~ "Low clustering",
    icc_val < 0.20 ~ "Moderate clustering",
    icc_val < 0.30 ~ "High clustering",
    TRUE ~ "Very high clustering"
  )
  
  icc_lifestyle_results <- rbind(icc_lifestyle_results, data.frame(
    Variable = lifestyle_labels[lv],
    ICC = round(icc_result$icc, 4),
    CI_95_Lower = round(icc_result$ci_lower, 4),
    CI_95_Upper = round(icc_result$ci_upper, 4),
    P_value = icc_result$p_value,
    P_value_Bonf = NA_real_,
    Interpretation = interp,
    stringsAsFactors = FALSE
  ))
  
  if (!is.na(icc_result$icc)) {
    cat(" ICC =", sprintf("%.4f", icc_result$icc), 
        "(", sprintf("%.4f", icc_result$ci_lower), "-", sprintf("%.4f", icc_result$ci_upper), ")",
        "P =", format_pvalue(icc_result$p_value), "\n")
  } else {
    cat(" ICC = NA\n")
  }
}

# Apply Bonferroni correction
icc_lifestyle_results$P_value_Bonf <- bonferroni_correct(icc_lifestyle_results$P_value)

# Format P-values
icc_lifestyle_results$P_value_formatted <- sapply(icc_lifestyle_results$P_value, format_pvalue)
icc_lifestyle_results$P_value_Bonf_formatted <- sapply(icc_lifestyle_results$P_value_Bonf, format_pvalue)

cat("\n=== ICC Results for Lifestyle Variables ===\n")
print(icc_lifestyle_results[, c("Variable", "ICC", "CI_95_Lower", "CI_95_Upper", 
                                 "P_value_formatted", "P_value_Bonf_formatted", "Interpretation")])

# -----------------------------------------------------------------------------
# 5.2 ICC for Outcome Variables
# -----------------------------------------------------------------------------

cat("\n\n=== ICC for Outcome Variables (All 5 PPC-MM Outcomes) ===\n")

outcome_vars <- c("event_ppcmm", "event_mm_phys_psych", "event_mm_phys_cog",
                  "event_mm_psych_cog", "event_mm_all_three")
outcome_labels <- c("Overall PPC-MM", "P1P2 (Physical-Psychological)",
                    "P1C (Physical-Cognitive)", "P2C (Psychological-Cognitive)",
                    "P1P2C (All Three)")
names(outcome_labels) <- outcome_vars

icc_outcome_results <- data.frame(
  Variable = character(),
  ICC = numeric(),
  CI_95_Lower = numeric(),
  CI_95_Upper = numeric(),
  P_value = numeric(),
  P_value_Bonf = numeric(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

n_outcome_tests <- length(outcome_vars)

for (ov in outcome_vars) {
  cat("  Processing:", outcome_labels[ov], "...")
  
  # Check if variable exists
  if (!ov %in% names(pooled_data)) {
    cat(" Variable not found, skipping\n")
    next
  }
  
  icc_result <- calculate_icc_complete(ov, "cohort", pooled_data)
  
  icc_val <- icc_result$icc
  interp <- dplyr::case_when(
    is.na(icc_val) ~ "Cannot calculate",
    icc_val < 0.05 ~ "Very low clustering",
    icc_val < 0.10 ~ "Low clustering",
    icc_val < 0.20 ~ "Moderate clustering",
    icc_val < 0.30 ~ "High clustering",
    TRUE ~ "Very high clustering"
  )
  
  icc_outcome_results <- rbind(icc_outcome_results, data.frame(
    Variable = outcome_labels[ov],
    ICC = round(icc_result$icc, 4),
    CI_95_Lower = round(icc_result$ci_lower, 4),
    CI_95_Upper = round(icc_result$ci_upper, 4),
    P_value = icc_result$p_value,
    P_value_Bonf = NA_real_,
    Interpretation = interp,
    stringsAsFactors = FALSE
  ))
  
  if (!is.na(icc_result$icc)) {
    cat(" ICC =", sprintf("%.4f", icc_result$icc), 
        "(", sprintf("%.4f", icc_result$ci_lower), "-", sprintf("%.4f", icc_result$ci_upper), ")",
        "P =", format_pvalue(icc_result$p_value), "\n")
  } else {
    cat(" ICC = NA\n")
  }
}

# Apply Bonferroni correction
icc_outcome_results$P_value_Bonf <- bonferroni_correct(icc_outcome_results$P_value)

# Format P-values
icc_outcome_results$P_value_formatted <- sapply(icc_outcome_results$P_value, format_pvalue)
icc_outcome_results$P_value_Bonf_formatted <- sapply(icc_outcome_results$P_value_Bonf, format_pvalue)

cat("\n=== ICC Results for Outcome Variables ===\n")
print(icc_outcome_results[, c("Variable", "ICC", "CI_95_Lower", "CI_95_Upper", 
                               "P_value_formatted", "P_value_Bonf_formatted", "Interpretation")])

# -----------------------------------------------------------------------------
# 5.3 ICC for Baseline Health Status
# -----------------------------------------------------------------------------

cat("\n\n=== ICC for Baseline Health Status ===\n")

baseline_vars <- c("physical_base", "psych_base", "cog_base")
baseline_labels <- c("Physical disease", "Depression", "Cognitive impairment")
names(baseline_labels) <- baseline_vars

icc_baseline_results <- data.frame(
  Variable = character(),
  ICC = numeric(),
  CI_95_Lower = numeric(),
  CI_95_Upper = numeric(),
  P_value = numeric(),
  P_value_Bonf = numeric(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

for (bv in baseline_vars) {
  cat("  Processing:", baseline_labels[bv], "...")
  
  if (!bv %in% names(pooled_data)) {
    cat(" Variable not found, skipping\n")
    next
  }
  
  icc_result <- calculate_icc_complete(bv, "cohort", pooled_data)
  
  icc_val <- icc_result$icc
  interp <- dplyr::case_when(
    is.na(icc_val) ~ "Cannot calculate",
    icc_val < 0.05 ~ "Very low clustering",
    icc_val < 0.10 ~ "Low clustering",
    icc_val < 0.20 ~ "Moderate clustering",
    icc_val < 0.30 ~ "High clustering",
    TRUE ~ "Very high clustering"
  )
  
  icc_baseline_results <- rbind(icc_baseline_results, data.frame(
    Variable = baseline_labels[bv],
    ICC = round(icc_result$icc, 4),
    CI_95_Lower = round(icc_result$ci_lower, 4),
    CI_95_Upper = round(icc_result$ci_upper, 4),
    P_value = icc_result$p_value,
    P_value_Bonf = NA_real_,
    Interpretation = interp,
    stringsAsFactors = FALSE
  ))
  
  if (!is.na(icc_result$icc)) {
    cat(" ICC =", sprintf("%.4f", icc_result$icc), 
        "(", sprintf("%.4f", icc_result$ci_lower), "-", sprintf("%.4f", icc_result$ci_upper), ")",
        "P =", format_pvalue(icc_result$p_value), "\n")
  } else {
    cat(" ICC = NA\n")
  }
}

# Apply Bonferroni correction
icc_baseline_results$P_value_Bonf <- bonferroni_correct(icc_baseline_results$P_value)
icc_baseline_results$P_value_formatted <- sapply(icc_baseline_results$P_value, format_pvalue)
icc_baseline_results$P_value_Bonf_formatted <- sapply(icc_baseline_results$P_value_Bonf, format_pvalue)

cat("\n=== ICC Results for Baseline Health Status ===\n")
print(icc_baseline_results[, c("Variable", "ICC", "CI_95_Lower", "CI_95_Upper", 
                                "P_value_formatted", "P_value_Bonf_formatted", "Interpretation")])

# =============================================================================
# PART 5B: SUBGROUP ANALYSIS - Phi Coefficients by Age and Sex
# =============================================================================

cat("\n\n")
cat("================================================================\n")
cat("         SUBGROUP: Phi Coefficients by Age and Sex              \n")
cat("================================================================\n")

# DEBUG: Check data before subgroup analysis
cat("\n[DEBUG] Before Phi subgroup analysis:\n")
cat("  Total N in pooled_data:", nrow(pooled_data), "\n")
cat("  age_baseline range:", min(pooled_data$age_baseline, na.rm=TRUE), "-", 
    max(pooled_data$age_baseline, na.rm=TRUE), "\n")

# Create age groups if not exists
if (!"age_group_4" %in% names(pooled_data)) {
  cat("  Creating age_group_4 variable...\n")
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
}

# Print distributions
cat("\n[DEBUG] Age group distribution:\n")
print(table(pooled_data$age_group_4, useNA = "ifany"))
cat("\n[DEBUG] Sex distribution:\n")
print(table(pooled_data$sex, useNA = "ifany"))

lifestyle_vars <- c("unhealthy_drink", "unhealthy_smoke", "unhealthy_pa", "unhealthy_soc")
lifestyle_labels <- c("Drinking", "Smoking", "Physical Inactivity", "Social Isolation")

# --- Phi by Age Group ---
cat("\n--- Phi Coefficients by Age Group ---\n")

phi_age_stratified <- list()
age_groups <- c("50-59", "60-69", "70-79", "80+")

for (ag in age_groups) {
  cat("\n  Age Group:", ag, "\n")
  
  df_age <- pooled_data %>% filter(age_group_4 == ag)
  n_ag <- nrow(df_age)
  cat("    N =", n_ag, "\n")
  
  if (n_ag < 100) {
    cat("    [SKIP] Insufficient sample size\n")
    next
  }
  
  for (i in 1:(length(lifestyle_vars) - 1)) {
    for (j in (i + 1):length(lifestyle_vars)) {
      var1 <- lifestyle_vars[i]
      var2 <- lifestyle_vars[j]
      
      if (!var1 %in% names(df_age) || !var2 %in% names(df_age)) next
      
      x <- as.numeric(df_age[[var1]])
      y <- as.numeric(df_age[[var2]])
      
      complete_idx <- complete.cases(x, y)
      x <- x[complete_idx]
      y <- y[complete_idx]
      
      if (length(x) < 50) next
      
      # Calculate Phi (simplified - no bootstrap for speed)
      tryCatch({
        ct <- table(x, y)
        chi_result <- chisq.test(ct, correct = FALSE)
        phi <- sqrt(chi_result$statistic / length(x))
        
        phi_age_stratified[[paste0(ag, "_", var1, "_", var2)]] <- data.frame(
          Age_Group = ag,
          Variable_1 = var1,
          Variable_2 = var2,
          Variable_1_Label = lifestyle_labels[i],
          Variable_2_Label = lifestyle_labels[j],
          N = length(x),
          Phi = round(as.numeric(phi), 4),
          Chi_Square = round(chi_result$statistic, 2),
          P_value = chi_result$p.value,
          P_value_formatted = format.pval(chi_result$p.value, digits = 4),
          stringsAsFactors = FALSE
        )
      }, error = function(e) NULL)
    }
  }
}

# Combine age-stratified Phi
if (length(phi_age_stratified) > 0) {
  phi_age_df <- bind_rows(phi_age_stratified)
  phi_age_df$P_value_Bonf <- bonferroni_correct(phi_age_df$P_value)
  phi_age_df$P_value_Bonf_formatted <- sapply(phi_age_df$P_value_Bonf, format_pvalue)
} else {
  phi_age_df <- data.frame()
}

# Save ALWAYS
cat("\n=== Age-Stratified Phi Summary ===\n")
cat("Total rows:", nrow(phi_age_df), "\n")

tryCatch({
  phi_age_file <- file.path(OUTPUT_DIR, "Phi_Age_Stratified.csv")
  if (nrow(phi_age_df) > 0) {
    print(phi_age_df %>% select(Age_Group, Variable_1_Label, Variable_2_Label, N, Phi, P_value_formatted))
    write.csv(phi_age_df, phi_age_file, row.names = FALSE)
  } else {
    placeholder <- data.frame(Note = "No age-stratified Phi results", Age_Groups = "50-59, 60-69, 70-79, 80+")
    write.csv(placeholder, phi_age_file, row.names = FALSE)
  }
  cat("\n  [OK] Saved: Phi_Age_Stratified.csv\n")
}, error = function(e) {
  cat("  [ERROR]:", e$message, "\n")
})

# --- Phi by Sex ---
cat("\n--- Phi Coefficients by Sex ---\n")

phi_sex_stratified <- list()
sex_groups <- c("Men", "Women")

for (sg in sex_groups) {
  cat("\n  Sex:", sg, "\n")
  
  df_sex <- pooled_data %>% filter(sex == sg)
  n_sg <- nrow(df_sex)
  cat("    N =", n_sg, "\n")
  
  if (n_sg < 100) {
    cat("    [SKIP] Insufficient sample size\n")
    next
  }
  
  for (i in 1:(length(lifestyle_vars) - 1)) {
    for (j in (i + 1):length(lifestyle_vars)) {
      var1 <- lifestyle_vars[i]
      var2 <- lifestyle_vars[j]
      
      if (!var1 %in% names(df_sex) || !var2 %in% names(df_sex)) next
      
      x <- as.numeric(df_sex[[var1]])
      y <- as.numeric(df_sex[[var2]])
      
      complete_idx <- complete.cases(x, y)
      x <- x[complete_idx]
      y <- y[complete_idx]
      
      if (length(x) < 50) next
      
      tryCatch({
        ct <- table(x, y)
        chi_result <- chisq.test(ct, correct = FALSE)
        phi <- sqrt(chi_result$statistic / length(x))
        
        phi_sex_stratified[[paste0(sg, "_", var1, "_", var2)]] <- data.frame(
          Sex = sg,
          Variable_1 = var1,
          Variable_2 = var2,
          Variable_1_Label = lifestyle_labels[i],
          Variable_2_Label = lifestyle_labels[j],
          N = length(x),
          Phi = round(as.numeric(phi), 4),
          Chi_Square = round(chi_result$statistic, 2),
          P_value = chi_result$p.value,
          P_value_formatted = format.pval(chi_result$p.value, digits = 4),
          stringsAsFactors = FALSE
        )
      }, error = function(e) NULL)
    }
  }
}

# Combine sex-stratified Phi
if (length(phi_sex_stratified) > 0) {
  phi_sex_df <- bind_rows(phi_sex_stratified)
  phi_sex_df$P_value_Bonf <- bonferroni_correct(phi_sex_df$P_value)
  phi_sex_df$P_value_Bonf_formatted <- sapply(phi_sex_df$P_value_Bonf, format_pvalue)
} else {
  phi_sex_df <- data.frame()
}

# Save ALWAYS
cat("\n=== Sex-Stratified Phi Summary ===\n")
cat("Total rows:", nrow(phi_sex_df), "\n")

tryCatch({
  phi_sex_file <- file.path(OUTPUT_DIR, "Phi_Sex_Stratified.csv")
  if (nrow(phi_sex_df) > 0) {
    print(phi_sex_df %>% select(Sex, Variable_1_Label, Variable_2_Label, N, Phi, P_value_formatted))
    write.csv(phi_sex_df, phi_sex_file, row.names = FALSE)
  } else {
    placeholder <- data.frame(Note = "No sex-stratified Phi results", Sex_Groups = "Men, Women")
    write.csv(placeholder, phi_sex_file, row.names = FALSE)
  }
  cat("\n  [OK] Saved: Phi_Sex_Stratified.csv\n")
}, error = function(e) {
  cat("  [ERROR]:", e$message, "\n")
})

# =============================================================================
# PART 5C: SUBGROUP ANALYSIS - ICC by Age and Sex
# =============================================================================

cat("\n\n")
cat("================================================================\n")
cat("         SUBGROUP: ICC by Age and Sex                           \n")
cat("================================================================\n")

# --- ICC by Age Group ---
cat("\n--- ICC by Age Group ---\n")

icc_age_stratified <- list()

for (ag in age_groups) {
  cat("\n  Age Group:", ag, "\n")
  
  df_age <- pooled_data %>% filter(age_group_4 == ag)
  n_ag <- nrow(df_age)
  
  if (n_ag < 100) {
    cat("    [SKIP] Insufficient sample size\n")
    next
  }
  
  for (i in seq_along(lifestyle_vars)) {
    lv <- lifestyle_vars[i]
    if (!lv %in% names(df_age)) next
    
    tryCatch({
      formula_str <- paste0(lv, " ~ 1 + (1|cohort)")
      model <- lme4::glmer(as.formula(formula_str), data = df_age, family = binomial,
                           control = lme4::glmerControl(optimizer = "bobyqa"))
      
      icc_result <- performance::icc(model)
      icc_val <- icc_result$ICC_adjusted
      
      icc_age_stratified[[paste0(ag, "_", lv)]] <- data.frame(
        Age_Group = ag,
        Variable = lv,
        Variable_Label = lifestyle_labels[i],
        N = n_ag,
        N_Cohorts = length(unique(df_age$cohort)),
        ICC = round(icc_val, 4),
        stringsAsFactors = FALSE
      )
      
      cat("    ", lifestyle_labels[i], ": ICC =", round(icc_val, 4), "\n")
    }, error = function(e) {
      cat("    ", lifestyle_labels[i], ": [ERROR]\n")
    })
  }
}

# Combine age-stratified ICC
if (length(icc_age_stratified) > 0) {
  icc_age_df <- bind_rows(icc_age_stratified)
} else {
  icc_age_df <- data.frame()
}

# Save ALWAYS
cat("\n=== Age-Stratified ICC Summary ===\n")
cat("Total rows:", nrow(icc_age_df), "\n")

tryCatch({
  icc_age_file <- file.path(OUTPUT_DIR, "ICC_Age_Stratified.csv")
  if (nrow(icc_age_df) > 0) {
    print(icc_age_df %>% select(Age_Group, Variable_Label, N, ICC))
    write.csv(icc_age_df, icc_age_file, row.names = FALSE)
  } else {
    placeholder <- data.frame(Note = "No age-stratified ICC results", Age_Groups = "50-59, 60-69, 70-79, 80+")
    write.csv(placeholder, icc_age_file, row.names = FALSE)
  }
  cat("\n  [OK] Saved: ICC_Age_Stratified.csv\n")
}, error = function(e) {
  cat("  [ERROR]:", e$message, "\n")
})

# --- ICC by Sex ---
cat("\n--- ICC by Sex ---\n")

icc_sex_stratified <- list()

for (sg in sex_groups) {
  cat("\n  Sex:", sg, "\n")
  
  df_sex <- pooled_data %>% filter(sex == sg)
  n_sg <- nrow(df_sex)
  
  if (n_sg < 100) {
    cat("    [SKIP] Insufficient sample size\n")
    next
  }
  
  for (i in seq_along(lifestyle_vars)) {
    lv <- lifestyle_vars[i]
    if (!lv %in% names(df_sex)) next
    
    tryCatch({
      formula_str <- paste0(lv, " ~ 1 + (1|cohort)")
      model <- lme4::glmer(as.formula(formula_str), data = df_sex, family = binomial,
                           control = lme4::glmerControl(optimizer = "bobyqa"))
      
      icc_result <- performance::icc(model)
      icc_val <- icc_result$ICC_adjusted
      
      icc_sex_stratified[[paste0(sg, "_", lv)]] <- data.frame(
        Sex = sg,
        Variable = lv,
        Variable_Label = lifestyle_labels[i],
        N = n_sg,
        N_Cohorts = length(unique(df_sex$cohort)),
        ICC = round(icc_val, 4),
        stringsAsFactors = FALSE
      )
      
      cat("    ", lifestyle_labels[i], ": ICC =", round(icc_val, 4), "\n")
    }, error = function(e) {
      cat("    ", lifestyle_labels[i], ": [ERROR]\n")
    })
  }
}

# Combine sex-stratified ICC
if (length(icc_sex_stratified) > 0) {
  icc_sex_df <- bind_rows(icc_sex_stratified)
} else {
  icc_sex_df <- data.frame()
}

# Save ALWAYS
cat("\n=== Sex-Stratified ICC Summary ===\n")
cat("Total rows:", nrow(icc_sex_df), "\n")

tryCatch({
  icc_sex_file <- file.path(OUTPUT_DIR, "ICC_Sex_Stratified.csv")
  if (nrow(icc_sex_df) > 0) {
    print(icc_sex_df %>% select(Sex, Variable_Label, N, ICC))
    write.csv(icc_sex_df, icc_sex_file, row.names = FALSE)
  } else {
    placeholder <- data.frame(Note = "No sex-stratified ICC results", Sex_Groups = "Men, Women")
    write.csv(placeholder, icc_sex_file, row.names = FALSE)
  }
  cat("\n  [OK] Saved: ICC_Sex_Stratified.csv\n")
}, error = function(e) {
  cat("  [ERROR]:", e$message, "\n")
})

# =============================================================================
# PART 6: Create Publication-Ready Tables
# =============================================================================

cat("\n\n")
cat("================================================================\n")
cat("         Creating Publication-Ready Tables                      \n")
cat("================================================================\n")

# -----------------------------------------------------------------------------
# 6.1 Phi Coefficient Table (Publication Format)
# -----------------------------------------------------------------------------

phi_publication <- phi_results_overall %>%
  mutate(
    `Phi (95% CI)` = sprintf("%.3f (%.3f, %.3f)", Phi, CI_95_Lower, CI_95_Upper),
    `χ²` = sprintf("%.2f", Chi_Square),
    `P-value` = P_value_formatted,
    `P-value (Bonf.)` = P_value_Bonf_formatted,
    `Sig.` = Significance
  ) %>%
  select(`Variable 1` = Variable_1, `Variable 2` = Variable_2, 
         N, `Phi (95% CI)`, `χ²`, `P-value`, `P-value (Bonf.)`, `Sig.`)

# -----------------------------------------------------------------------------
# 6.2 ICC Tables (Publication Format)
# -----------------------------------------------------------------------------

# Combine all ICC results
icc_all_combined <- bind_rows(
  icc_lifestyle_results %>% mutate(Category = "Lifestyle Factors"),
  icc_baseline_results %>% mutate(Category = "Baseline Health Status"),
  icc_outcome_results %>% mutate(Category = "PPC-MM Outcomes")
)

icc_publication <- icc_all_combined %>%
  mutate(
    `ICC (95% CI)` = ifelse(
      is.na(CI_95_Lower) | is.na(CI_95_Upper),
      sprintf("%.4f", ICC),
      sprintf("%.4f (%.4f, %.4f)", ICC, CI_95_Lower, CI_95_Upper)
    ),
    `P-value` = P_value_formatted,
    `P-value (Bonf.)` = P_value_Bonf_formatted
  ) %>%
  select(Category, Variable, `ICC (95% CI)`, `P-value`, `P-value (Bonf.)`, Interpretation)

# =============================================================================
# PART 7: Save Results
# =============================================================================

cat("\n--- Saving Phi and ICC analysis results ---\n")

# Debug: Check that all required objects exist
cat("\n  Checking required objects...\n")
cat("    - phi_results_overall:", ifelse(exists("phi_results_overall"), paste("OK,", nrow(phi_results_overall), "rows"), "MISSING"), "\n")
cat("    - phi_matrix_overall:", ifelse(exists("phi_matrix_overall"), "OK", "MISSING"), "\n")
cat("    - phi_all_cohorts:", ifelse(exists("phi_all_cohorts"), paste("OK,", nrow(phi_all_cohorts), "rows"), "MISSING"), "\n")
cat("    - icc_lifestyle_results:", ifelse(exists("icc_lifestyle_results"), paste("OK,", nrow(icc_lifestyle_results), "rows"), "MISSING"), "\n")
cat("    - icc_outcome_results:", ifelse(exists("icc_outcome_results"), paste("OK,", nrow(icc_outcome_results), "rows"), "MISSING"), "\n")
cat("    - icc_baseline_results:", ifelse(exists("icc_baseline_results"), paste("OK,", nrow(icc_baseline_results), "rows"), "MISSING"), "\n")
cat("    - phi_publication:", ifelse(exists("phi_publication"), paste("OK,", nrow(phi_publication), "rows"), "MISSING"), "\n")
cat("    - icc_publication:", ifelse(exists("icc_publication"), paste("OK,", nrow(icc_publication), "rows"), "MISSING"), "\n")

# Prepare interpretation guide
interpretation_guide <- data.frame(
  Term = c(
    "Phi Coefficient",
    "ICC (Intraclass Correlation)",
    "Chi-square (χ²)",
    "Bonferroni Correction",
    "Bootstrap CI",
    "LRT P-value",
    "Interpretation (ICC)"
  ),
  Definition = c(
    "Correlation coefficient for two binary variables (-1 to 1). |Phi| > 0.25 = strong, 0.15-0.25 = moderate, 0.05-0.15 = weak",
    "Proportion of variance attributable to clustering (cohort differences)",
    "Test statistic for independence of two categorical variables",
    "Multiple comparison adjustment: P × number of tests. Controls family-wise error rate",
    "95% confidence interval from 1000 bootstrap resamples",
    "P-value from Likelihood Ratio Test comparing random intercept model to null model",
    "<0.05 = Very low, 0.05-0.10 = Low, 0.10-0.20 = Moderate, 0.20-0.30 = High, >0.30 = Very high"
  ),
  stringsAsFactors = FALSE
)

# Statistical notes
stat_notes <- data.frame(
  Note = c(
    "Phi coefficient P-values are from chi-square tests of independence",
    "ICC P-values are from likelihood ratio tests (rptR package)",
    "Bootstrap confidence intervals based on 1000 resamples",
    paste0("Bonferroni threshold for Phi: ", sprintf("%.4f", 0.05/n_pairs), " (", n_pairs, " comparisons)"),
    paste0("Bonferroni threshold for ICC (Lifestyle): ", sprintf("%.4f", 0.05/n_lifestyle_tests), " (", n_lifestyle_tests, " tests)"),
    paste0("Bonferroni threshold for ICC (Outcomes): ", sprintf("%.4f", 0.05/n_outcome_tests), " (", n_outcome_tests, " tests)"),
    "Significance: **** P<0.0001, *** P<0.001, ** P<0.01, * P<0.05, ns = not significant"
  ),
  stringsAsFactors = FALSE
)

# Prepare results list for Excel (including subgroup analyses)
results_list <- list(
  # Interpretation guide
  "Interpretation_Guide" = interpretation_guide,
  "Statistical_Notes" = stat_notes,
  
  # Phi coefficient results - Overall
  "Phi_Summary" = phi_publication,
  "Phi_Detailed" = phi_results_overall %>%
    select(Variable_1, Variable_2, N, Phi, CI_95_Lower, CI_95_Upper, 
           Chi_Square, P_value, P_value_Bonf, Significance),
  "Phi_Matrix" = as.data.frame(phi_matrix_overall),
  "Phi_Pvalue_Matrix" = as.data.frame(p_matrix_overall),
  "Phi_by_Cohort" = phi_all_cohorts,
  
  # Phi coefficient results - Subgroups
  "Phi_Age_Stratified" = if(exists("phi_age_df") && nrow(phi_age_df) > 0) phi_age_df else data.frame(Note = "No data"),
  "Phi_Sex_Stratified" = if(exists("phi_sex_df") && nrow(phi_sex_df) > 0) phi_sex_df else data.frame(Note = "No data"),
  
  # ICC results - Overall
  "ICC_Summary" = icc_publication,
  "ICC_Lifestyle" = icc_lifestyle_results %>%
    select(Variable, ICC, CI_95_Lower, CI_95_Upper, P_value, P_value_Bonf, Interpretation),
  "ICC_Outcomes" = icc_outcome_results %>%
    select(Variable, ICC, CI_95_Lower, CI_95_Upper, P_value, P_value_Bonf, Interpretation),
  "ICC_Baseline" = icc_baseline_results %>%
    select(Variable, ICC, CI_95_Lower, CI_95_Upper, P_value, P_value_Bonf, Interpretation),
  
  # ICC results - Subgroups
  "ICC_Age_Stratified" = if(exists("icc_age_df") && nrow(icc_age_df) > 0) icc_age_df else data.frame(Note = "No data"),
  "ICC_Sex_Stratified" = if(exists("icc_sex_df") && nrow(icc_sex_df) > 0) icc_sex_df else data.frame(Note = "No data")
)

# Save to Excel with error handling
tryCatch({
  save_to_excel(results_list, "Phi_ICC_Analysis_Results.xlsx")
  cat("  [OK] Excel file saved successfully!\n")
}, error = function(e) {
  cat("  [ERROR] Failed to save Excel:", e$message, "\n")
  cat("  Attempting to save individual sheets as CSV...\n")
})

# Save CSV for easy integration
tryCatch({
  phi_csv <- phi_results_overall %>%
    select(Variable_1, Variable_2, N, Phi, CI_95_Lower, CI_95_Upper, 
           Chi_Square, P_value, P_value_Bonf, Significance)
  save_to_csv(phi_csv, "Phi_Coefficients_Complete.csv")
  cat("  [OK] Phi CSV saved!\n")
}, error = function(e) {
  cat("  [ERROR] Failed to save Phi CSV:", e$message, "\n")
})

tryCatch({
  icc_csv <- icc_all_combined %>%
    select(Category, Variable, ICC, CI_95_Lower, CI_95_Upper, P_value, P_value_Bonf, Interpretation)
  save_to_csv(icc_csv, "ICC_Results_Complete.csv")
  cat("  [OK] ICC CSV saved!\n")
}, error = function(e) {
  cat("  [ERROR] Failed to save ICC CSV:", e$message, "\n")
})

# Save as RDS for later use
tryCatch({
  saveRDS(list(
    phi_overall = phi_results_overall,
    phi_matrix = phi_matrix_overall,
    phi_pvalue_matrix = p_matrix_overall,
    phi_by_cohort = phi_all_cohorts,
    icc_lifestyle = icc_lifestyle_results,
    icc_outcomes = icc_outcome_results,
    icc_baseline = icc_baseline_results,
    icc_all = icc_all_combined
  ), file.path(OUTPUT_DIR, "Phi_ICC_results.rds"))
  cat("  [OK] RDS file saved!\n")
}, error = function(e) {
  cat("  [ERROR] Failed to save RDS:", e$message, "\n")
})

# =============================================================================
# PART 8: Visualization
# =============================================================================

cat("\n--- Generating correlation heatmap with P-values ---\n")

# Close any open devices
while (dev.cur() != 1) {
  dev.off()
}

tryCatch({
  heatmap_pdf <- file.path(OUTPUT_DIR, "Phi_Correlation_Heatmap.pdf")
  heatmap_png <- file.path(OUTPUT_DIR, "Phi_Correlation_Heatmap.png")
  
  # PDF version
  pdf(heatmap_pdf, width = 9, height = 7, useDingbats = FALSE, onefile = TRUE)
  
  # Create custom corrplot with significance markers
  corrplot(phi_matrix_overall, 
           method = "color",
           type = "lower",
           addCoef.col = "black",
           number.cex = 1.0,
           tl.col = "black",
           tl.srt = 45,
           tl.cex = 0.9,
           title = "Phi Coefficients between Lifestyle Variables\n(Bonferroni-corrected P-values shown)",
           mar = c(0, 0, 3, 0),
           col = colorRampPalette(c("#2166AC", "#67A9CF", "#D1E5F0", 
                                    "white", 
                                    "#FDDBC7", "#EF8A62", "#B2182B"))(100),
           diag = FALSE)
  
  # Add significance markers in upper triangle
  for (i in 1:(nrow(p_matrix_overall) - 1)) {
    for (j in (i + 1):ncol(p_matrix_overall)) {
      p <- p_matrix_overall[i, j]
      if (!is.na(p)) {
        sig <- ifelse(p < 0.0001, "****",
                      ifelse(p < 0.001, "***",
                             ifelse(p < 0.01, "**",
                                    ifelse(p < 0.05, "*", ""))))
        if (sig != "") {
          # Note: corrplot coordinates are flipped
          text(j, nrow(p_matrix_overall) - i + 1, sig, cex = 0.8, col = "red")
        }
      }
    }
  }
  
  dev.off()
  
  if (file.exists(heatmap_pdf) && file.info(heatmap_pdf)$size > 100) {
    cat("  [OK] Heatmap PDF saved (", round(file.info(heatmap_pdf)$size/1024, 1), "KB)\n")
  }
  
  # PNG version
  png(heatmap_png, width = 9, height = 7, units = "in", res = 300)
  
  corrplot(phi_matrix_overall, 
           method = "color",
           type = "lower",
           addCoef.col = "black",
           number.cex = 1.0,
           tl.col = "black",
           tl.srt = 45,
           tl.cex = 0.9,
           title = "Phi Coefficients between Lifestyle Variables\n(Bonferroni-corrected P-values shown)",
           mar = c(0, 0, 3, 0),
           col = colorRampPalette(c("#2166AC", "#67A9CF", "#D1E5F0", 
                                    "white", 
                                    "#FDDBC7", "#EF8A62", "#B2182B"))(100),
           diag = FALSE)
  
  dev.off()
  
  if (file.exists(heatmap_png) && file.info(heatmap_png)$size > 100) {
    cat("  [OK] Heatmap PNG saved (", round(file.info(heatmap_png)$size/1024, 1), "KB)\n")
  }
  
}, error = function(e) {
  cat("  [ERROR] Heatmap generation failed:", e$message, "\n")
  while (dev.cur() != 1) {
    tryCatch(dev.off(), error = function(e) break)
  }
})

# =============================================================================
# PART 9: Summary
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("         Phi Coefficient and ICC Analysis Complete              \n")
cat("================================================================\n\n")

cat("SUMMARY OF RESULTS:\n")
cat("-------------------\n\n")

cat("1. Phi Coefficient Analysis:\n")
cat("   - ", nrow(phi_results_overall), " pairwise correlations analyzed\n")
cat("   - Bonferroni-corrected significance threshold: ", sprintf("%.4f", 0.05/n_pairs), "\n")
n_sig <- sum(phi_results_overall$P_value_Bonf < 0.05, na.rm = TRUE)
cat("   - Significant correlations: ", n_sig, "/", nrow(phi_results_overall), "\n\n")

cat("2. ICC Analysis (Lifestyle Variables):\n")
for (i in 1:nrow(icc_lifestyle_results)) {
  row <- icc_lifestyle_results[i, ]
  cat(sprintf("   - %s: ICC = %.4f, P = %s (%s)\n",
              row$Variable, row$ICC, row$P_value_Bonf_formatted, row$Interpretation))
}

cat("\n3. ICC Analysis (PPC-MM Outcomes):\n")
for (i in 1:nrow(icc_outcome_results)) {
  row <- icc_outcome_results[i, ]
  cat(sprintf("   - %s: ICC = %.4f, P = %s (%s)\n",
              row$Variable, row$ICC, row$P_value_Bonf_formatted, row$Interpretation))
}

cat("\n4. Output Files:\n")
cat("   - Phi_ICC_Analysis_Results.xlsx (comprehensive Excel workbook)\n")
cat("   - Phi_Coefficients_Complete.csv (for supplementary materials)\n")
cat("   - ICC_Results_Complete.csv (for supplementary materials)\n")
cat("   - Phi_Correlation_Heatmap.pdf/png\n")

cat("\nResults saved to:", OUTPUT_DIR, "\n")

# =============================================================================
# CLEANUP: Stop Parallel Cluster
# =============================================================================

cat("\n--- Stopping parallel cluster ---\n")
tryCatch({
  stopCluster(cl)
  cat("  Parallel cluster stopped successfully.\n")
}, error = function(e) {
  cat("  Note: Cluster cleanup:", e$message, "\n")
})

cat("\n================================================================\n")
cat("         Analysis Complete!                                      \n")
cat("================================================================\n")
