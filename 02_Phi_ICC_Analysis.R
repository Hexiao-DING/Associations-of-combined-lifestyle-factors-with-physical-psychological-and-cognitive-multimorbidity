###############################################################################
# 02_Phi_ICC_Analysis.R
# ============================================================================
# Project: Lifestyle Factors and PPC-MM Association Study
# Purpose: Phi Coefficient and ICC (Intraclass Correlation) Analysis
# ============================================================================
# This script includes:
#   1. Phi coefficient calculation - correlation between lifestyle variables
#      - Cumulative effect correlation (0 vs 2 vs 3 vs 4 unhealthy lifestyles)
#      - Specific effect correlation (pairwise lifestyle factors)
#      - Stratified by age groups
#   2. ICC calculation - assess clustering by cohort
#      - Lifestyle variables ICC
#      - PPC-MM outcome ICC
###############################################################################

# =============================================================================
# PART 1: Environment Setup
# =============================================================================

source("C:/Users/user/Desktop/Scientific Research/UKB Project/Jung Sun Lab/Project1_Lifestyle_PPM_Hexiao,Hongtao/Code/00_Functions_and_Setup.R")

cat("\n")
cat("================================================================\n")
cat("         Phi Coefficient and ICC Analysis                       \n")
cat("================================================================\n\n")

# =============================================================================
# PART 2: Load Pooled Data
# =============================================================================

cat("--- Loading pooled data ---\n")

pooled_data <- readRDS(file.path(OUTPUT_DIR, "Pooled_main_data.rds"))

cat("Pooled data loaded, sample size:", nrow(pooled_data), "\n")

# =============================================================================
# PART 3: Phi Coefficient - Pairwise Lifestyle Correlation
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("         Phi Coefficient Analysis: Lifestyle Correlations       \n")
cat("================================================================\n")

# Lifestyle variables
lifestyle_vars <- c("unhealthy_drink", "unhealthy_smoke", 
                    "unhealthy_pa", "unhealthy_soc")

lifestyle_labels <- c("Drinking" = "unhealthy_drink",
                      "Smoking" = "unhealthy_smoke",
                      "Physical Inactivity" = "unhealthy_pa",
                      "Social Isolation" = "unhealthy_soc")

# -----------------------------------------------------------------------------
# 3.1 Overall Phi Coefficient Matrix
# -----------------------------------------------------------------------------

cat("\n=== Overall Phi Coefficient Matrix (Pooled Data) ===\n")

phi_matrix_overall <- matrix(NA, 
                             nrow = length(lifestyle_vars), 
                             ncol = length(lifestyle_vars),
                             dimnames = list(names(lifestyle_labels), 
                                             names(lifestyle_labels)))

for (i in 1:length(lifestyle_vars)) {
  for (j in 1:length(lifestyle_vars)) {
    if (i == j) {
      phi_matrix_overall[i, j] <- 1.000
    } else {
      phi_val <- calculate_phi(
        pooled_data[[lifestyle_vars[i]]],
        pooled_data[[lifestyle_vars[j]]]
      )
      phi_matrix_overall[i, j] <- round(phi_val, 3)
    }
  }
}

cat("\nLifestyle variable Phi coefficient matrix:\n")
print(phi_matrix_overall)

# Correlation strength interpretation
cat("\nCorrelation strength interpretation (Phi > 0.25 = strong correlation):\n")
for (i in 1:(length(lifestyle_vars) - 1)) {
  for (j in (i + 1):length(lifestyle_vars)) {
    phi_val <- phi_matrix_overall[i, j]
    strength <- ifelse(abs(phi_val) > 0.25, "Strong", 
                       ifelse(abs(phi_val) > 0.15, "Moderate",
                              ifelse(abs(phi_val) > 0.05, "Weak", "Very weak")))
    cat(sprintf("  %s vs %s: Phi = %.3f (%s)\n",
                names(lifestyle_labels)[i], 
                names(lifestyle_labels)[j],
                phi_val, strength))
  }
}

# -----------------------------------------------------------------------------
# 3.2 Phi Coefficients Stratified by Age Group
# -----------------------------------------------------------------------------

cat("\n\n=== Phi Coefficients Stratified by Age Group ===\n")

age_groups <- c("50-59", "60-69", "70-79", "80+")

phi_by_age <- list()

for (ag in age_groups) {
  cat("\n--- Age group:", ag, "---\n")
  
  # Filter by age group
  df_age <- pooled_data %>% filter(age_group_4 == ag)
  
  if (nrow(df_age) < 50) {
    cat("Sample size too small (N =", nrow(df_age), "), skipping\n")
    next
  }
  
  cat("Sample size: N =", nrow(df_age), "\n")
  
  # Calculate Phi matrix
  phi_mat_age <- matrix(NA, 
                        nrow = length(lifestyle_vars), 
                        ncol = length(lifestyle_vars),
                        dimnames = list(names(lifestyle_labels), 
                                        names(lifestyle_labels)))
  
  for (i in 1:length(lifestyle_vars)) {
    for (j in 1:length(lifestyle_vars)) {
      if (i == j) {
        phi_mat_age[i, j] <- 1.000
      } else {
        phi_val <- calculate_phi(
          df_age[[lifestyle_vars[i]]],
          df_age[[lifestyle_vars[j]]]
        )
        phi_mat_age[i, j] <- round(phi_val, 3)
      }
    }
  }
  
  phi_by_age[[ag]] <- phi_mat_age
  print(phi_mat_age)
}

# -----------------------------------------------------------------------------
# 3.3 Phi Coefficients Stratified by Cohort
# -----------------------------------------------------------------------------

cat("\n\n=== Phi Coefficients Stratified by Cohort ===\n")

cohort_names <- unique(pooled_data$cohort)

phi_by_cohort <- list()

for (coh in cohort_names) {
  cat("\n--- Cohort:", coh, "---\n")
  
  # Filter by cohort
  df_cohort <- pooled_data %>% filter(cohort == coh)
  
  if (nrow(df_cohort) < 50) {
    cat("Sample size too small (N =", nrow(df_cohort), "), skipping\n")
    next
  }
  
  cat("Sample size: N =", nrow(df_cohort), "\n")
  
  # Calculate Phi matrix
  phi_mat_cohort <- matrix(NA, 
                           nrow = length(lifestyle_vars), 
                           ncol = length(lifestyle_vars),
                           dimnames = list(names(lifestyle_labels), 
                                           names(lifestyle_labels)))
  
  for (i in 1:length(lifestyle_vars)) {
    for (j in 1:length(lifestyle_vars)) {
      if (i == j) {
        phi_mat_cohort[i, j] <- 1.000
      } else {
        phi_val <- calculate_phi(
          df_cohort[[lifestyle_vars[i]]],
          df_cohort[[lifestyle_vars[j]]]
        )
        phi_mat_cohort[i, j] <- round(phi_val, 3)
      }
    }
  }
  
  phi_by_cohort[[coh]] <- phi_mat_cohort
  print(phi_mat_cohort)
  
  # Interpretation
  cat("\nCorrelation strength for", coh, ":\n")
  for (i in 1:(length(lifestyle_vars) - 1)) {
    for (j in (i + 1):length(lifestyle_vars)) {
      phi_val <- phi_mat_cohort[i, j]
      strength <- ifelse(abs(phi_val) > 0.25, "Strong", 
                         ifelse(abs(phi_val) > 0.15, "Moderate",
                                ifelse(abs(phi_val) > 0.05, "Weak", "Very weak")))
      cat(sprintf("  %s vs %s: Phi = %.3f (%s)\n",
                  names(lifestyle_labels)[i], 
                  names(lifestyle_labels)[j],
                  phi_val, strength))
    }
  }
}

# -----------------------------------------------------------------------------
# 3.3 Cumulative Effect Correlation Analysis
# -----------------------------------------------------------------------------

cat("\n\n=== Cumulative Effect Correlation Analysis ===\n")
cat("(Correlation between unhealthy lifestyle count and specific factors)\n")

# Create dummy variables for lifestyle score categories
pooled_data <- pooled_data %>%
  mutate(
    score_0 = as.numeric(unhealthy_score == 0),
    score_1 = as.numeric(unhealthy_score == 1),
    score_2 = as.numeric(unhealthy_score == 2),
    score_3 = as.numeric(unhealthy_score == 3),
    score_4 = as.numeric(unhealthy_score == 4)
  )

score_vars <- c("score_0", "score_1", "score_2", "score_3", "score_4")
score_labels <- c("0 unhealthy", "1 unhealthy", "2 unhealthy", "3 unhealthy", "4 unhealthy")

# Correlation between score categories and specific lifestyle factors
cat("\nPhi correlation between score categories and lifestyle factors:\n")

for (sv in lifestyle_vars) {
  cat("\n", names(lifestyle_labels)[lifestyle_labels == sv], ":\n", sep = "")
  for (i in 1:length(score_vars)) {
    phi_val <- calculate_phi(
      pooled_data[[sv]],
      pooled_data[[score_vars[i]]]
    )
    cat(sprintf("  vs %s: Phi = %.3f\n", score_labels[i], round(phi_val, 3)))
  }
}

# =============================================================================
# PART 4: ICC (Intraclass Correlation Coefficient) Calculation
# =============================================================================

cat("\n\n")
cat("================================================================\n")
cat("         ICC Analysis: Clustering by Cohort                     \n")
cat("================================================================\n")

# -----------------------------------------------------------------------------
# 4.1 ICC for Lifestyle Variables
# -----------------------------------------------------------------------------

cat("\n=== ICC for Lifestyle Variables ===\n")

icc_lifestyle_results <- data.frame(
  Variable = character(),
  ICC_Adjusted = numeric(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

for (lv in lifestyle_vars) {
  cat("\nCalculating ICC for", names(lifestyle_labels)[lifestyle_labels == lv], "...\n")
  
  icc_result <- calculate_icc(lv, "cohort", pooled_data)
  
  # ICC interpretation
  icc_val <- icc_result$ICC_adjusted
  interpretation <- case_when(
    is.na(icc_val) ~ "Cannot calculate",
    icc_val < 0.05 ~ "Very low clustering",
    icc_val < 0.10 ~ "Low clustering",
    icc_val < 0.20 ~ "Moderate clustering",
    icc_val < 0.30 ~ "High clustering",
    TRUE ~ "Very high clustering"
  )
  
  icc_lifestyle_results <- rbind(icc_lifestyle_results, data.frame(
    Variable = names(lifestyle_labels)[lifestyle_labels == lv],
    ICC_Adjusted = round(icc_val, 4),
    Interpretation = interpretation,
    stringsAsFactors = FALSE
  ))
  
  cat("  ICC (Adjusted):", round(icc_val, 4), "-", interpretation, "\n")
}

cat("\nLifestyle variables ICC summary:\n")
print(icc_lifestyle_results)

# -----------------------------------------------------------------------------
# 4.2 ICC for Unhealthy Lifestyle Score
# -----------------------------------------------------------------------------

cat("\n=== ICC for Unhealthy Lifestyle Score ===\n")

icc_score <- calculate_icc("unhealthy_score", "cohort", pooled_data)

cat("Unhealthy lifestyle score (0-4):\n")
cat("  ICC (Adjusted):", round(icc_score$ICC_adjusted, 4), "\n")

# -----------------------------------------------------------------------------
# 4.3 ICC for Outcome Variables (All 5 Outcomes)
# -----------------------------------------------------------------------------

cat("\n=== ICC for Outcome Variables (All 5 Outcomes) ===\n")

icc_outcome_results <- data.frame(
  Variable = character(),
  ICC_Adjusted = numeric(),
  Interpretation = character(),
  stringsAsFactors = FALSE
)

# Include all 5 PPC-MM outcomes
outcome_vars <- c("physical_base", "psych_base", "cog_base", 
                  "event_ppcmm", "event_mm_phys_psych", "event_mm_phys_cog",
                  "event_mm_psych_cog", "event_mm_all_three")
outcome_labels <- c("Physical disease (baseline)", "Depression (baseline)", 
                    "Cognitive impairment (baseline)", 
                    "Overall PPC-MM (Primary)", "P1P2 (Physical-Psychological)",
                    "P1C (Physical-Cognitive)", "P2C (Psychological-Cognitive)",
                    "P1P2C (All Three)")

for (i in 1:length(outcome_vars)) {
  cat("\nCalculating ICC for", outcome_labels[i], "...\n")
  
  icc_result <- calculate_icc(outcome_vars[i], "cohort", pooled_data)
  
  icc_val <- icc_result$ICC_adjusted
  interpretation <- case_when(
    is.na(icc_val) ~ "Cannot calculate",
    icc_val < 0.05 ~ "Very low clustering",
    icc_val < 0.10 ~ "Low clustering",
    icc_val < 0.20 ~ "Moderate clustering",
    icc_val < 0.30 ~ "High clustering",
    TRUE ~ "Very high clustering"
  )
  
  icc_outcome_results <- rbind(icc_outcome_results, data.frame(
    Variable = outcome_labels[i],
    ICC_Adjusted = round(icc_val, 4),
    Interpretation = interpretation,
    stringsAsFactors = FALSE
  ))
  
  cat("  ICC:", round(icc_val, 4), "-", interpretation, "\n")
}

cat("\nOutcome variables ICC summary:\n")
print(icc_outcome_results)

# -----------------------------------------------------------------------------
# 4.4 ICC for Each Individual Cohort
# -----------------------------------------------------------------------------

cat("\n\n=== ICC for Each Individual Cohort ===\n")
cat("(Using age group as clustering variable for within-cohort ICC)\n")

# Cohort settings
cohort_files <- c("CHARLS", "ELSA", "HRS", "SHARE", "MHAS")
baseline_waves <- list(CHARLS = 1, ELSA = 7, HRS = 9, SHARE = 4, MHAS = 3)

# For SHARE, we will use country; for others, use age groups
icc_by_cohort <- list()

for (coh in cohort_files) {
  cat("\n")
  cat("================================================================\n")
  cat("         ICC Analysis for:", coh, "\n")
  cat("================================================================\n")
  
  file_path <- file.path(OUTPUT_DIR, paste0(coh, "_main_data.rds"))
  
  if (!file.exists(file_path)) {
    cat("  ERROR: Data file not found at:", file_path, "\n")
    cat("  Skipping", coh, "\n")
    next
  }
  
  df_cohort <- readRDS(file_path)
  cat("  Sample size: N =", nrow(df_cohort), "\n")
  cat("  Available columns:", paste(head(names(df_cohort), 15), collapse = ", "), "...\n")
  
  bw <- baseline_waves[[coh]]
  
  # Define variable names for this cohort
  lifestyle_vars_cohort <- c(
    paste0("w", bw, "_unhealthy_drink"),
    paste0("w", bw, "_unhealthy_smoke"),
    paste0("w", bw, "_unhealthy_pa"),
    paste0("w", bw, "_unhealthy_soc")
  )
  lifestyle_labels_cohort <- c("Drinking", "Smoking", "Physical Inactivity", "Social Isolation")
  score_var <- paste0("w", bw, "_unhealthy_score")
  age_var <- paste0("w", bw, "_age")
  
  # -------------------------------------------------------------------------
  # Determine clustering variable
  # -------------------------------------------------------------------------
  cluster_var <- NULL
  cluster_type <- NULL
  
  # For SHARE: prefer country
  if (coh == "SHARE" && "country" %in% names(df_cohort)) {
    cluster_var <- "country"
    cluster_type <- "Country"
    n_clusters <- length(unique(na.omit(df_cohort[[cluster_var]])))
    cat("  Clustering variable: country (", n_clusters, " countries)\n")
  }
  
  # If no cluster_var yet, try to create age groups
  if (is.null(cluster_var)) {
    # Find age variable
    if (age_var %in% names(df_cohort)) {
      df_cohort$age_group_icc <- cut(df_cohort[[age_var]], 
                                     breaks = c(0, 60, 70, 80, Inf),
                                     labels = c("50-59", "60-69", "70-79", "80+"),
                                     right = FALSE)
      cluster_var <- "age_group_icc"
      cluster_type <- "Age Group"
      n_clusters <- length(unique(na.omit(df_cohort[[cluster_var]])))
      cat("  Clustering variable: age_group (", n_clusters, " groups)\n")
    } else {
      # Try alternative age variable names
      alt_age_vars <- c("age", "age_baseline", "baseline_age")
      for (av in alt_age_vars) {
        if (av %in% names(df_cohort)) {
          df_cohort$age_group_icc <- cut(df_cohort[[av]], 
                                         breaks = c(0, 60, 70, 80, Inf),
                                         labels = c("50-59", "60-69", "70-79", "80+"),
                                         right = FALSE)
          cluster_var <- "age_group_icc"
          cluster_type <- "Age Group"
          n_clusters <- length(unique(na.omit(df_cohort[[cluster_var]])))
          cat("  Using", av, "for age grouping (", n_clusters, " groups)\n")
          break
        }
      }
    }
  }
  
  # If still no cluster variable, skip this cohort
  if (is.null(cluster_var)) {
    cat("  WARNING: No suitable clustering variable found for", coh, "\n")
    cat("  Available variables:", paste(names(df_cohort), collapse = ", "), "\n")
    next
  }
  
  # Check if we have enough clusters
  if (n_clusters < 2) {
    cat("  WARNING: Only", n_clusters, "cluster(s), cannot calculate ICC\n")
    next
  }
  
  # -------------------------------------------------------------------------
  # Calculate ICC for this cohort
  # -------------------------------------------------------------------------
  
  icc_cohort_results <- data.frame(
    Cohort = character(),
    Cluster_Type = character(),
    N_Clusters = integer(),
    Variable = character(),
    ICC_Adjusted = numeric(),
    Interpretation = character(),
    stringsAsFactors = FALSE
  )
  
  # Helper function to add ICC result
  add_icc_result <- function(var_name, var_label) {
    if (!var_name %in% names(df_cohort)) {
      cat("    Variable", var_name, "not found\n")
      return(NULL)
    }
    
    result <- tryCatch({
      icc_res <- calculate_icc(var_name, cluster_var, df_cohort)
      icc_val <- icc_res$ICC_adjusted
      
      interp <- case_when(
        is.na(icc_val) ~ "Cannot calculate",
        icc_val < 0.05 ~ "Very low",
        icc_val < 0.10 ~ "Low",
        icc_val < 0.20 ~ "Moderate",
        icc_val < 0.30 ~ "High",
        TRUE ~ "Very high"
      )
      
      cat("    ", var_label, ": ICC =", sprintf("%.4f", icc_val), "(", interp, ")\n")
      
      data.frame(
        Cohort = coh,
        Cluster_Type = cluster_type,
        N_Clusters = n_clusters,
        Variable = var_label,
        ICC_Adjusted = round(icc_val, 4),
        Interpretation = interp,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      cat("    Error for", var_label, ":", e$message, "\n")
      return(NULL)
    })
    
    return(result)
  }
  
  cat("\n  --- Lifestyle Variables ICC ---\n")
  for (i in seq_along(lifestyle_vars_cohort)) {
    res <- add_icc_result(lifestyle_vars_cohort[i], lifestyle_labels_cohort[i])
    if (!is.null(res)) {
      icc_cohort_results <- rbind(icc_cohort_results, res)
    }
  }
  
  cat("\n  --- Unhealthy Score ICC ---\n")
  res <- add_icc_result(score_var, "Unhealthy Score (0-4)")
  if (!is.null(res)) {
    icc_cohort_results <- rbind(icc_cohort_results, res)
  }
  
  cat("\n  --- Outcome Variables ICC ---\n")
  res <- add_icc_result("event_ppcmm", "PPC-MM Event")
  if (!is.null(res)) {
    icc_cohort_results <- rbind(icc_cohort_results, res)
  }
  
  # Store results for this cohort
  if (nrow(icc_cohort_results) > 0) {
    icc_by_cohort[[coh]] <- icc_cohort_results
    cat("\n  Total ICC results for", coh, ":", nrow(icc_cohort_results), "variables\n")
  } else {
    cat("\n  WARNING: No ICC results generated for", coh, "\n")
  }
}

# -------------------------------------------------------------------------
# Combine all cohort ICC results
# -------------------------------------------------------------------------

cat("\n\n================================================================\n")
cat("         Summary: ICC Results by Cohort                         \n")
cat("================================================================\n")

if (length(icc_by_cohort) > 0) {
  icc_all_cohorts <- bind_rows(icc_by_cohort)
  
  cat("\nTotal cohorts with ICC results:", length(icc_by_cohort), "\n")
  cat("Cohorts:", paste(names(icc_by_cohort), collapse = ", "), "\n")
  cat("Total rows:", nrow(icc_all_cohorts), "\n\n")
  
  # Print summary by cohort
  for (coh in names(icc_by_cohort)) {
    cat("\n--- ", coh, " ---\n")
    print(icc_by_cohort[[coh]] %>% select(Variable, ICC_Adjusted, Interpretation))
  }
} else {
  icc_all_cohorts <- data.frame()
  cat("\nWARNING: No cohort-specific ICC results generated!\n")
  cat("Please check that cohort data files exist in:", OUTPUT_DIR, "\n")
}

# -----------------------------------------------------------------------------
# 4.5 ICC Interpretation
# -----------------------------------------------------------------------------

cat("\n\n=== ICC Interpretation Guide ===\n")
cat("
ICC (Intraclass Correlation Coefficient) measures the degree of similarity 
among participants within the same cluster (cohort or country).

Interpretation:
  - ICC < 0.05: Very low clustering, cluster effect is negligible
  - ICC 0.05-0.10: Low clustering, minor cluster effect
  - ICC 0.10-0.20: Moderate clustering, some cluster effect present
  - ICC 0.20-0.30: High clustering, notable cluster effect
  - ICC > 0.30: Very high clustering, significant heterogeneity between clusters

Implications for this study:
  - Pooled ICC: Measures heterogeneity across different cohorts
  - SHARE ICC (by country): Measures heterogeneity across European countries
  - Other cohort ICC: Measures within-cohort clustering (if applicable)
  - High ICC suggests need for multilevel modeling or stratified analysis
")

# =============================================================================
# PART 5: Save Results
# =============================================================================

cat("\n--- Saving Phi and ICC analysis results ---\n")

# Prepare results list
phi_results_list <- list(
  # Phi coefficient results
  "Phi_Overall" = as.data.frame(phi_matrix_overall),
  "Phi_Age_50_59" = if (!is.null(phi_by_age[["50-59"]])) as.data.frame(phi_by_age[["50-59"]]) else data.frame(),
  "Phi_Age_60_69" = if (!is.null(phi_by_age[["60-69"]])) as.data.frame(phi_by_age[["60-69"]]) else data.frame(),
  "Phi_Age_70_79" = if (!is.null(phi_by_age[["70-79"]])) as.data.frame(phi_by_age[["70-79"]]) else data.frame(),
  "Phi_Age_80_plus" = if (!is.null(phi_by_age[["80+"]])) as.data.frame(phi_by_age[["80+"]]) else data.frame(),
  "Phi_CHARLS" = if (!is.null(phi_by_cohort[["CHARLS"]])) as.data.frame(phi_by_cohort[["CHARLS"]]) else data.frame(),
  "Phi_ELSA" = if (!is.null(phi_by_cohort[["ELSA"]])) as.data.frame(phi_by_cohort[["ELSA"]]) else data.frame(),
  "Phi_HRS" = if (!is.null(phi_by_cohort[["HRS"]])) as.data.frame(phi_by_cohort[["HRS"]]) else data.frame(),
  "Phi_SHARE" = if (!is.null(phi_by_cohort[["SHARE"]])) as.data.frame(phi_by_cohort[["SHARE"]]) else data.frame(),
  "Phi_MHAS" = if (!is.null(phi_by_cohort[["MHAS"]])) as.data.frame(phi_by_cohort[["MHAS"]]) else data.frame(),
  # ICC results - Pooled data
  "ICC_Pooled_Lifestyle" = icc_lifestyle_results,
  "ICC_Pooled_Outcomes" = icc_outcome_results,
  # ICC results - All cohorts combined
  "ICC_All_Cohorts" = if (exists("icc_all_cohorts") && nrow(icc_all_cohorts) > 0) icc_all_cohorts else data.frame(),
  # ICC results - Individual cohorts
  "ICC_CHARLS" = if (!is.null(icc_by_cohort[["CHARLS"]])) icc_by_cohort[["CHARLS"]] else data.frame(),
  "ICC_ELSA" = if (!is.null(icc_by_cohort[["ELSA"]])) icc_by_cohort[["ELSA"]] else data.frame(),
  "ICC_HRS" = if (!is.null(icc_by_cohort[["HRS"]])) icc_by_cohort[["HRS"]] else data.frame(),
  "ICC_SHARE" = if (!is.null(icc_by_cohort[["SHARE"]])) icc_by_cohort[["SHARE"]] else data.frame(),
  "ICC_MHAS" = if (!is.null(icc_by_cohort[["MHAS"]])) icc_by_cohort[["MHAS"]] else data.frame()
)

# Save to Excel
save_to_excel(phi_results_list, "Phi_ICC_Analysis_Results.xlsx")

# Save as RDS for later use
saveRDS(list(
  phi_overall = phi_matrix_overall,
  phi_by_age = phi_by_age,
  phi_by_cohort = phi_by_cohort,
  icc_pooled_lifestyle = icc_lifestyle_results,
  icc_pooled_outcomes = icc_outcome_results,
  icc_all_cohorts = if (exists("icc_all_cohorts")) icc_all_cohorts else data.frame(),
  icc_by_cohort_list = if (exists("icc_by_cohort")) icc_by_cohort else list()
), file.path(OUTPUT_DIR, "Phi_ICC_results.rds"))

# =============================================================================
# PART 6: Visualization (Optional)
# =============================================================================

cat("\n--- Generating correlation heatmap ---\n")

# Close any open devices
while (dev.cur() != 1) {
  dev.off()
}

tryCatch({
  heatmap_pdf <- file.path(OUTPUT_DIR, "Phi_Correlation_Heatmap.pdf")
  heatmap_png <- file.path(OUTPUT_DIR, "Phi_Correlation_Heatmap.png")
  
  # PDF version
  pdf(heatmap_pdf, width = 8, height = 6, useDingbats = FALSE, onefile = TRUE)
  
  if (dev.cur() == 1) {
    stop("Failed to open PDF device")
  }
  
  corrplot(phi_matrix_overall, 
           method = "color",
           type = "upper",
           addCoef.col = "black",
           tl.col = "black",
           tl.srt = 45,
           title = "Phi Coefficients between Lifestyle Variables (Pooled)",
           mar = c(0, 0, 2, 0),
           col = colorRampPalette(c("#4575B4", "white", "#D73027"))(100))
  
  dev.off()
  
  Sys.sleep(0.1)
  if (file.exists(heatmap_pdf) && file.info(heatmap_pdf)$size > 100) {
    cat("  [OK] Heatmap PDF saved (", round(file.info(heatmap_pdf)$size/1024, 1), "KB)\n")
  } else {
    cat("  [ERROR] Heatmap PDF not created\n")
  }
  
  # PNG version
  png(heatmap_png, width = 8, height = 6, units = "in", res = 300)
  
  if (dev.cur() == 1) {
    stop("Failed to open PNG device")
  }
  
  corrplot(phi_matrix_overall, 
           method = "color",
           type = "upper",
           addCoef.col = "black",
           tl.col = "black",
           tl.srt = 45,
           title = "Phi Coefficients between Lifestyle Variables (Pooled)",
           mar = c(0, 0, 2, 0),
           col = colorRampPalette(c("#4575B4", "white", "#D73027"))(100))
  
  dev.off()
  
  Sys.sleep(0.1)
  if (file.exists(heatmap_png) && file.info(heatmap_png)$size > 100) {
    cat("  [OK] Heatmap PNG saved (", round(file.info(heatmap_png)$size/1024, 1), "KB)\n")
  } else {
    cat("  [ERROR] Heatmap PNG not created\n")
  }
  
}, error = function(e) {
  cat("  [ERROR] Heatmap generation failed:", e$message, "\n")
  while (dev.cur() != 1) {
    tryCatch(dev.off(), error = function(e) break)
  }
})

# =============================================================================
# COMPLETE
# =============================================================================

cat("\n================================================================\n")
cat("         Phi Coefficient and ICC Analysis Complete              \n")
cat("================================================================\n")
cat("Results saved to:", OUTPUT_DIR, "\n")
