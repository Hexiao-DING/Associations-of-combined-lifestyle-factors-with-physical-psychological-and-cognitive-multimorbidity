###############################################################################
# 05_Subgroup_Analysis.R
# ============================================================================
# Project: Lifestyle Factors and PPC-MM Association Study
# Purpose: Age/Sex stratified Cox and PAF analysis, Interaction tests
# ============================================================================
#
# OUTPUT FILES:
#   - Cox_Age_Stratified.csv
#   - Cox_Sex_Stratified.csv
#   - Cox_Interaction_Analysis.csv
#   - PAF_Age_Stratified_Individual.csv
#   - PAF_Sex_Stratified_Individual.csv
#   - Subgroup_Analysis_Complete.xlsx
#
###############################################################################

# Load environment (do NOT use rm() - breaks main_analysis.R pipeline)
source("C:/Users/user/Desktop/Scientific Research/UKB Project/Jung Sun Lab/Project1_Lifestyle_PPM_Hexiao,Hongtao/Code/00_Functions_and_Setup.R")

library(tidyverse)
library(survival)
library(writexl)

cat("\n")
cat("================================================================\n")
cat("   Subgroup Analysis: Age/Sex Stratification                     \n")
cat("================================================================\n")
cat("\nOutput directory:", OUTPUT_DIR, "\n")

# =============================================================================
# Load Data
# =============================================================================

pooled_rds <- file.path(OUTPUT_DIR, "Pooled_main_data.rds")
pooled_csv <- file.path(OUTPUT_DIR, "Pooled_main_data.csv")

if (file.exists(pooled_rds)) {
  pooled_data <- readRDS(pooled_rds)
  cat("  Loaded from RDS, N =", nrow(pooled_data), "\n")
} else if (file.exists(pooled_csv)) {
  pooled_data <- read.csv(pooled_csv, stringsAsFactors = FALSE)
  cat("  Loaded from CSV, N =", nrow(pooled_data), "\n")
} else {
  stop("Error: Pooled_main_data not found. Run 01_Pooled_Descriptive.R first.")
}

# =============================================================================
# Prepare Analysis Data
# =============================================================================

cat("\n--- Preparing analysis data ---\n")

pooled_analysis <- pooled_data %>%
  filter(!is.na(age_baseline) & !is.na(sex) & !is.na(edu) &
         !is.na(unhealthy_drink) & !is.na(unhealthy_smoke) &
         !is.na(unhealthy_pa) & !is.na(unhealthy_soc) &
         !is.na(event_ppcmm) & !is.na(time_ppcmm_months) &
         time_ppcmm_months > 0) %>%
  mutate(
    sex = factor(sex),
    edu = factor(edu),
    cohort = factor(cohort),
    age_group_4 = case_when(
      age_baseline >= 50 & age_baseline < 60 ~ "50-59",
      age_baseline >= 60 & age_baseline < 70 ~ "60-69",
      age_baseline >= 70 & age_baseline < 80 ~ "70-79",
      age_baseline >= 80 ~ "80+",
      TRUE ~ NA_character_
    ),
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

if ("region" %in% names(pooled_analysis)) {
  pooled_analysis$region <- factor(pooled_analysis$region)
}

cat("  Analysis sample: N =", nrow(pooled_analysis), "\n")
cat("  Total events:", sum(pooled_analysis$event_ppcmm), "\n")

cat("\nAge group distribution:\n")
print(table(pooled_analysis$age_group_4))

cat("\nSex distribution:\n")
print(table(pooled_analysis$sex))

# =============================================================================
# PAF Helper Functions
# =============================================================================

calc_paf <- function(p_case, hr) {
  if (is.na(p_case) || is.na(hr) || p_case == 0) return(NA)
  paf <- p_case * (hr - 1) / hr
  return(paf)
}

calc_paf_with_ci <- function(p_case, hr, hr_lower, hr_upper) {
  paf <- calc_paf(p_case, hr)
  paf_lower <- calc_paf(p_case, hr_upper)
  paf_upper <- calc_paf(p_case, hr_lower)
  return(list(paf = paf, paf_lower = paf_lower, paf_upper = paf_upper))
}

# =============================================================================
# Cox: Age-Stratified Analysis
# =============================================================================

cat("\n--- Cox: Age-Stratified Analysis ---\n")

lf_vars <- c("unhealthy_drink", "unhealthy_smoke", "unhealthy_pa", "unhealthy_soc")
lf_labels <- c("Drinking", "Smoking", "Physical Inactivity", "Social Isolation")

age_groups <- c("50-59", "60-69", "70-79", "80+")
age_results <- list()

for (ag in age_groups) {
  cat("\n  Age Group:", ag, "\n")
  
  df_age <- pooled_analysis %>% filter(age_group_4 == ag)
  n_total <- nrow(df_age)
  n_events <- sum(df_age$event_ppcmm, na.rm = TRUE)
  
  cat("    N =", n_total, ", Events =", n_events, "\n")
  
  if (n_events < 10) {
    cat("    [SKIP] Insufficient events\n")
    next
  }
  
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
    cat("    [ERROR]", e$message, "\n")
    NULL
  })
  
  if (!is.null(fit)) {
    s <- summary(fit)
    for (i in seq_along(lf_vars)) {
      lv <- lf_vars[i]
      if (!lv %in% rownames(s$coefficients)) next
      
      age_results[[paste0(ag, "_", lv)]] <- data.frame(
        Analysis = "Age-Stratified",
        Age_Group = ag,
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
    cat("    [OK]\n")
  }
}

age_stratified_df <- if (length(age_results) > 0) bind_rows(age_results) else data.frame(Note = "No results")

cat("\n  Age-stratified results:", nrow(age_stratified_df), "rows\n")
write.csv(age_stratified_df, file.path(OUTPUT_DIR, "Cox_Age_Stratified.csv"), row.names = FALSE)
cat("  [SAVED] Cox_Age_Stratified.csv\n")

# =============================================================================
# Cox: Sex-Stratified Analysis
# =============================================================================

cat("\n--- Cox: Sex-Stratified Analysis ---\n")

sex_groups <- c("Men", "Women")
sex_results <- list()

for (sg in sex_groups) {
  cat("\n  Sex:", sg, "\n")
  
  df_sex <- pooled_analysis %>% filter(sex == sg)
  n_total <- nrow(df_sex)
  n_events <- sum(df_sex$event_ppcmm, na.rm = TRUE)
  
  cat("    N =", n_total, ", Events =", n_events, "\n")
  
  if (n_events < 10) {
    cat("    [SKIP] Insufficient events\n")
    next
  }
  
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
    cat("    [ERROR]", e$message, "\n")
    NULL
  })
  
  if (!is.null(fit)) {
    s <- summary(fit)
    for (i in seq_along(lf_vars)) {
      lv <- lf_vars[i]
      if (!lv %in% rownames(s$coefficients)) next
      
      sex_results[[paste0(sg, "_", lv)]] <- data.frame(
        Analysis = "Sex-Stratified",
        Sex = sg,
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
    cat("    [OK]\n")
  }
}

sex_stratified_df <- if (length(sex_results) > 0) bind_rows(sex_results) else data.frame(Note = "No results")

cat("\n  Sex-stratified results:", nrow(sex_stratified_df), "rows\n")
write.csv(sex_stratified_df, file.path(OUTPUT_DIR, "Cox_Sex_Stratified.csv"), row.names = FALSE)
cat("  [SAVED] Cox_Sex_Stratified.csv\n")

# =============================================================================
# Cox: Interaction Analysis
# =============================================================================

cat("\n--- Cox: Interaction Analysis ---\n")

interaction_results <- list()

for (i in seq_along(lf_vars)) {
  lv <- lf_vars[i]
  
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
    
    interaction_results[[lv]] <- data.frame(
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

interaction_df <- if (length(interaction_results) > 0) bind_rows(interaction_results) else data.frame(Note = "No results")

cat("\n  Interaction results:", nrow(interaction_df), "rows\n")
write.csv(interaction_df, file.path(OUTPUT_DIR, "Cox_Interaction_Analysis.csv"), row.names = FALSE)
cat("  [SAVED] Cox_Interaction_Analysis.csv\n")

# =============================================================================
# PAF: Age-Stratified Analysis
# =============================================================================

cat("\n--- PAF: Age-Stratified Analysis ---\n")

age_paf_results <- list()

for (ag in age_groups) {
  cat("\n  Age Group:", ag, "\n")
  
  df_age <- pooled_analysis %>% filter(age_group_4 == ag)
  n_total <- nrow(df_age)
  n_events <- sum(df_age$event_ppcmm, na.rm = TRUE)
  
  cat("    N =", n_total, ", Events =", n_events, "\n")
  
  if (n_events < 20) {
    cat("    [SKIP] Insufficient events\n")
    next
  }
  
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
  }, error = function(e) NULL)
  
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
      
      age_paf_results[[paste0(ag, "_", lv)]] <- data.frame(
        Analysis = "Age-Stratified PAF",
        Age_Group = ag,
        Variable = lv,
        Variable_Label = lf_labels[i],
        HR = round(hr, 3),
        HR_CI = sprintf("%.2f (%.2f-%.2f)", hr, hr_lower, hr_upper),
        P_value = p_value,
        P_case = round(p_case, 4),
        P_case_pct = sprintf("%.1f%%", p_case * 100),
        PAF = round(paf_result$paf, 4),
        PAF_pct = sprintf("%.1f%%", paf_result$paf * 100),
        N = n_total,
        N_Events = n_events,
        stringsAsFactors = FALSE
      )
    }
    cat("    [OK]\n")
  }
}

age_paf_df <- if (length(age_paf_results) > 0) bind_rows(age_paf_results) else data.frame(Note = "No results")

cat("\n  Age-stratified PAF results:", nrow(age_paf_df), "rows\n")
write.csv(age_paf_df, file.path(OUTPUT_DIR, "PAF_Age_Stratified_Individual.csv"), row.names = FALSE)
cat("  [SAVED] PAF_Age_Stratified_Individual.csv\n")

# =============================================================================
# PAF: Sex-Stratified Analysis
# =============================================================================

cat("\n--- PAF: Sex-Stratified Analysis ---\n")

sex_paf_results <- list()

for (sg in sex_groups) {
  cat("\n  Sex:", sg, "\n")
  
  df_sex <- pooled_analysis %>% filter(sex == sg)
  n_total <- nrow(df_sex)
  n_events <- sum(df_sex$event_ppcmm, na.rm = TRUE)
  
  cat("    N =", n_total, ", Events =", n_events, "\n")
  
  if (n_events < 20) {
    cat("    [SKIP] Insufficient events\n")
    next
  }
  
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
  }, error = function(e) NULL)
  
  if (!is.null(fit)) {
    s <- summary(fit)
    cases_only <- df_sex %>% filter(event_ppcmm == 1)
    
    for (i in seq_along(lf_vars)) {
      lv <- lf_vars[i]
      if (!lv %in% rownames(s$coefficients)) next
      
      hr <- s$coefficients[lv, "exp(coef)"]
      hr_lower <- s$conf.int[lv, "lower .95"]
      hr_upper <- s$conf.int[lv, "upper .95"]
      p_value <- s$coefficients[lv, "Pr(>|z|)"]
      p_case <- mean(cases_only[[lv]] == 1, na.rm = TRUE)
      
      paf_result <- calc_paf_with_ci(p_case, hr, hr_lower, hr_upper)
      
      sex_paf_results[[paste0(sg, "_", lv)]] <- data.frame(
        Analysis = "Sex-Stratified PAF",
        Sex = sg,
        Variable = lv,
        Variable_Label = lf_labels[i],
        HR = round(hr, 3),
        HR_CI = sprintf("%.2f (%.2f-%.2f)", hr, hr_lower, hr_upper),
        P_value = p_value,
        P_case = round(p_case, 4),
        P_case_pct = sprintf("%.1f%%", p_case * 100),
        PAF = round(paf_result$paf, 4),
        PAF_pct = sprintf("%.1f%%", paf_result$paf * 100),
        N = n_total,
        N_Events = n_events,
        stringsAsFactors = FALSE
      )
    }
    cat("    [OK]\n")
  }
}

sex_paf_df <- if (length(sex_paf_results) > 0) bind_rows(sex_paf_results) else data.frame(Note = "No results")

cat("\n  Sex-stratified PAF results:", nrow(sex_paf_df), "rows\n")
write.csv(sex_paf_df, file.path(OUTPUT_DIR, "PAF_Sex_Stratified_Individual.csv"), row.names = FALSE)
cat("  [SAVED] PAF_Sex_Stratified_Individual.csv\n")

# =============================================================================
# Save Combined Excel
# =============================================================================

cat("\n--- Saving Combined Excel ---\n")

tryCatch({
  excel_list <- list(
    "Cox_Age_Stratified" = age_stratified_df,
    "Cox_Sex_Stratified" = sex_stratified_df,
    "Cox_Interaction" = interaction_df,
    "PAF_Age_Stratified" = age_paf_df,
    "PAF_Sex_Stratified" = sex_paf_df
  )
  
  writexl::write_xlsx(excel_list, file.path(OUTPUT_DIR, "Subgroup_Analysis_Complete.xlsx"))
  cat("  [SAVED] Subgroup_Analysis_Complete.xlsx\n")
}, error = function(e) {
  cat("  [ERROR] Excel save failed:", e$message, "\n")
})

# =============================================================================
# Output Verification
# =============================================================================

cat("\n--- Output Verification ---\n")

output_files <- c(
  "Cox_Age_Stratified.csv",
  "Cox_Sex_Stratified.csv",
  "Cox_Interaction_Analysis.csv",
  "PAF_Age_Stratified_Individual.csv",
  "PAF_Sex_Stratified_Individual.csv",
  "Subgroup_Analysis_Complete.xlsx"
)

for (f in output_files) {
  fp <- file.path(OUTPUT_DIR, f)
  if (file.exists(fp)) {
    info <- file.info(fp)
    cat("  [OK]", f, "-", info$size, "bytes\n")
  } else {
    cat("  [MISSING]", f, "\n")
  }
}

cat("\n>>> Subgroup Analysis Complete\n")
