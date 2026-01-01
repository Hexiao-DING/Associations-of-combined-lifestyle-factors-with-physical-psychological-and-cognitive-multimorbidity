###############################################################################
# 08_Methods_Parameters.R
# ============================================================================
# Project: Lifestyle Factors and PPC-MM Association Study
# Purpose: Output all methodology parameters for Methods section
# ============================================================================
#
# OUTPUT:
#   Methods_Parameters.xlsx with sheets:
#   - Study_Overview: Study design, data sources, sample sizes
#   - Cohort_Details: Baseline waves, countries, regions, sample sizes
#   - MICE_Parameters: 2-level MICE imputation parameters
#   - Cox_Parameters: Cox regression covariates, stratification, model formula
#   - Exposure_Definitions: Exposure variable coding definitions
#   - Outcome_Definitions: Outcome variable definitions
#   - PhiICC_Parameters: Phi and ICC analysis parameters
#   - Meta_Parameters: Meta-analysis parameters
#   - PAF_Parameters: PAF calculation parameters
#   - Sensitivity_Analyses: List of sensitivity analyses
#
###############################################################################

# =============================================================================
# PART 1: Environment Setup
# =============================================================================

source("C:/Users/user/Desktop/Scientific Research/UKB Project/Jung Sun Lab/Project1_Lifestyle_PPM_Hexiao,Hongtao/Code/00_Functions_and_Setup.R")

cat("\n")
cat("================================================================\n")
cat("   Methods Parameters Export                                    \n")
cat("================================================================\n\n")

# =============================================================================
# PART 2: Study Overview
# =============================================================================

study_overview <- data.frame(
  Parameter = c(
    "Study Title",
    "Study Design",
    "Data Sources",
    "Number of Cohorts",
    "Total Countries",
    "Total Regions",
    "Baseline Period",
    "Follow-up Period",
    "Primary Outcome",
    "Secondary Outcomes",
    "Primary Exposure",
    "Software",
    "Analysis Date"
  ),
  Value = c(
    "Lifestyle Factors and Physical-Psychological-Cognitive Multi-Morbidity (PPC-MM) Association Study",
    "Prospective cohort study (pooled analysis of harmonized data)",
    "CHARLS, ELSA, HRS, MHAS, SHARE",
    "5",
    "20",
    "7",
    "Varies by cohort (2010-2015)",
    "Up to 10 years",
    "Overall PPC-MM incidence",
    "P1P2 (Physical-Psychological), P1C (Physical-Cognitive), P2C (Psychological-Cognitive), P1P2C (All three)",
    "Cumulative unhealthy lifestyle score (0-4)",
    "R version 4.x with survival, mice, meta packages",
    format(Sys.Date(), "%Y-%m-%d")
  ),
  stringsAsFactors = FALSE
)

# =============================================================================
# PART 3: Cohort Details
# =============================================================================

cohort_details <- data.frame(
  Cohort = c("CHARLS", "ELSA", "HRS", "MHAS", "SHARE"),
  Full_Name = c(
    "China Health and Retirement Longitudinal Study",
    "English Longitudinal Study of Ageing",
    "Health and Retirement Study",
    "Mexican Health and Aging Study",
    "Survey of Health, Ageing and Retirement in Europe"
  ),
  Country = c("China", "England", "United States", "Mexico", "16 European countries"),
  Region = c("Eastern Asia", "Northern Europe", "Northern America", "Central America", "Western/Northern/Southern/Eastern Europe"),
  Baseline_Wave = c(1, 7, 10, 3, 4),
  Baseline_Year = c("2011", "2014", "2010", "2012", "2011"),
  N_Countries = c(1, 1, 1, 1, 16),
  stringsAsFactors = FALSE
)

# Add sample sizes from data
for (i in 1:nrow(cohort_details)) {
  coh <- cohort_details$Cohort[i]
  if (file.exists(DATA_PATHS[[coh]])) {
    df <- read.csv(DATA_PATHS[[coh]], stringsAsFactors = FALSE)
    cohort_details$N_Main[i] <- nrow(df)
  } else {
    cohort_details$N_Main[i] <- NA
  }
  if (file.exists(MICE_PATHS[[coh]])) {
    df <- read.csv(MICE_PATHS[[coh]], stringsAsFactors = FALSE)
    cohort_details$N_MICE[i] <- nrow(df)
  } else {
    cohort_details$N_MICE[i] <- NA
  }
}

# =============================================================================
# PART 4: MICE Parameters
# =============================================================================

mice_parameters <- data.frame(
  Parameter = c(
    "Method",
    "Level 2 Clustering Variable",
    "Level 1",
    "Number of Imputations (m)",
    "Maximum Iterations (maxit)",
    "Random Seed",
    "R Package",
    "Imputation Method",
    "Variables to Impute",
    "Predictor Variables",
    "Auxiliary Variables",
    "Excluded Variables",
    "Number of Level 2 Clusters",
    "Analysis Role",
    "Diagnostics",
    "Pooling Method"
  ),
  Value = c(
    "2-Level Multiple Imputation by Chained Equations (MICE)",
    "cohort",
    "Individual",
    "20",
    "30",
    "12345",
    "mice + miceadds",
    "2l.pmm (Multilevel Predictive Mean Matching)",
    "edu (education level)",
    "age_baseline, sex, marital_binary, physical_base, psych_base, cog_base, unhealthy_score",
    "event_ppcmm, time_ppcmm_months",
    "rural (ELSA 100% missing), employment (excluded by design)",
    "5 (CHARLS, ELSA, HRS, SHARE, MHAS)",
    "Sensitivity Analysis S3",
    "Convergence trace plots, density comparison plots, summary statistics",
    "Rubin's rules for combining multiple imputation results"
  ),
  Description = c(
    "Multilevel imputation accounting for cohort-level clustering",
    "5 cohorts used as Level 2 clustering units",
    "Individual participants within cohorts",
    "Standard recommendation for complex analyses",
    "Sufficient for convergence assessment",
    "For reproducibility of imputation results",
    "R packages for multilevel imputation",
    "Preserves distribution of categorical variables",
    "Education level has ~1-9% missing by cohort",
    "Complete variables used to predict missing edu values",
    "Outcome variables included to satisfy MAR assumption",
    "Variables not imputed due to high missingness or study design",
    "Cohorts as natural clustering structure",
    "MICE results used for sensitivity analysis, not primary analysis",
    "Visual inspection of MCMC convergence and imputation quality",
    "Combined HR estimates across m=20 imputations"
  ),
  stringsAsFactors = FALSE
)

# =============================================================================
# PART 5: Cox Model Parameters
# =============================================================================

cox_parameters <- data.frame(
  Parameter = c(
    "Model Type",
    "Time Variable",
    "Event Variable",
    "Covariates (Pooled Analysis)",
    "Stratification Variable",
    "Covariates (Cohort-specific)",
    "Covariates (SHARE only)",
    "Reference Categories",
    "Confidence Interval",
    "Significance Level",
    "R Package",
    "Model Formula (Pooled)"
  ),
  Value = c(
    "Cox Proportional Hazards Regression",
    "time_ppcmm_months (months from baseline)",
    "event_ppcmm (1 = PPC-MM incident case)",
    paste(COX_COVARIATES$pooled, collapse = " + "),
    paste0("strata(", COX_COVARIATES$strata_var, ")"),
    paste(COX_COVARIATES$cohort_specific, collapse = " + "),
    paste(COX_COVARIATES$share_specific, collapse = " + "),
    "sex=Men, edu=Primary, region=Eastern Asia, n_lifestyle_cat=0",
    "95% CI",
    "α = 0.05 (two-sided)",
    "survival::coxph()",
    "Surv(time, event) ~ exposure + age_baseline + sex + edu + region + strata(cohort)"
  ),
  Description = c(
    "Accounts for time-to-event and censoring",
    "Time from baseline to event or censoring",
    "PPC-MM incidence (≥2 of 3 domains impaired)",
    "Adjusted for demographic and geographic factors",
    "Controls for cohort-specific baseline hazards",
    "For meta-analysis of single-country cohorts",
    "SHARE has multiple countries, adjust for region",
    "Reference levels for categorical variables",
    "Hazard ratio confidence interval",
    "For hypothesis testing",
    "Standard R package for survival analysis",
    "Complete model specification"
  ),
  stringsAsFactors = FALSE
)

# =============================================================================
# PART 6: Exposure Definitions
# =============================================================================

exposure_definitions <- data.frame(
  Variable = c(
    "unhealthy_drink",
    "unhealthy_smoke",
    "unhealthy_pa",
    "unhealthy_soc",
    "heavy_drink",
    "unhealthy_score",
    "unhealthy_score_heavy",
    "n_lifestyle_cat",
    "n_lifestyle_5cat"
  ),
  Definition = c(
    "Any current alcohol consumption",
    "Current smoker",
    "Physical inactivity (no moderate/vigorous activity)",
    "Social isolation (living alone, no social participation)",
    "Heavy drinking (>14 drinks/week men, >7 women)",
    "Sum of 4 unhealthy lifestyle factors (0-4)",
    "Sum using heavy_drink instead of unhealthy_drink",
    "Categorical: 0, 1, 2, 3+ (reference = 0)",
    "Categorical: 0, 1, 2, 3, 4 (reference = 0)"
  ),
  Coding = c(
    "0 = No/Former, 1 = Current drinker",
    "0 = Never/Former, 1 = Current smoker",
    "0 = Active, 1 = Inactive",
    "0 = Not isolated, 1 = Isolated",
    "0 = No/Moderate, 1 = Heavy drinker",
    "Continuous 0-4",
    "Continuous 0-4",
    "Factor with 4 levels",
    "Factor with 5 levels"
  ),
  stringsAsFactors = FALSE
)

# =============================================================================
# PART 7: Outcome Definitions
# =============================================================================

outcome_definitions <- data.frame(
  Outcome = c("Overall", "P1P2", "P1C", "P2C", "P1P2C"),
  Label = c(
    "Overall PPC-MM",
    "Physical-Psychological (P1P2)",
    "Physical-Cognitive (P1C)",
    "Psychological-Cognitive (P2C)",
    "All Three Domains (P1P2C)"
  ),
  Type = c("Primary", "Secondary", "Secondary", "Secondary", "Secondary"),
  Event_Variable = c(
    "event_ppcmm",
    "event_mm_phys_psych",
    "event_mm_phys_cog",
    "event_mm_psych_cog",
    "event_mm_all_three"
  ),
  Time_Variable = c(
    "time_ppcmm_months",
    "time_mm_phys_psych",
    "time_mm_phys_cog",
    "time_mm_psych_cog",
    "time_mm_all_three"
  ),
  Definition = c(
    "≥2 of 3 domains impaired (Physical, Psychological, Cognitive)",
    "Both Physical disease and Psychological condition present",
    "Both Physical disease and Cognitive impairment present",
    "Both Psychological condition and Cognitive impairment present",
    "All three domains impaired"
  ),
  stringsAsFactors = FALSE
)

# =============================================================================
# PART 8: Phi and ICC Parameters
# =============================================================================

phi_icc_parameters <- data.frame(
  Analysis = c(
    "Phi Coefficient",
    "Phi Coefficient",
    "Phi Coefficient",
    "Phi Coefficient",
    "ICC",
    "ICC",
    "ICC",
    "ICC"
  ),
  Parameter = c(
    "Method",
    "Bootstrap Iterations",
    "Confidence Interval",
    "Multiple Comparison Correction",
    "Method",
    "Cluster Variable",
    "P-value Method",
    "CI Method"
  ),
  Value = c(
    "Pearson correlation for binary variables",
    "1000",
    "95% CI (Bootstrap percentile)",
    "Bonferroni (6 pairwise comparisons, α = 0.0083)",
    "Mixed-effects logistic regression",
    "cohort",
    "Likelihood Ratio Test (LRT)",
    "Delta method approximation"
  ),
  stringsAsFactors = FALSE
)

# =============================================================================
# PART 9: Meta-Analysis Parameters
# =============================================================================

meta_parameters <- data.frame(
  Parameter = c(
    "Effect Measure",
    "Fixed-Effects Model",
    "Random-Effects Model",
    "Heterogeneity Statistics",
    "Model Selection Criterion",
    "Sensitivity Analysis",
    "Publication Bias Tests",
    "R Packages"
  ),
  Value = c(
    "Hazard Ratio (HR) with 95% CI",
    "Inverse variance weighting",
    "DerSimonian-Laird (DL) estimator",
    "Cochran's Q, I², τ²",
    "Random-effects if I² > 50%, otherwise Fixed-effects",
    "Leave-one-out analysis",
    "Egger's regression test, Begg's rank correlation test",
    "meta, metafor"
  ),
  Description = c(
    "Primary effect estimate for time-to-event outcomes",
    "Assumes homogeneous effects across studies",
    "Accounts for between-study heterogeneity",
    "Quantifies heterogeneity and between-study variance",
    "Standard threshold for substantial heterogeneity",
    "Assess influence of individual studies",
    "Assess small-study effects/publication bias",
    "R packages for meta-analysis"
  ),
  stringsAsFactors = FALSE
)

# =============================================================================
# PART 10: PAF Parameters
# =============================================================================

paf_parameters <- data.frame(
  Parameter = c(
    "PAF Formula",
    "HR Source",
    "Prevalence Definition",
    "Confidence Interval",
    "Interpretation"
  ),
  Value = c(
    "PAF = P_case × (HR - 1) / HR (Miettinen formula)",
    "Cox regression adjusted HR",
    "Proportion of cases exposed (P_case)",
    "Bootstrap (1000 iterations)",
    "Proportion of cases attributable to exposure if causal"
  ),
  Description = c(
    "Population Attributable Fraction for cohort studies",
    "Hazard ratios from primary Cox analysis",
    "Exposure prevalence among incident cases",
    "Non-parametric percentile method",
    "Assumes causal relationship and no unmeasured confounding"
  ),
  stringsAsFactors = FALSE
)

# =============================================================================
# PART 11: Sensitivity Analyses
# =============================================================================

sensitivity_analyses <- data.frame(
  Analysis_Code = c("S1", "S2", "S3", "S4"),
  Analysis_Name = c(
    "5-Level Lifestyle Categories",
    "Heavy Drinking Definition",
    "MICE Imputed Data",
    "Exclude First Follow-up Wave"
  ),
  Description = c(
    "Use 5-level categories (0/1/2/3/4) instead of 4-level (0/1/2/3+)",
    "Replace any drinking with heavy drinking definition",
    "Use multiply imputed data (m=20) with Rubin's rules",
    "Exclude events in first follow-up wave to reduce reverse causality"
  ),
  Rationale = c(
    "Test dose-response with finer granularity",
    "Alternative exposure definition for alcohol",
    "Account for missing data uncertainty",
    "Address potential reverse causality"
  ),
  stringsAsFactors = FALSE
)

# =============================================================================
# PART 12: Sankey Diagram Parameters
# =============================================================================

sankey_parameters <- data.frame(
  Parameter = c(
    "Baseline States",
    "Follow-up States",
    "Color Scheme",
    "Label Format",
    "R Package"
  ),
  Value = c(
    "Healthy, P1 (Physical), P2 (Psychological), C (Cognitive)",
    "Healthy, P1, P2, C, P1P2, P1C, P2C, P1P2C",
    "Custom 8-color scheme (green=Healthy, red=P1P2C)",
    "State (n=xxx, xx.x%)",
    "ggalluvial"
  ),
  Description = c(
    "Single-condition states at baseline",
    "All possible multi-morbidity states at follow-up",
    "Color-coded by severity",
    "Shows sample size and percentage for each state",
    "R package for alluvial/Sankey diagrams"
  ),
  stringsAsFactors = FALSE
)

# =============================================================================
# PART 13: Save to Excel
# =============================================================================

cat("--- Saving Methods Parameters to Excel ---\n\n")

output_file <- file.path(OUTPUT_DIR, "Methods_Parameters.xlsx")

writexl::write_xlsx(
  list(
    Study_Overview = study_overview,
    Cohort_Details = cohort_details,
    MICE_Parameters = mice_parameters,
    Cox_Parameters = cox_parameters,
    Exposure_Definitions = exposure_definitions,
    Outcome_Definitions = outcome_definitions,
    PhiICC_Parameters = phi_icc_parameters,
    Meta_Parameters = meta_parameters,
    PAF_Parameters = paf_parameters,
    Sensitivity_Analyses = sensitivity_analyses,
    Sankey_Parameters = sankey_parameters
  ),
  output_file
)

cat("  Saved: Methods_Parameters.xlsx\n")
cat("  Location:", output_file, "\n\n")

cat("  Sheets included:\n")
cat("    1. Study_Overview\n")
cat("    2. Cohort_Details\n")
cat("    3. MICE_Parameters\n")
cat("    4. Cox_Parameters\n")
cat("    5. Exposure_Definitions\n")
cat("    6. Outcome_Definitions\n")
cat("    7. PhiICC_Parameters\n")
cat("    8. Meta_Parameters\n")
cat("    9. PAF_Parameters\n")
cat("   10. Sensitivity_Analyses\n")
cat("   11. Sankey_Parameters\n")

cat("\n")
cat("================================================================\n")
cat("   Methods Parameters Export Complete                           \n")
cat("================================================================\n")

