###############################################################################
# main_analysis.R
# ============================================================================
# Project: Lifestyle Factors and PPC-MM Association Study
# Purpose: Master script to run all analysis (Lancet Standard)
# ============================================================================
#
# SCRIPT STRUCTURE:
# =================
#   00_Functions_and_Setup.R  - Common functions and environment setup
#   01_Pooled_Descriptive.R   - Pooled data generation and descriptive stats
#   02_Phi_ICC_Analysis.R     - Phi coefficient and ICC analysis
#   03_Pooled_Cox_Analysis.R  - Cox regression (pooled)
#   04_Sankey_Diagram.R       - Health state transition diagrams
#   05_PAF_Analysis.R         - Population Attributable Fraction
#   06_Meta_Analysis.R        - Meta-analysis across cohorts
#
# ANALYSIS FRAMEWORK (Lancet Standard):
# =====================================
#
# PRIMARY ANALYSIS:
#   - Outcomes: Overall PPC-MM (Primary) + 4 subtypes (Secondary)
#   - Exposure A: Individual lifestyle factors (mutually adjusted)
#   - Exposure B: Cumulative effect (4-level: 0/1/2/3+)
#
# SENSITIVITY ANALYSES:
#   S1: 5-level cumulative effect (0/1/2/3/4)
#   S2: Heavy drinking definition
#   S3: MICE imputed data
#   S4: Drop1st (exclude first follow-up wave)
#
###############################################################################

# =============================================================================
# PART 1: Setup and Initialization
# =============================================================================

cat("\n")
cat("################################################################\n")
cat("#                                                              #\n")
cat("#     Lifestyle Factors and PPC-MM Association Study           #\n")
cat("#           Lancet Standard Analysis Pipeline                  #\n")
cat("#                                                              #\n")
cat("################################################################\n")
cat("\n")

start_time <- Sys.time()

# Set working directory to Code folder
CODE_DIR <- "C:/Users/user/Desktop/Scientific Research/UKB Project/Jung Sun Lab/Project1_Lifestyle_PPM_Hexiao,Hongtao/Code"
setwd(CODE_DIR)

# Initialize analysis status tracking
analysis_status <- list()

# =============================================================================
# PART 2: Load Functions and Setup
# =============================================================================

cat("\n================================================================\n")
cat("         Loading Functions and Setup                            \n")
cat("================================================================\n")

cat("\n>>> Loading 00_Functions_and_Setup.R...\n")
tryCatch({
  source("00_Functions_and_Setup.R", local = FALSE)
  analysis_status[["00_Setup"]] <- "SUCCESS"
  cat("\n>>> Setup completed successfully!\n")
  
  # Verify data paths exist
  cat("\n--- Verifying data files ---\n")
  for (coh in names(DATA_PATHS)) {
    if (file.exists(DATA_PATHS[[coh]])) {
      cat("  [OK]", coh, "main data found\n")
    } else {
      cat("  [MISSING]", coh, ":", DATA_PATHS[[coh]], "\n")
    }
    # Check MICE data
    if (file.exists(MICE_PATHS[[coh]])) {
      cat("  [OK]", coh, "MICE data found\n")
    }
  }
  
}, error = function(e) {
  analysis_status[["00_Setup"]] <- paste("FAILED:", e$message)
  cat("\n>>> Setup failed:", e$message, "\n")
  stop("Cannot proceed without setup. Please fix the error above.")
})

# =============================================================================
# PART 3: Pooled Data Generation and Descriptive Statistics
# =============================================================================

cat("\n================================================================\n")
cat("         Running Pooled Data Generation                         \n")
cat("================================================================\n")

cat("\n>>> Running 01_Pooled_Descriptive.R...\n")
tryCatch({
  source("01_Pooled_Descriptive.R", local = TRUE)
  analysis_status[["01_Descriptive"]] <- "SUCCESS"
  cat("\n>>> Pooled Descriptive Analysis completed successfully!\n")
}, error = function(e) {
  analysis_status[["01_Descriptive"]] <- paste("FAILED:", e$message)
  cat("\n>>> Pooled Descriptive Analysis failed:", e$message, "\n")
})

# =============================================================================
# PART 4: Phi Coefficient and ICC Analysis
# =============================================================================

cat("\n================================================================\n")
cat("         Running Phi Coefficient and ICC Analysis               \n")
cat("================================================================\n")

cat("\n>>> Running 02_Phi_ICC_Analysis.R...\n")
tryCatch({
  source("02_Phi_ICC_Analysis.R", local = TRUE)
  analysis_status[["02_Phi_ICC"]] <- "SUCCESS"
  cat("\n>>> Phi and ICC Analysis completed successfully!\n")
}, error = function(e) {
  analysis_status[["02_Phi_ICC"]] <- paste("FAILED:", e$message)
  cat("\n>>> Phi and ICC Analysis failed:", e$message, "\n")
})

# =============================================================================
# PART 5: Pooled Cox Regression Analysis
# =============================================================================

cat("\n================================================================\n")
cat("         Running Pooled Cox Regression Analysis                 \n")
cat("================================================================\n")

cat("\n>>> Running 03_Pooled_Cox_Analysis.R...\n")
tryCatch({
  source("03_Pooled_Cox_Analysis.R", local = TRUE)
  analysis_status[["03_Cox"]] <- "SUCCESS"
  cat("\n>>> Pooled Cox Analysis completed successfully!\n")
}, error = function(e) {
  analysis_status[["03_Cox"]] <- paste("FAILED:", e$message)
  cat("\n>>> Pooled Cox Analysis failed:", e$message, "\n")
})

# =============================================================================
# PART 6: Sankey Diagrams
# =============================================================================

cat("\n================================================================\n")
cat("         Running Sankey Diagram Generation                      \n")
cat("================================================================\n")

cat("\n>>> Running 04_Sankey_Diagram.R...\n")
tryCatch({
  source("04_Sankey_Diagram.R", local = TRUE)
  analysis_status[["04_Sankey"]] <- "SUCCESS"
  cat("\n>>> Sankey Diagram Generation completed successfully!\n")
}, error = function(e) {
  analysis_status[["04_Sankey"]] <- paste("FAILED:", e$message)
  cat("\n>>> Sankey Diagram Generation failed:", e$message, "\n")
})

# =============================================================================
# PART 7: PAF Analysis
# =============================================================================

cat("\n================================================================\n")
cat("         Running PAF Analysis                                   \n")
cat("================================================================\n")

cat("\n>>> Running 05_PAF_Analysis.R...\n")
tryCatch({
  source("05_PAF_Analysis.R", local = TRUE)
  analysis_status[["05_PAF"]] <- "SUCCESS"
  cat("\n>>> PAF Analysis completed successfully!\n")
}, error = function(e) {
  analysis_status[["05_PAF"]] <- paste("FAILED:", e$message)
  cat("\n>>> PAF Analysis failed:", e$message, "\n")
})

# =============================================================================
# PART 8: Meta-Analysis
# =============================================================================

cat("\n================================================================\n")
cat("         Running Meta-Analysis                                  \n")
cat("================================================================\n")

cat("\n>>> Running 06_Meta_Analysis.R...\n")
tryCatch({
  source("06_Meta_Analysis.R", local = TRUE)
  analysis_status[["06_Meta"]] <- "SUCCESS"
  cat("\n>>> Meta-Analysis completed successfully!\n")
}, error = function(e) {
  analysis_status[["06_Meta"]] <- paste("FAILED:", e$message)
  cat("\n>>> Meta-Analysis failed:", e$message, "\n")
})

# =============================================================================
# PART 9: Summary Report
# =============================================================================

end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("################################################################\n")
cat("#                                                              #\n")
cat("#              Analysis Pipeline Completed!                    #\n")
cat("#                                                              #\n")
cat("################################################################\n")
cat("\n")

# Print timing
cat("Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("End time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total duration:", round(as.numeric(duration), 2), "minutes\n")
cat("\n")

# Print analysis status
cat("================================================================\n")
cat("         Analysis Status Summary                                \n")
cat("================================================================\n\n")

n_success <- sum(sapply(analysis_status, function(x) x == "SUCCESS"))
n_failed <- length(analysis_status) - n_success

for (analysis_name in names(analysis_status)) {
  status <- analysis_status[[analysis_name]]
  if (status == "SUCCESS") {
    cat("  [OK]", analysis_name, "\n")
  } else {
    cat("  [FAILED]", analysis_name, "-", status, "\n")
  }
}

cat("\n")
cat("Summary:", n_success, "succeeded,", n_failed, "failed\n")
cat("\n")

# =============================================================================
# PART 10: Output Files Summary
# =============================================================================

cat("================================================================\n")
cat("         Output Files Summary                                   \n")
cat("================================================================\n\n")

cat("Output directory:", OUTPUT_DIR, "\n\n")

cat("=== MAIN RESULTS ===\n\n")

cat("1. Descriptive Statistics:\n")
cat("   - Descriptive_Statistics.xlsx\n")
cat("     * Lifestyle and outcome distributions by cohort\n")
cat("     * All_Outcomes_Summary: Overall + 4 subtypes event counts\n")
cat("   - Table1_by_Cohort.csv\n")
cat("\n")

cat("2. Phi and ICC Analysis:\n")
cat("   - Phi_ICC_Analysis_Results.xlsx\n")
cat("     * Phi coefficients for lifestyle variable correlations\n")
cat("     * ICC for clustering by cohort (all 5 outcomes)\n")
cat("   - Phi_Correlation_Heatmap.pdf/png\n")
cat("\n")

cat("3. Cox Regression (Pooled):\n")
cat("   - Pooled_Cox_Results_Comprehensive.xlsx\n")
cat("     * Primary_Individual: 4 factors mutually adjusted (5 outcomes)\n")
cat("     * Primary_Cumulative: 4-level (0/1/2/3+) with P-trend\n")
cat("     * S1_5Level: 5-level (0/1/2/3/4)\n")
cat("     * S2_HeavyDrink: Heavy drinking definition\n")
cat("     * S3_MICE: MICE imputed data\n")
cat("     * S4_Drop1st: Excluding first follow-up wave\n")
cat("   - DoseResponse_*.png: Dose-response plots\n")
cat("\n")

cat("4. Sankey Diagrams:\n")
cat("   - Sankey_Comprehensive_Results.xlsx\n")
cat("     * Interpretation_Guide: Term definitions\n")
cat("     * All_Transitions_Detail: N, row %, total % for all transitions\n")
cat("     * Transition_Matrices: Wide format matrices with N (%)\n")
cat("     * PPCMM_Summary: PPC-MM incidence by cohort and lifestyle\n")
cat("     * Pooled_Detail/Summary: Pooled data analysis\n")
cat("     * By_Cohort: Each cohort overall summary\n")
cat("   - Sankey_All_Transitions.csv, Sankey_PPCMM_Summary.csv\n")
cat("   - Figures/Sankey/*.pdf/png\n")
cat("     * Sankey_[Cohort]_Overall, Sankey_[Cohort]_Cat[0/1/2/3plus]\n")
cat("     * Sankey_Pooled_Overall, Sankey_Pooled_Cat[0/1/2/3plus]\n")
cat("     * Sankey_Legend.pdf/png (standalone legend file)\n")
cat("\n")

cat("5. PAF Analysis:\n")
cat("   - PAF_Analysis_Comprehensive.xlsx\n")
cat("     * Primary_Individual: Each factor's PAF\n")
cat("     * Primary_Cumulative: Combined PAF by category\n")
cat("     * S2/S3/S4: Sensitivity analyses\n")
cat("\n")

cat("6. Meta-Analysis:\n")
cat("   - Meta_Analysis_Comprehensive.xlsx\n")
cat("     * Study_Characteristics: Sample sizes, demographics by cohort\n")
cat("     * MA_Summary_All: Fixed + Random HR, I2, Q, Tau2\n")
cat("     * Primary/S1_5Level/S2_HeavyDrink/S3_MICE/S4_Drop1st sheets\n")
cat("     * Leave_One_Out: Sensitivity analysis results\n")
cat("     * Eggers_Test: Publication bias assessment\n")
cat("   - Study_Specific_HR_Results.csv: Cohort-specific Cox results\n")
cat("   - Meta_Analysis_Summary.csv\n")
cat("   - Figures/Forest/*.pdf/png: Forest plots (by outcome × level)\n")
cat("   - Figures/Funnel/*.pdf/png: Funnel plots (for k≥3 analyses)\n")
cat("\n")

# =============================================================================
# PART 11: Analysis Framework Summary
# =============================================================================

cat("================================================================\n")
cat("         Analysis Framework (Lancet Standard)                   \n")
cat("================================================================\n\n")

cat("OUTCOMES:\n")
cat("  Primary: Overall PPC-MM (any >=2 conditions)\n")
cat("  Secondary: P1P2, P1C, P2C, P1P2C (4 subtypes)\n")
cat("\n")

cat("EXPOSURES:\n")
cat("  A. Individual Factors (mutually adjusted):\n")
cat("     - Drinking (any drink / heavy drink for sensitivity)\n")
cat("     - Smoking\n")
cat("     - Physical Inactivity\n")
cat("     - Social Isolation\n")
cat("\n")
cat("  B. Cumulative Effect:\n")
cat("     - Main: 0 (ref) / 1 / 2 / 3+\n")
cat("     - Sensitivity: 0 / 1 / 2 / 3 / 4\n")
cat("\n")

cat("SENSITIVITY ANALYSES:\n")
cat("  S1: 5-level categories (0/1/2/3/4)\n")
cat("     - Purpose: Finer dose-response assessment\n")
cat("\n")
cat("  S2: Heavy drinking definition\n")
cat("     - Purpose: Test robustness when using stricter drinking threshold\n")
cat("     - Note: Main analysis uses any drinking, S2 uses heavy drinking\n")
cat("\n")
cat("  S3: MICE imputed data\n")
cat("     - Purpose: Assess impact of missing data on results\n")
cat("\n")
cat("  S4: Drop first follow-up wave\n")
cat("     - Purpose: Address potential reverse causality\n")
cat("\n")

cat("STATISTICAL METHODS:\n")
cat("  - Cox proportional hazards regression (pooled)\n")
cat("  - Meta-analysis: Fixed + Random effects (DerSimonian-Laird)\n")
cat("  - Heterogeneity: Q, I2, Tau2\n")
cat("  - PAF: Miettinen formula\n")
cat("  - Sensitivity: Leave-one-out, Egger's test\n")
cat("\n")

cat("================================================================\n")
cat("         Pipeline Complete                                      \n")
cat("================================================================\n")
cat("\nAll analyses completed.\n")
cat("Check output files in:", OUTPUT_DIR, "\n")
cat("\n")
