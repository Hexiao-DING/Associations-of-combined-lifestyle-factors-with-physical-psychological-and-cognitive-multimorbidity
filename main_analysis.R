###############################################################################
# main_analysis.R
# ============================================================================
# Project: Lifestyle Factors and PPC-MM Association Study
# Purpose: Master script to run all analysis (Lancet Standard)
# ============================================================================
#
# SCRIPT STRUCTURE (Numbered 00-08):
# ==================================
#   00_Functions_and_Setup.R  - Environment setup, paths, covariate configuration
#   01_Pooled_Descriptive.R   - Pooled data generation and descriptive statistics
#   02_MICE_Imputation.R      - 2-Level MICE imputation (optional, run separately)
#   03_Phi_ICC_Analysis.R     - Phi coefficient and ICC analysis
#   04_Pooled_Cox_Analysis.R  - Cox proportional hazards regression (pooled)
#   05_Sankey_Diagram.R       - Health state transition diagrams
#   06_PAF_Analysis.R         - Population Attributable Fraction
#   07_Meta_Analysis.R        - Meta-analysis across cohorts
#   08_Methods_Parameters.R   - Export all analysis parameters
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
# PART 2: Load Functions and Setup (00_Functions_and_Setup.R)
# =============================================================================

cat("\n================================================================\n")
cat("  [00] Loading Functions and Environment Setup                   \n")
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
  }
  
}, error = function(e) {
  analysis_status[["00_Setup"]] <- paste("FAILED:", e$message)
  cat("\n>>> Setup failed:", e$message, "\n")
  stop("Cannot proceed without setup. Please fix the error above.")
})

# =============================================================================
# PART 3: Pooled Data Generation (01_Pooled_Descriptive.R)
# =============================================================================

cat("\n================================================================\n")
cat("  [01] Pooled Data Generation and Descriptive Statistics         \n")
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
# PART 4: MICE Imputation (02_MICE_Imputation.R) - OPTIONAL
# =============================================================================

# Note: MICE imputation is computationally intensive (may take 30+ minutes)
# Uncomment below to run as part of pipeline, or run 02_MICE_Imputation.R separately

# cat("\n================================================================\n")
# cat("  [02] 2-Level MICE Imputation                                   \n")
# cat("================================================================\n")
# 
# cat("\n>>> Running 02_MICE_Imputation.R...\n")
# cat(">>> Note: This may take 30+ minutes...\n")
# tryCatch({
#   source("02_MICE_Imputation.R", local = TRUE)
#   analysis_status[["02_MICE"]] <- "SUCCESS"
#   cat("\n>>> MICE Imputation completed successfully!\n")
# }, error = function(e) {
#   analysis_status[["02_MICE"]] <- paste("FAILED:", e$message)
#   cat("\n>>> MICE Imputation failed:", e$message, "\n")
# })

# =============================================================================
# PART 5: Phi Coefficient and ICC Analysis (03_Phi_ICC_Analysis.R)
# =============================================================================

cat("\n================================================================\n")
cat("  [03] Phi Coefficient and ICC Analysis                          \n")
cat("================================================================\n")

cat("\n>>> Running 03_Phi_ICC_Analysis.R...\n")
tryCatch({
  source("03_Phi_ICC_Analysis.R", local = TRUE)
  analysis_status[["03_Phi_ICC"]] <- "SUCCESS"
  cat("\n>>> Phi and ICC Analysis completed successfully!\n")
}, error = function(e) {
  analysis_status[["03_Phi_ICC"]] <- paste("FAILED:", e$message)
  cat("\n>>> Phi and ICC Analysis failed:", e$message, "\n")
})

# =============================================================================
# PART 6: Pooled Cox Regression Analysis (04_Pooled_Cox_Analysis.R)
# =============================================================================

cat("\n================================================================\n")
cat("  [04] Pooled Cox Proportional Hazards Regression                \n")
cat("================================================================\n")

cat("\n>>> Running 04_Pooled_Cox_Analysis.R...\n")
tryCatch({
  source("04_Pooled_Cox_Analysis.R", local = TRUE)
  analysis_status[["04_Cox"]] <- "SUCCESS"
  cat("\n>>> Pooled Cox Analysis completed successfully!\n")
}, error = function(e) {
  analysis_status[["04_Cox"]] <- paste("FAILED:", e$message)
  cat("\n>>> Pooled Cox Analysis failed:", e$message, "\n")
})

# =============================================================================
# PART 7: Sankey Diagrams (05_Sankey_Diagram.R)
# =============================================================================

cat("\n================================================================\n")
cat("  [05] Health State Transition (Sankey) Diagrams                 \n")
cat("================================================================\n")

cat("\n>>> Running 05_Sankey_Diagram.R...\n")
tryCatch({
  source("05_Sankey_Diagram.R", local = TRUE)
  analysis_status[["05_Sankey"]] <- "SUCCESS"
  cat("\n>>> Sankey Diagram Generation completed successfully!\n")
}, error = function(e) {
  analysis_status[["05_Sankey"]] <- paste("FAILED:", e$message)
  cat("\n>>> Sankey Diagram Generation failed:", e$message, "\n")
})

# =============================================================================
# PART 8: PAF Analysis (06_PAF_Analysis.R)
# =============================================================================

cat("\n================================================================\n")
cat("  [06] Population Attributable Fraction (PAF) Analysis           \n")
cat("================================================================\n")

cat("\n>>> Running 06_PAF_Analysis.R...\n")
tryCatch({
  source("06_PAF_Analysis.R", local = TRUE)
  analysis_status[["06_PAF"]] <- "SUCCESS"
  cat("\n>>> PAF Analysis completed successfully!\n")
}, error = function(e) {
  analysis_status[["06_PAF"]] <- paste("FAILED:", e$message)
  cat("\n>>> PAF Analysis failed:", e$message, "\n")
})

# =============================================================================
# PART 9: Meta-Analysis (07_Meta_Analysis.R)
# =============================================================================

cat("\n================================================================\n")
cat("  [07] Meta-Analysis Across Cohorts                              \n")
cat("================================================================\n")

cat("\n>>> Running 07_Meta_Analysis.R...\n")
tryCatch({
  source("07_Meta_Analysis.R", local = TRUE)
  analysis_status[["07_Meta"]] <- "SUCCESS"
  cat("\n>>> Meta-Analysis completed successfully!\n")
}, error = function(e) {
  analysis_status[["07_Meta"]] <- paste("FAILED:", e$message)
  cat("\n>>> Meta-Analysis failed:", e$message, "\n")
})

# =============================================================================
# PART 10: Methods Parameters Export (08_Methods_Parameters.R)
# =============================================================================

cat("\n================================================================\n")
cat("  [08] Methods Parameters Export                                 \n")
cat("================================================================\n")

cat("\n>>> Running 08_Methods_Parameters.R...\n")
tryCatch({
  source("08_Methods_Parameters.R", local = TRUE)
  analysis_status[["08_Methods"]] <- "SUCCESS"
  cat("\n>>> Methods Parameters Export completed successfully!\n")
}, error = function(e) {
  analysis_status[["08_Methods"]] <- paste("FAILED:", e$message)
  cat("\n>>> Methods Parameters Export failed:", e$message, "\n")
})

# =============================================================================
# PART 11: Summary Report
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
# PART 12: Output Files Summary
# =============================================================================

cat("================================================================\n")
cat("         Output Files Summary                                   \n")
cat("================================================================\n\n")

cat("Output directory:", OUTPUT_DIR, "\n\n")

cat("=== MAIN RESULTS (Excel Workbooks) ===\n\n")

cat("1. Descriptive Statistics (01_Pooled_Descriptive.R):\n")
cat("   [Excel] Descriptive_Statistics.xlsx\n")
cat("   [CSV]   Table1_by_Cohort.csv\n")
cat("   [CSV]   Pooled_main_data.csv\n")
cat("   [RDS]   Pooled_main_data.rds\n")
cat("\n")

cat("2. Phi and ICC Analysis (03_Phi_ICC_Analysis.R):\n")
cat("   [Excel] Phi_ICC_Analysis_Results.xlsx\n")
cat("   [CSV]   Phi_Coefficients_Complete.csv\n")
cat("   [CSV]   ICC_Results_Complete.csv\n")
cat("   [Plot]  Phi_Correlation_Heatmap.pdf/png\n")
cat("\n")

cat("3. Cox Regression (04_Pooled_Cox_Analysis.R):\n")
cat("   [Excel] Pooled_Cox_Results_Comprehensive.xlsx\n")
cat("   [CSV]   Pooled_Cox_All_Results.csv (All results combined)\n")
cat("   [CSV]   Cox_Primary_Individual_Factors.csv\n")
cat("   [CSV]   Cox_Primary_Cumulative_4Level.csv\n")
cat("   [CSV]   Cox_S1_5Level_Categories.csv\n")
cat("   [CSV]   Cox_S2a_HeavyDrink_Individual.csv\n")
cat("   [CSV]   Cox_S2b_HeavyDrink_Cumulative.csv\n")
cat("   [CSV]   Cox_S3_MICE_Individual.csv (if MICE data available)\n")
cat("   [CSV]   Cox_S3_MICE_Cumulative.csv (if MICE data available)\n")
cat("   [CSV]   Cox_S4_Drop1st.csv\n")
cat("   [CSV]   Cox_Summary_HR_CI.csv (Publication format)\n")
cat("   [Plot]  DoseResponse_*.png/pdf\n")
cat("   [RDS]   Pooled_cox_results.rds\n")
cat("\n")

cat("4. Sankey Diagrams (05_Sankey_Diagram.R):\n")
cat("   [Excel] Sankey_Comprehensive_Results.xlsx\n")
cat("   [CSV]   Sankey_All_Transitions.csv\n")
cat("   [CSV]   Sankey_PPCMM_Summary.csv\n")
cat("   [Plot]  Figures/Sankey/Sankey_*.pdf/png\n")
cat("   [Plot]  Figures/Sankey/Sankey_Legend.pdf/png\n")
cat("\n")

cat("5. PAF Analysis (06_PAF_Analysis.R):\n")
cat("   [Excel] PAF_Analysis_Comprehensive.xlsx\n")
cat("   [CSV]   PAF_Analysis_All_Results.csv (All results combined)\n")
cat("   [CSV]   PAF_Primary_Individual.csv\n")
cat("   [CSV]   PAF_Primary_Cumulative.csv\n")
cat("   [CSV]   PAF_S2_HeavyDrink_Individual.csv\n")
cat("   [CSV]   PAF_S2_HeavyDrink_Cumulative.csv\n")
cat("   [CSV]   PAF_S4_Drop1st.csv\n")
cat("   [CSV]   PAF_Summary_Overall.csv (Publication format)\n")
cat("   [RDS]   PAF_results.rds\n")
cat("\n")

cat("6. Meta-Analysis (07_Meta_Analysis.R):\n")
cat("   [Excel] Meta_Analysis_Comprehensive.xlsx\n")
cat("   [CSV]   Meta_Analysis_Summary.csv (All results combined)\n")
cat("   [CSV]   Meta_Study_Characteristics.csv\n")
cat("   [CSV]   Meta_Study_Specific_HR.csv (Cohort-specific HRs)\n")
cat("   [CSV]   Meta_Primary_Summary.csv\n")
cat("   [CSV]   Meta_S1_5Level_Summary.csv\n")
cat("   [CSV]   Meta_S2_HeavyDrink_Summary.csv\n")
cat("   [CSV]   Meta_LeaveOneOut.csv\n")
cat("   [CSV]   Meta_EggersTest.csv\n")
cat("   [CSV]   Meta_Summary_Publication.csv (Publication format)\n")
cat("   [Plot]  Figures/Forest/Forest_*.pdf/png\n")
cat("   [Plot]  Figures/Funnel/Funnel_*.pdf/png\n")
cat("   [RDS]   Meta_analysis_objects.rds\n")
cat("\n")

cat("7. Methods Parameters (08_Methods_Parameters.R):\n")
cat("   [Excel] Methods_Parameters.xlsx (11 sheets)\n")
cat("\n")

cat("8. MICE Imputation (02_MICE_Imputation.R, if run):\n")
cat("   [Excel] Missing_Data_Summary.xlsx\n")
cat("   [RDS]   Pooled_mice_imputed.rds\n")
cat("   [Plot]  MICE_Convergence_*.png\n")
cat("   [Plot]  MICE_Density_*.png\n")
cat("\n")

cat("================================================================\n")
cat("         Pipeline Complete                                      \n")
cat("================================================================\n")
cat("\nAll analyses completed.\n")
cat("Check output files in:", OUTPUT_DIR, "\n")
cat("\n")
