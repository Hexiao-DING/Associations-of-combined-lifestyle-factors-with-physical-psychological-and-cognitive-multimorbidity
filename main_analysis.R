###############################################################################
#                                                                             #
#                           main_analysis.R                                   #
#                                                                             #
#  Project: Lifestyle Factors and PPC-MM Association Study                    #
#  Purpose: Master script to run all analyses (Lancet Standard)               #
#                                                                             #
###############################################################################
#
#  SCRIPT ORDER:
#  =============
#    00_Functions_and_Setup.R   - Environment, paths, configuration
#    01_Pooled_Descriptive.R    - Data pooling, descriptive statistics
#    02_MICE_Imputation.R       - 2-Level MICE imputation [Optional]
#    03_Phi_ICC_Analysis.R      - Phi coefficient and ICC analysis
#    04_Pooled_Cox_Analysis.R   - Cox proportional hazards regression
#    05_Subgroup_Analysis.R     - Age/Sex stratification, Interaction
#    06_Sankey_Diagram.R        - Health state transition diagrams
#    07_PAF_Analysis.R          - Population Attributable Fraction
#    08_Meta_Analysis.R         - Meta-analysis across cohorts
#    09_Methods_Parameters.R    - Export analysis parameters
#
###############################################################################

# =============================================================================
# INITIALIZATION
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

CODE_DIR <- "C:/Users/user/Desktop/Scientific Research/UKB Project/Jung Sun Lab/Project1_Lifestyle_PPM_Hexiao,Hongtao/Code"
setwd(CODE_DIR)

analysis_status <- list()

# Helper function
run_script <- function(script_name, description) {
  cat("\n================================================================\n")
  cat(" ", description, "\n")
  cat("================================================================\n\n")
  
  result <- tryCatch({
    source(script_name, local = FALSE)
    list(success = TRUE, message = "SUCCESS")
  }, error = function(e) {
    list(success = FALSE, message = e$message)
  })
  
  if (result$success) {
    cat("\n>>> [OK]", script_name, "\n")
  } else {
    cat("\n>>> [FAILED]", script_name, ":", result$message, "\n")
  }
  
  return(result)
}

# =============================================================================
# STEP 00: Environment Setup
# =============================================================================

cat("\n================================================================\n")
cat("  [00] Environment Setup                                         \n")
cat("================================================================\n\n")

result <- tryCatch({
  source("00_Functions_and_Setup.R", local = FALSE)
  analysis_status[["00_Setup"]] <- "SUCCESS"
  cat("\n>>> [OK] Setup completed\n")
}, error = function(e) {
  analysis_status[["00_Setup"]] <- paste("FAILED:", e$message)
  stop("Setup failed: ", e$message)
})

# =============================================================================
# STEP 01: Pooled Data Generation
# =============================================================================

result <- run_script("01_Pooled_Descriptive.R", "[01] Pooled Data & Descriptive")
analysis_status[["01_Descriptive"]] <- ifelse(result$success, "SUCCESS", result$message)

# =============================================================================
# STEP 02: MICE Imputation (Optional)
# =============================================================================

# Uncomment to run MICE (takes ~30-45 minutes)
# result <- run_script("02_MICE_Imputation.R", "[02] MICE Imputation")
# analysis_status[["02_MICE"]] <- ifelse(result$success, "SUCCESS", result$message)

# =============================================================================
# STEP 03: Phi and ICC Analysis
# =============================================================================

result <- run_script("03_Phi_ICC_Analysis.R", "[03] Phi & ICC Analysis")
analysis_status[["03_Phi_ICC"]] <- ifelse(result$success, "SUCCESS", result$message)

# =============================================================================
# STEP 04: Cox Regression Analysis
# =============================================================================

result <- run_script("04_Pooled_Cox_Analysis.R", "[04] Cox Regression")
analysis_status[["04_Cox"]] <- ifelse(result$success, "SUCCESS", result$message)

# =============================================================================
# STEP 05: Subgroup Analysis
# =============================================================================

result <- run_script("05_Subgroup_Analysis.R", "[05] Subgroup Analysis")
analysis_status[["05_Subgroup"]] <- ifelse(result$success, "SUCCESS", result$message)

# =============================================================================
# STEP 06: Sankey Diagrams
# =============================================================================

result <- run_script("06_Sankey_Diagram.R", "[06] Sankey Diagrams")
analysis_status[["06_Sankey"]] <- ifelse(result$success, "SUCCESS", result$message)

# =============================================================================
# STEP 07: PAF Analysis
# =============================================================================

result <- run_script("07_PAF_Analysis.R", "[07] PAF Analysis")
analysis_status[["07_PAF"]] <- ifelse(result$success, "SUCCESS", result$message)

# =============================================================================
# STEP 08: Meta-Analysis
# =============================================================================

result <- run_script("08_Meta_Analysis.R", "[08] Meta-Analysis")
analysis_status[["08_Meta"]] <- ifelse(result$success, "SUCCESS", result$message)

# =============================================================================
# STEP 09: Methods Parameters
# =============================================================================

result <- run_script("09_Methods_Parameters.R", "[09] Methods Parameters")
analysis_status[["09_Methods"]] <- ifelse(result$success, "SUCCESS", result$message)

# =============================================================================
# SUMMARY
# =============================================================================

end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("################################################################\n")
cat("#              ANALYSIS PIPELINE COMPLETED                     #\n")
cat("################################################################\n")
cat("\n")

cat("Duration:", round(as.numeric(duration), 2), "minutes\n\n")

cat("================================================================\n")
cat("  STATUS SUMMARY\n")
cat("================================================================\n\n")

n_success <- sum(sapply(analysis_status, function(x) x == "SUCCESS"))
n_failed <- length(analysis_status) - n_success

for (name in names(analysis_status)) {
  status <- analysis_status[[name]]
  if (status == "SUCCESS") {
    cat("  [OK]  ", name, "\n")
  } else {
    cat("  [FAIL]", name, ":", status, "\n")
  }
}

cat("\n  Total:", n_success, "succeeded,", n_failed, "failed\n")

# Output verification
cat("\n================================================================\n")
cat("  OUTPUT VERIFICATION\n")
cat("================================================================\n\n")

key_files <- c(
  # Excel
  "Descriptive_Statistics.xlsx",
  "Pooled_Cox_Results_Comprehensive.xlsx",
  "Subgroup_Analysis_Complete.xlsx",
  "PAF_Analysis_Comprehensive.xlsx",
  "Meta_Analysis_Comprehensive.xlsx",
  "Methods_Parameters.xlsx",
  # Subgroup CSVs
  "Cox_Age_Stratified.csv",
  "Cox_Sex_Stratified.csv",
  "Cox_Interaction_Analysis.csv",
  "PAF_Age_Stratified_Individual.csv",
  "PAF_Sex_Stratified_Individual.csv"
)

for (f in key_files) {
  fp <- file.path(OUTPUT_DIR, f)
  if (file.exists(fp)) {
    cat("  [OK]", f, "\n")
  } else {
    cat("  [--]", f, "\n")
  }
}

cat("\n================================================================\n")
cat("  COMPLETE\n")
cat("================================================================\n")
cat("\nOutput directory:", OUTPUT_DIR, "\n\n")
