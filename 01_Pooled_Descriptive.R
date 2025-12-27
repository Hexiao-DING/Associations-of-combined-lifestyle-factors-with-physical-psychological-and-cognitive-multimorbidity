###############################################################################
# 01_Pooled_Descriptive.R
# ============================================================================
# Project: Lifestyle Factors and PPC-MM Association Study
# Purpose: Generate pooled data and descriptive statistics
# ============================================================================
# This script:
#   1. Reads processed data from each cohort
#   2. Standardizes variable names across cohorts
#   3. Creates pooled dataset combining all cohorts
#   4. Generates baseline characteristics tables (by cohort and overall)
#   5. Summarizes lifestyle and outcome distributions
###############################################################################

# =============================================================================
# PART 1: Environment Setup
# =============================================================================

source("C:/Users/user/Desktop/Scientific Research/UKB Project/Jung Sun Lab/Project1_Lifestyle_PPM_Hexiao,Hongtao/Code/00_Functions_and_Setup.R")

cat("\n")
cat("================================================================\n")
cat("         Pooled Data Generation and Descriptive Statistics      \n")
cat("================================================================\n\n")

# =============================================================================
# PART 2: Read Cohort Data
# =============================================================================

cat("--- Reading cohort analysis data ---\n")

# Use DATA_PATHS from 00_Functions_and_Setup.R (new CSV data structure)
cohort_data <- list()

for (cohort_name in names(DATA_PATHS)) {
  file_path <- DATA_PATHS[[cohort_name]]
  if (!is.null(file_path) && file.exists(file_path)) {
    cohort_data[[cohort_name]] <- read.csv(file_path, stringsAsFactors = FALSE)
    cat("  [OK]", cohort_name, "loaded, N =", nrow(cohort_data[[cohort_name]]), "\n")
  } else {
    cat("  [MISSING]", cohort_name, "file not found:", file_path, "\n")
  }
}

# =============================================================================
# PART 3: Standardize Variable Names and Create Pooled Data
# =============================================================================

cat("\n--- Standardizing variable names and creating pooled data ---\n")

# Standard variable names (extended for new outcome subtypes and heavy drink sensitivity)
standard_vars <- c(
  "cohort", "country",
  "age_baseline", "sex", "edu",
  "unhealthy_drink", "unhealthy_smoke", "unhealthy_pa", "unhealthy_soc", 
  "unhealthy_score", "unhealthy_cat",
  "heavy_drink", "unhealthy_score_heavy",  # For heavy drink sensitivity analysis
  "physical_base", "psych_base", "cog_base", "ppcmm_base",
  "time_ppcmm_months", "event_ppcmm",
  "event_ppcmm_drop1st", "time_ppcmm_drop1st",  # For drop1st sensitivity
  "event_mm_phys_psych", "event_mm_phys_cog", "event_mm_psych_cog", "event_mm_all_three",
  "time_mm_phys_psych", "time_mm_phys_cog", "time_mm_psych_cog", "time_mm_all_three"
)

# Function to standardize cohort data
standardize_cohort <- function(df, cohort_name) {
  baseline_wave <- BASELINE_WAVES[[cohort_name]]
  
  drink_var <- paste0("w", baseline_wave, "_unhealthy_drink")
  smoke_var <- paste0("w", baseline_wave, "_unhealthy_smoke")
  pa_var <- paste0("w", baseline_wave, "_unhealthy_pa")
  soc_var <- paste0("w", baseline_wave, "_unhealthy_soc")
  score_var <- paste0("w", baseline_wave, "_unhealthy_score")
  cat_var <- paste0("w", baseline_wave, "_unhealthy_cat")
  heavy_drink_var <- paste0("w", baseline_wave, "_heavy_drink")
  score_heavy_var <- paste0("w", baseline_wave, "_unhealthy_score_heavy")
  
  df_std <- df %>%
    mutate(
      cohort = cohort_name,
      country = COHORT_COUNTRIES[[cohort_name]]
    )
  
  # Rename lifestyle variables using base R (compatible with all dplyr versions)
  rename_if_exists <- function(data, old_name, new_name) {
    if (old_name %in% names(data)) {
      names(data)[names(data) == old_name] <- new_name
    }
    return(data)
  }
  
  df_std <- rename_if_exists(df_std, drink_var, "unhealthy_drink")
  df_std <- rename_if_exists(df_std, smoke_var, "unhealthy_smoke")
  df_std <- rename_if_exists(df_std, pa_var, "unhealthy_pa")
  df_std <- rename_if_exists(df_std, soc_var, "unhealthy_soc")
  df_std <- rename_if_exists(df_std, score_var, "unhealthy_score")
  df_std <- rename_if_exists(df_std, cat_var, "unhealthy_cat")
  df_std <- rename_if_exists(df_std, heavy_drink_var, "heavy_drink")
  df_std <- rename_if_exists(df_std, score_heavy_var, "unhealthy_score_heavy")
  
  df_std <- df_std %>% select(any_of(standard_vars))
  
  return(df_std)
}

# Standardize each cohort
cohort_data_std <- list()

for (cohort_name in names(cohort_data)) {
  cohort_data_std[[cohort_name]] <- standardize_cohort(
    cohort_data[[cohort_name]],
    cohort_name
  )
}

# Combine all cohorts into pooled dataset
pooled_data <- bind_rows(cohort_data_std)

# Add age_group_4 if not present
if (!"age_group_4" %in% names(pooled_data) && "age_baseline" %in% names(pooled_data)) {
  pooled_data <- pooled_data %>%
    mutate(
      age_group_4 = case_when(
        age_baseline >= 50 & age_baseline <= 59 ~ "50-59",
        age_baseline >= 60 & age_baseline <= 69 ~ "60-69",
        age_baseline >= 70 & age_baseline <= 79 ~ "70-79",
        age_baseline >= 80 ~ "80+",
        TRUE ~ NA_character_
      )
    )
}

cat("\nPooled dataset created successfully\n")
cat("Total sample size:", nrow(pooled_data), "\n")
cat("\nSample size by cohort:\n")
print(table(pooled_data$cohort))

# =============================================================================
# PART 4: Save Pooled Data and Individual Cohort Data
# =============================================================================

cat("\n--- Saving pooled data and individual cohort data ---\n")

# Save pooled data as RDS
saveRDS(pooled_data, file.path(OUTPUT_DIR, "Pooled_main_data.rds"))

# Save pooled data as CSV
save_to_csv(pooled_data, "Pooled_main_data.csv")

# Save each cohort's standardized data as RDS (for downstream analyses)
for (cohort_name in names(cohort_data_std)) {
  # Save standardized data with original wave-specific variables too
  df_original <- cohort_data[[cohort_name]]
  df_std <- cohort_data_std[[cohort_name]]
  
  # Merge standardized names with original data
  df_full <- df_original %>%
    mutate(
      cohort = cohort_name,
      country = COHORT_COUNTRIES[[cohort_name]]
    )
  
  # Add standardized variable names
  bw <- BASELINE_WAVES[[cohort_name]]
  drink_var <- paste0("w", bw, "_unhealthy_drink")
  smoke_var <- paste0("w", bw, "_unhealthy_smoke")
  pa_var <- paste0("w", bw, "_unhealthy_pa")
  soc_var <- paste0("w", bw, "_unhealthy_soc")
  score_var <- paste0("w", bw, "_unhealthy_score")
  cat_var <- paste0("w", bw, "_unhealthy_cat")
  
  if (drink_var %in% names(df_full)) df_full$unhealthy_drink <- df_full[[drink_var]]
  if (smoke_var %in% names(df_full)) df_full$unhealthy_smoke <- df_full[[smoke_var]]
  if (pa_var %in% names(df_full)) df_full$unhealthy_pa <- df_full[[pa_var]]
  if (soc_var %in% names(df_full)) df_full$unhealthy_soc <- df_full[[soc_var]]
  if (score_var %in% names(df_full)) df_full$unhealthy_score <- df_full[[score_var]]
  if (cat_var %in% names(df_full)) df_full$unhealthy_cat <- df_full[[cat_var]]
  
  # Save as RDS
  saveRDS(df_full, file.path(OUTPUT_DIR, paste0(cohort_name, "_main_data.rds")))
  cat("  Saved", cohort_name, "data: N =", nrow(df_full), "\n")
}

cat("Pooled data and individual cohort data saved successfully\n")

# =============================================================================
# PART 5: Baseline Characteristics Tables
# =============================================================================

cat("\n--- Generating baseline characteristics tables ---\n")

# Variables for descriptive statistics
desc_vars <- c("age_baseline", "sex", "edu",
               "unhealthy_drink", "unhealthy_smoke", 
               "unhealthy_pa", "unhealthy_soc", "unhealthy_score",
               "physical_base", "psych_base", "cog_base",
               "event_ppcmm", "time_ppcmm_months")

cat_vars <- c("sex", "edu", "unhealthy_drink", "unhealthy_smoke",
              "unhealthy_pa", "unhealthy_soc", "unhealthy_score",
              "physical_base", "psych_base", "cog_base", "event_ppcmm")

# Overall table
cat("\n=== Overall Baseline Characteristics ===\n")
table1_overall <- CreateTableOne(
  vars = desc_vars[desc_vars %in% names(pooled_data)],
  data = pooled_data,
  factorVars = cat_vars[cat_vars %in% names(pooled_data)]
)
print(table1_overall, showAllLevels = TRUE, quote = FALSE)

# By cohort table
cat("\n=== Baseline Characteristics by Cohort ===\n")
table1_by_cohort <- CreateTableOne(
  vars = desc_vars[desc_vars %in% names(pooled_data)],
  strata = "cohort",
  data = pooled_data,
  factorVars = cat_vars[cat_vars %in% names(pooled_data)]
)
print(table1_by_cohort, showAllLevels = TRUE, quote = FALSE)

# =============================================================================
# PART 6: Lifestyle Distribution Summary
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("         Lifestyle Variable Distribution Summary                \n")
cat("================================================================\n")

# Lifestyle summary by cohort (using 0/1/2/3+ categories)
lifestyle_summary <- pooled_data %>%
  group_by(cohort) %>%
  summarise(
    N = n(),
    pct_drink = mean(unhealthy_drink == 1, na.rm = TRUE) * 100,
    pct_smoke = mean(unhealthy_smoke == 1, na.rm = TRUE) * 100,
    pct_pa = mean(unhealthy_pa == 1, na.rm = TRUE) * 100,
    pct_soc = mean(unhealthy_soc == 1, na.rm = TRUE) * 100,
    mean_score = mean(unhealthy_score, na.rm = TRUE),
    sd_score = sd(unhealthy_score, na.rm = TRUE),
    pct_cat_0 = mean(unhealthy_cat == "0", na.rm = TRUE) * 100,
    pct_cat_1 = mean(unhealthy_cat == "1", na.rm = TRUE) * 100,
    pct_cat_2 = mean(unhealthy_cat == "2", na.rm = TRUE) * 100,
    pct_cat_3plus = mean(unhealthy_cat == "3+", na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat("\nLifestyle distribution by cohort:\n")
print(as.data.frame(lifestyle_summary), digits = 2)

# Overall lifestyle summary
lifestyle_overall <- pooled_data %>%
  summarise(
    cohort = "Overall",
    N = n(),
    pct_drink = mean(unhealthy_drink == 1, na.rm = TRUE) * 100,
    pct_smoke = mean(unhealthy_smoke == 1, na.rm = TRUE) * 100,
    pct_pa = mean(unhealthy_pa == 1, na.rm = TRUE) * 100,
    pct_soc = mean(unhealthy_soc == 1, na.rm = TRUE) * 100,
    mean_score = mean(unhealthy_score, na.rm = TRUE),
    sd_score = sd(unhealthy_score, na.rm = TRUE),
    pct_cat_0 = mean(unhealthy_cat == "0", na.rm = TRUE) * 100,
    pct_cat_1 = mean(unhealthy_cat == "1", na.rm = TRUE) * 100,
    pct_cat_2 = mean(unhealthy_cat == "2", na.rm = TRUE) * 100,
    pct_cat_3plus = mean(unhealthy_cat == "3+", na.rm = TRUE) * 100
  )

cat("\nOverall lifestyle distribution:\n")
print(as.data.frame(lifestyle_overall), digits = 2)

# =============================================================================
# PART 7: Event and Follow-up Summary (All Outcomes)
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("         Event and Follow-up Summary (All Outcomes)             \n")
cat("================================================================\n")

# Define all outcome variables
outcome_vars <- list(
  Overall = list(event = "event_ppcmm", time = "time_ppcmm_months", 
                 label = "Overall PPC-MM", type = "Primary"),
  P1P2 = list(event = "event_mm_phys_psych", time = "time_mm_phys_psych",
              label = "Physical-Psychological", type = "Secondary"),
  P1C = list(event = "event_mm_phys_cog", time = "time_mm_phys_cog",
             label = "Physical-Cognitive", type = "Secondary"),
  P2C = list(event = "event_mm_psych_cog", time = "time_mm_psych_cog",
             label = "Psychological-Cognitive", type = "Secondary"),
  P1P2C = list(event = "event_mm_all_three", time = "time_mm_all_three",
               label = "All Three", type = "Secondary")
)

# Event summary by cohort for Overall (primary)
event_summary <- pooled_data %>%
  group_by(cohort) %>%
  summarise(
    N = n(),
    n_events = sum(event_ppcmm, na.rm = TRUE),
    event_rate_pct = mean(event_ppcmm, na.rm = TRUE) * 100,
    mean_followup_months = mean(time_ppcmm_months, na.rm = TRUE),
    median_followup_months = median(time_ppcmm_months, na.rm = TRUE),
    total_person_years = sum(time_ppcmm_months, na.rm = TRUE) / 12,
    incidence_rate_1000py = n_events / total_person_years * 1000,
    .groups = "drop"
  )

cat("\n--- Overall PPC-MM (Primary Outcome) by Cohort ---\n")
print(as.data.frame(event_summary), digits = 2)

# Overall event summary
event_overall <- pooled_data %>%
  summarise(
    cohort = "Overall",
    N = n(),
    n_events = sum(event_ppcmm, na.rm = TRUE),
    event_rate_pct = mean(event_ppcmm, na.rm = TRUE) * 100,
    mean_followup_months = mean(time_ppcmm_months, na.rm = TRUE),
    median_followup_months = median(time_ppcmm_months, na.rm = TRUE),
    total_person_years = sum(time_ppcmm_months, na.rm = TRUE) / 12,
    incidence_rate_1000py = n_events / total_person_years * 1000
  )

cat("\n--- Overall Summary ---\n")
print(as.data.frame(event_overall), digits = 2)

# Summary for all outcome subtypes
cat("\n--- All Outcome Subtypes Summary ---\n")

all_outcomes_summary <- list()

for (out_name in names(outcome_vars)) {
  out_info <- outcome_vars[[out_name]]
  event_var <- out_info$event
  time_var <- out_info$time
  
  if (event_var %in% names(pooled_data) && time_var %in% names(pooled_data)) {
    outcome_sum <- pooled_data %>%
      summarise(
        Outcome = out_name,
        Outcome_Label = out_info$label,
        Outcome_Type = out_info$type,
        N = n(),
        N_Events = sum(.data[[event_var]], na.rm = TRUE),
        Event_Rate_Pct = round(mean(.data[[event_var]], na.rm = TRUE) * 100, 2),
        Mean_Followup_Months = round(mean(.data[[time_var]], na.rm = TRUE), 1),
        Person_Years = round(sum(.data[[time_var]], na.rm = TRUE) / 12, 0),
        IR_per_1000py = round(N_Events / Person_Years * 1000, 2)
      )
    all_outcomes_summary[[out_name]] <- outcome_sum
  }
}

all_outcomes_df <- bind_rows(all_outcomes_summary)
print(all_outcomes_df)

# Event counts by cohort and outcome
cat("\n--- Event Counts by Cohort and Outcome ---\n")

event_by_cohort_outcome <- pooled_data %>%
  group_by(cohort) %>%
  summarise(
    N = n(),
    Overall = sum(event_ppcmm, na.rm = TRUE),
    P1P2 = sum(event_mm_phys_psych, na.rm = TRUE),
    P1C = sum(event_mm_phys_cog, na.rm = TRUE),
    P2C = sum(event_mm_psych_cog, na.rm = TRUE),
    P1P2C = sum(event_mm_all_three, na.rm = TRUE),
    .groups = "drop"
  )

print(event_by_cohort_outcome)

# =============================================================================
# PART 8: Save Descriptive Results
# =============================================================================

cat("\n--- Saving descriptive statistics results ---\n")

# Combine summaries
descriptive_results <- list(
  "Lifestyle_by_Cohort" = as.data.frame(lifestyle_summary),
  "Lifestyle_Overall" = as.data.frame(lifestyle_overall),
  "Events_by_Cohort" = as.data.frame(event_summary),
  "Events_Overall" = as.data.frame(event_overall),
  "All_Outcomes_Summary" = all_outcomes_df,
  "Events_by_Cohort_Outcome" = event_by_cohort_outcome
)

# Save to Excel
save_to_excel(descriptive_results, "Descriptive_Statistics.xlsx")

# Save tableone results
table1_print <- print(table1_by_cohort, showAllLevels = TRUE, quote = FALSE, printToggle = FALSE)
write.csv(table1_print, file.path(OUTPUT_DIR, "Table1_by_Cohort.csv"))

# =============================================================================
# COMPLETE
# =============================================================================

cat("\n================================================================\n")
cat("         Pooled Data and Descriptive Statistics Complete        \n")
cat("================================================================\n")
cat("Total pooled sample size:", nrow(pooled_data), "\n")
cat("Total events:", sum(pooled_data$event_ppcmm), "\n")
cat("Number of cohorts:", length(unique(pooled_data$cohort)), "\n")
cat("Results saved to:", OUTPUT_DIR, "\n")
