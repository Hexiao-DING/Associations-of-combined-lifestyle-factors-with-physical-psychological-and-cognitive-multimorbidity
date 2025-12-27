###############################################################################
# 05_Sankey_Diagram.R
# ============================================================================
# Project: Lifestyle Factors and PPC-MM Association Study
# Purpose: Create Sankey diagrams showing health state transitions
# ============================================================================
#
# ANALYSIS OVERVIEW:
# ------------------
# BASELINE STATES (4 categories):
#   - Healthy: No physical, psychological, or cognitive conditions
#   - P1: Physical condition only
#   - P2: Psychological condition only  
#   - C: Cognitive condition only
#
# FOLLOW-UP STATES (8 categories):
#   - Healthy, P1, P2, C: Single/no conditions
#   - P1P2, P1C, P2C: Two conditions (PPC-MM subtypes)
#   - P1P2C: All three conditions (PPC-MM subtype)
#
# OUTPUT:
#   - Sankey diagrams with N and percentages for each state
#   - Detailed transition tables in Excel format
#   - Summary statistics for public health interpretation
#
###############################################################################

# =============================================================================
# PART 1: ENVIRONMENT SETUP
# =============================================================================

cat("\n")
cat("###############################################################################\n")
cat("#              SANKEY DIAGRAM: HEALTH STATE TRANSITIONS                       #\n")
cat("###############################################################################\n")
cat("\n")

source("C:/Users/user/Desktop/Scientific Research/UKB Project/Jung Sun Lab/Project1_Lifestyle_PPM_Hexiao,Hongtao/Code/00_Functions_and_Setup.R")

# Install and load required packages
required_pkgs <- c("ggalluvial", "ggplot2", "scales", "ggrepel", "patchwork")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cat("[INFO] All required packages loaded.\n\n")

# =============================================================================
# PART 2: CONFIGURATION
# =============================================================================

# State definitions
BASELINE_STATES <- c("Healthy", "P1", "P2", "C")
FOLLOWUP_STATES <- c("Healthy", "P1", "P2", "C", "P1P2", "P1C", "P2C", "P1P2C")

# Color palette
STATE_COLORS <- c(
  "Healthy" = "#2ECC71",   # Green
  "P1"      = "#E74C3C",   # Red - Physical
  "P2"      = "#3498DB",   # Blue - Psychological
  "C"       = "#F39C12",   # Orange - Cognitive
  "P1P2"    = "#9B59B6",   # Purple
  "P1C"     = "#E67E22",   # Dark Orange
  "P2C"     = "#1ABC9C",   # Turquoise
  "P1P2C"   = "#C0392B"    # Dark Red
)

# Cohort baseline waves
COHORT_WAVES <- list(
  CHARLS = list(baseline = 1, followup = c(2, 3, 4)),
  ELSA   = list(baseline = 7, followup = c(8, 9)),
  HRS    = list(baseline = 10, followup = c(11, 12, 13, 14)),  # Updated: HRS uses w10 as baseline
  SHARE  = list(baseline = 4, followup = c(5, 6, 7, 8)),
  MHAS   = list(baseline = 3, followup = c(4, 5))
)

# Create output directories
fig_dir <- file.path(OUTPUT_DIR, "Figures", "Sankey")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
cat("[INFO] Output directory:", fig_dir, "\n\n")

# =============================================================================
# PART 3: HELPER FUNCTIONS
# =============================================================================

cat("=============================================================================\n")
cat("DEFINING HELPER FUNCTIONS\n")
cat("=============================================================================\n\n")

#' Get baseline state based on indicator variables
get_baseline_state <- function(healthy, physical_only, psych_only, cog_only) {
  dplyr::case_when(
    healthy == 1 ~ "Healthy",
    physical_only == 1 ~ "P1",
    psych_only == 1 ~ "P2",
    cog_only == 1 ~ "C",
    TRUE ~ NA_character_
  )
}

#' Get follow-up state based on indicator variables
get_followup_state <- function(healthy, physical_only, psych_only, cog_only,
                                mm_phys_psych, mm_phys_cog, mm_psych_cog, mm_all_three) {
  dplyr::case_when(
    mm_all_three == 1 ~ "P1P2C",
    mm_phys_psych == 1 ~ "P1P2",
    mm_phys_cog == 1 ~ "P1C",
    mm_psych_cog == 1 ~ "P2C",
    physical_only == 1 ~ "P1",
    psych_only == 1 ~ "P2",
    cog_only == 1 ~ "C",
    healthy == 1 ~ "Healthy",
    TRUE ~ NA_character_
  )
}

#' Find last available follow-up state for each participant
find_last_followup_state <- function(df, cohort) {
  waves_info <- COHORT_WAVES[[cohort]]
  if (is.null(waves_info)) return(df)
  
  followup_waves <- waves_info$followup
  df$followup_state <- NA_character_
  df$followup_wave <- NA_integer_
  
  # Check from last wave to first (prefer later observations)
  for (w in rev(followup_waves)) {
    # Variable names for this wave
    h_var <- paste0("healthy_w", w)
    p1_var <- paste0("physical_only_w", w)
    p2_var <- paste0("psych_only_w", w)
    c_var <- paste0("cog_only_w", w)
    pp_var <- paste0("mm_phys_psych_w", w)
    pc_var <- paste0("mm_phys_cog_w", w)
    p2c_var <- paste0("mm_psych_cog_w", w)
    all3_var <- paste0("mm_all_three_w", w)
    
    # Check if all variables exist
    all_vars <- c(h_var, p1_var, p2_var, c_var, pp_var, pc_var, p2c_var, all3_var)
    if (!all(all_vars %in% names(df))) next
    
    # Process each row
    for (i in which(is.na(df$followup_state))) {
      vals <- c(df[[h_var]][i], df[[p1_var]][i], df[[p2_var]][i], df[[c_var]][i],
                df[[pp_var]][i], df[[pc_var]][i], df[[p2c_var]][i], df[[all3_var]][i])
      
      if (!all(is.na(vals))) {
        state <- get_followup_state(
          ifelse(is.na(vals[1]), 0, vals[1]),
          ifelse(is.na(vals[2]), 0, vals[2]),
          ifelse(is.na(vals[3]), 0, vals[3]),
          ifelse(is.na(vals[4]), 0, vals[4]),
          ifelse(is.na(vals[5]), 0, vals[5]),
          ifelse(is.na(vals[6]), 0, vals[6]),
          ifelse(is.na(vals[7]), 0, vals[7]),
          ifelse(is.na(vals[8]), 0, vals[8])
        )
        if (!is.na(state)) {
          df$followup_state[i] <- state
          df$followup_wave[i] <- w
        }
      }
    }
  }
  return(df)
}

#' Prepare Sankey data for a cohort
prepare_sankey_data <- function(df, cohort) {
  waves_info <- COHORT_WAVES[[cohort]]
  if (is.null(waves_info)) {
    cat("    [ERROR] No wave info for cohort:", cohort, "\n")
    return(NULL)
  }
  
  bw <- waves_info$baseline
  
  # Variable names
  h_base <- paste0("healthy_w", bw)
  p1_base <- paste0("physical_only_w", bw)
  p2_base <- paste0("psych_only_w", bw)
  c_base <- paste0("cog_only_w", bw)
  score_base <- paste0("w", bw, "_unhealthy_score")
  
  # Check required variables
  base_vars <- c(h_base, p1_base, p2_base, c_base)
  missing_vars <- base_vars[!base_vars %in% names(df)]
  if (length(missing_vars) > 0) {
    cat("    [ERROR] Missing vars:", paste(missing_vars, collapse = ", "), "\n")
    return(NULL)
  }
  
  # Add baseline state
  df$baseline_state <- get_baseline_state(
    df[[h_base]], df[[p1_base]], df[[p2_base]], df[[c_base]]
  )
  
  # Find follow-up state
  df <- find_last_followup_state(df, cohort)
  
  # Add lifestyle category
  if (score_base %in% names(df)) {
    df$lifestyle_cat <- dplyr::case_when(
      df[[score_base]] == 0 ~ "0",
      df[[score_base]] == 1 ~ "1",
      df[[score_base]] == 2 ~ "2",
      df[[score_base]] >= 3 ~ "3+",
      TRUE ~ NA_character_
    )
    df$unhealthy_score <- df[[score_base]]
  }
  
  # Filter valid transitions
  df_valid <- df[!is.na(df$baseline_state) & !is.na(df$followup_state), ]
  
  return(df_valid)
}

#' Create detailed transition table with N and percentages
create_transition_table <- function(df, group_name = "Overall") {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  total_n <- nrow(df)
  
  # Count transitions
  trans <- df %>%
    dplyr::group_by(baseline_state, followup_state) %>%
    dplyr::summarise(N = dplyr::n(), .groups = "drop")
  
  # Add percentages
  trans <- trans %>%
    dplyr::group_by(baseline_state) %>%
    dplyr::mutate(
      N_baseline = sum(N),
      Pct_row = round(N / N_baseline * 100, 1)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Pct_total = round(N / total_n * 100, 2),
      Group = group_name,
      Total_N = total_n
    )
  
  return(trans)
}

#' Create wide-format transition matrix
create_transition_matrix <- function(trans_table) {
  if (is.null(trans_table) || nrow(trans_table) == 0) return(NULL)
  
  # Create cell values with N (%)
  trans_table$Cell <- paste0(trans_table$N, " (", trans_table$Pct_row, "%)")
  
  # Pivot to wide format
  matrix_wide <- trans_table %>%
    dplyr::select(Group, baseline_state, followup_state, Cell) %>%
    tidyr::pivot_wider(
      names_from = followup_state, 
      values_from = Cell, 
      values_fill = "0 (0%)"
    )
  
  return(matrix_wide)
}

#' Create standalone legend plot - compact layout for separate file
#' 2 columns x 4 rows
create_standalone_legend <- function() {
  # 2 columns x 4 rows layout (left: single states, right: MM states)
  legend_data <- data.frame(
    abbrev = c("Healthy", "P1", "P2", "C", "P1P2", "P1C", "P2C", "P1P2C"),
    full_name = c(
      "No condition",
      "Physical condition",
      "Psychological condition",
      "Cognitive impairment",
      "Physical + Psychological",
      "Physical + Cognitive",
      "Psychological + Cognitive",
      "All three conditions"
    ),
    color = c("#3D9970", "#0074D9", "#7FDBFF", "#FFDC00", 
              "#FF851B", "#FF4136", "#F4A582", "#B10DC9"),
    x = c(0, 0, 0, 0, 1.8, 1.8, 1.8, 1.8),  # 2 columns
    y = c(4, 3, 2, 1, 4, 3, 2, 1),          # 4 rows
    stringsAsFactors = FALSE
  )
  
  # Create label text
  legend_data$label <- paste0(legend_data$abbrev, ": ", legend_data$full_name)
  
  # Create compact legend
  p_legend <- ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(color = abbrev), size = 4, shape = 15) +
    ggplot2::geom_text(ggplot2::aes(label = label), hjust = 0, nudge_x = 0.08, 
                       size = 3.5, color = "grey20") +
    ggplot2::scale_color_manual(
      values = setNames(legend_data$color, legend_data$abbrev),
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(limits = c(-0.2, 4.2)) +
    ggplot2::scale_y_continuous(limits = c(0.5, 4.5)) +
    ggplot2::labs(title = "Legend: Health State Definitions") +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5,
                                         margin = ggplot2::margin(b = 10)),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(15, 20, 15, 20)
    )
  
  return(p_legend)
}

#' Create Sankey plot with labels adjacent to strata
#' Format: "State (n=xxx, xx.x%)" directly next to each color block
#' Clean design with legend below showing abbreviations and full names
create_sankey_plot <- function(trans_table, title, subtitle = "") {
  if (is.null(trans_table) || nrow(trans_table) == 0) return(NULL)
  
  # Prepare data
  plot_data <- trans_table %>%
    dplyr::mutate(
      baseline_state = factor(baseline_state, levels = BASELINE_STATES),
      followup_state = factor(followup_state, levels = FOLLOWUP_STATES)
    ) %>%
    dplyr::filter(!is.na(baseline_state) & !is.na(followup_state))
  
  if (nrow(plot_data) == 0) return(NULL)
  
  # Calculate total N
  total_n <- sum(plot_data$N)
  
  # Calculate baseline summaries with labels "State (n=xxx, xx.x%)"
  baseline_summary <- plot_data %>%
    dplyr::group_by(baseline_state) %>%
    dplyr::summarise(n = sum(N), .groups = "drop") %>%
    dplyr::mutate(
      pct = round(n / total_n * 100, 1),
      label = paste0(baseline_state, " (n=", format(n, big.mark = ","), ", ", pct, "%)")
    ) %>%
    dplyr::arrange(match(baseline_state, BASELINE_STATES))
  
  # Calculate follow-up summaries with labels "State (n=xxx, xx.x%)"
  followup_summary <- plot_data %>%
    dplyr::group_by(followup_state) %>%
    dplyr::summarise(n = sum(N), .groups = "drop") %>%
    dplyr::mutate(
      pct = round(n / total_n * 100, 1),
      label = paste0(followup_state, " (n=", format(n, big.mark = ","), ", ", pct, "%)")
    ) %>%
    dplyr::arrange(match(followup_state, FOLLOWUP_STATES))
  
  # Create label mappings
  baseline_label_map <- setNames(baseline_summary$label, as.character(baseline_summary$baseline_state))
  followup_label_map <- setNames(followup_summary$label, as.character(followup_summary$followup_state))
  
  # Replace state names with labels in data
  plot_data <- plot_data %>%
    dplyr::mutate(
      baseline_label = factor(baseline_label_map[as.character(baseline_state)], 
                              levels = baseline_summary$label),
      followup_label = factor(followup_label_map[as.character(followup_state)], 
                              levels = followup_summary$label)
    )
  
  # Define state colors
  state_colors <- c(
    "Healthy" = "#3D9970",  # Green
    "P1" = "#0074D9",       # Blue  
    "P2" = "#7FDBFF",       # Light blue
    "C" = "#FFDC00",        # Yellow
    "P1P2" = "#FF851B",     # Orange
    "P1C" = "#FF4136",      # Red
    "P2C" = "#F4A582",      # Salmon
    "P1P2C" = "#B10DC9"     # Purple
  )
  
  # Create color vector for labeled factors
  color_vec <- c(
    # Baseline colors
    setNames(
      state_colors[as.character(baseline_summary$baseline_state)],
      baseline_summary$label
    ),
    # Followup colors
    setNames(
      state_colors[as.character(followup_summary$followup_state)],
      followup_summary$label
    )
  )
  
  # Create plot - labels are embedded in the factor names
  # Wider stratum width so text fits inside the color blocks
  stratum_width <- 1/4  # Wider blocks for text to fit inside
  
  p <- ggplot2::ggplot(plot_data, 
                       ggplot2::aes(y = N, axis1 = baseline_label, axis2 = followup_label)) +
    ggalluvial::geom_alluvium(
      ggplot2::aes(fill = followup_label),
      width = stratum_width,
      curve_type = "sigmoid",
      alpha = 0.45
    ) +
    ggalluvial::geom_stratum(
      ggplot2::aes(fill = ggplot2::after_stat(stratum)),
      width = stratum_width,
      color = "white",
      linewidth = 0.5
    ) +
    ggplot2::scale_fill_manual(values = color_vec, guide = "none") +
    ggplot2::scale_x_discrete(expand = c(0.35, 0.05)) +
    ggplot2::scale_y_continuous(expand = c(0.02, 0)) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0, 
                                         margin = ggplot2::margin(b = 10)),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(15, 15, 15, 15)
    )
  
  # Add text labels inside the strata
  p <- p + ggplot2::geom_text(
    stat = ggalluvial::StatStratum,
    ggplot2::aes(label = ggplot2::after_stat(stratum)),
    size = 2.8,
    color = "black",
    fontface = "plain"
  )
  
  # Return main plot only (no legend - legend is saved separately)
  return(p)
}

#' Safe plot saving function
save_plot <- function(plot, filepath_base, width = 12, height = 8) {
  if (is.null(plot)) {
    cat("      [SKIP] Plot is NULL\n")
    return(FALSE)
  }
  
  success <- FALSE
  
  # Save PDF
  pdf_path <- paste0(filepath_base, ".pdf")
  tryCatch({
    ggplot2::ggsave(pdf_path, plot, width = width, height = height, device = "pdf")
    if (file.exists(pdf_path) && file.info(pdf_path)$size > 1000) {
      cat("      [OK] PDF saved\n")
      success <- TRUE
    }
  }, error = function(e) {
    cat("      [ERROR] PDF:", conditionMessage(e), "\n")
  })
  
  # Save PNG
  png_path <- paste0(filepath_base, ".png")
  tryCatch({
    ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = 300, device = "png")
    if (file.exists(png_path) && file.info(png_path)$size > 1000) {
      cat("      [OK] PNG saved\n")
      success <- TRUE
    }
  }, error = function(e) {
    cat("      [ERROR] PNG:", conditionMessage(e), "\n")
  })
  
  return(success)
}

cat("[INFO] Helper functions defined.\n\n")

# =============================================================================
# PART 4: PROCESS ALL COHORTS
# =============================================================================

cat("=============================================================================\n")
cat("PROCESSING ALL COHORTS\n")
cat("=============================================================================\n\n")

# Storage for results
all_transitions <- list()
all_matrices <- list()
pooled_data <- list()

# Process each cohort
for (cohort in names(DATA_PATHS)) {
  cat("\n--------------------------------------------------------------------------------\n")
  cat("Processing:", cohort, "\n")
  cat("--------------------------------------------------------------------------------\n")
  
  file_path <- DATA_PATHS[[cohort]]
  
  if (is.null(file_path) || !file.exists(file_path)) {
    cat("  [SKIP] File not found:", file_path, "\n")
    next
  }
  
  # Load data
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  cat("  Loaded N =", format(nrow(df), big.mark = ","), "\n")
  
  # Prepare Sankey data
  sankey_df <- prepare_sankey_data(df, cohort)
  
  if (is.null(sankey_df) || nrow(sankey_df) == 0) {
    cat("  [SKIP] No valid transitions\n")
    next
  }
  
  cat("  Valid transitions N =", format(nrow(sankey_df), big.mark = ","), "\n")
  
  # Store for pooled analysis
  sankey_df$cohort <- cohort
  pooled_data[[cohort]] <- sankey_df
  
  # === OVERALL ===
  cat("\n  Creating Overall Sankey...\n")
  trans_overall <- create_transition_table(sankey_df, paste0(cohort, "_Overall"))
  
  if (!is.null(trans_overall)) {
    all_transitions[[paste0(cohort, "_Overall")]] <- trans_overall
    all_matrices[[paste0(cohort, "_Overall")]] <- create_transition_matrix(trans_overall)
    
    p <- create_sankey_plot(trans_overall, 
                            title = paste0(cohort, " - Health State Transitions"),
                            subtitle = "Overall")
    cat("  Saving:\n")
    save_plot(p, file.path(fig_dir, paste0("Sankey_", cohort, "_Overall")))
  }
  
  # === STRATIFIED BY LIFESTYLE ===
  if ("lifestyle_cat" %in% names(sankey_df)) {
    cat("\n  Lifestyle distribution:\n")
    print(table(sankey_df$lifestyle_cat, useNA = "always"))
    
    for (ls_cat in c("0", "1", "2", "3+")) {
      df_sub <- sankey_df[sankey_df$lifestyle_cat == ls_cat & !is.na(sankey_df$lifestyle_cat), ]
      
      if (nrow(df_sub) < 30) {
        cat("    Lifestyle", ls_cat, ": [SKIP] N =", nrow(df_sub), "\n")
        next
      }
      
      cat("\n    Lifestyle", ls_cat, ": N =", nrow(df_sub), "\n")
      
      key <- paste0(cohort, "_Cat", gsub("\\+", "plus", ls_cat))
      trans_ls <- create_transition_table(df_sub, key)
      
      if (!is.null(trans_ls)) {
        all_transitions[[key]] <- trans_ls
        all_matrices[[key]] <- create_transition_matrix(trans_ls)
        
        p <- create_sankey_plot(trans_ls,
                                title = paste0(cohort, " - Health State Transitions"),
                                subtitle = paste0("Lifestyle: ", ls_cat, " unhealthy factors"))
        cat("    Saving:\n")
        save_plot(p, file.path(fig_dir, paste0("Sankey_", key)))
      }
    }
  }
}

# =============================================================================
# PART 5: POOLED ANALYSIS
# =============================================================================

cat("\n=============================================================================\n")
cat("POOLED ANALYSIS\n")
cat("=============================================================================\n\n")

cat("DEBUG: pooled_data length =", length(pooled_data), "\n")
cat("DEBUG: pooled_data names =", paste(names(pooled_data), collapse = ", "), "\n")

if (length(pooled_data) > 0) {
  
  # Combine all cohorts - select only common columns
  cat("DEBUG: Combining pooled data...\n")
  
  # Get common columns across all data frames
  common_cols <- Reduce(intersect, lapply(pooled_data, names))
  cat("DEBUG: Common columns =", length(common_cols), "\n")
  
  # Subset each data frame to common columns before binding
  pooled_data_subset <- lapply(pooled_data, function(df) {
    df[, common_cols, drop = FALSE]
  })
  
  pooled_df <- tryCatch({
    do.call(rbind, pooled_data_subset)
  }, error = function(e) {
    cat("ERROR in rbind:", e$message, "\n")
    NULL
  })
  
  if (is.null(pooled_df) || nrow(pooled_df) == 0) {
    cat("ERROR: Failed to create pooled data!\n")
  } else {
    cat("Pooled data: N =", format(nrow(pooled_df), big.mark = ","), "\n")
    cat("Cohorts included:", paste(names(pooled_data), collapse = ", "), "\n")
    
    # Check lifestyle_cat
    if ("lifestyle_cat" %in% names(pooled_df)) {
      cat("\nPooled lifestyle distribution:\n")
      print(table(pooled_df$lifestyle_cat, useNA = "always"))
    }
    
    # === POOLED OVERALL ===
    cat("\n--- Pooled Overall ---\n")
    trans_pooled <- create_transition_table(pooled_df, "Pooled_Overall")
    
    if (!is.null(trans_pooled)) {
      all_transitions[["Pooled_Overall"]] <- trans_pooled
      all_matrices[["Pooled_Overall"]] <- create_transition_matrix(trans_pooled)
      
      p <- create_sankey_plot(trans_pooled,
                              title = "Pooled - Health State Transitions",
                              subtitle = "All Cohorts Combined")
      cat("Saving Pooled Overall:\n")
      save_plot(p, file.path(fig_dir, "Sankey_Pooled_Overall"))
    } else {
      cat("ERROR: trans_pooled is NULL!\n")
    }
    
    # === POOLED STRATIFIED BY LIFESTYLE ===
    if ("lifestyle_cat" %in% names(pooled_df)) {
      cat("\n--- Pooled by Lifestyle Category ---\n")
      
      for (ls_cat in c("0", "1", "2", "3+")) {
        df_sub <- pooled_df[pooled_df$lifestyle_cat == ls_cat & !is.na(pooled_df$lifestyle_cat), ]
        
        if (nrow(df_sub) < 50) {
          cat("  Lifestyle", ls_cat, ": [SKIP] N =", nrow(df_sub), "\n")
          next
        }
        
        cat("\n  Lifestyle", ls_cat, ": N =", format(nrow(df_sub), big.mark = ","), "\n")
        
        key <- paste0("Pooled_Cat", gsub("\\+", "plus", ls_cat))
        trans_ls <- create_transition_table(df_sub, key)
        
        if (!is.null(trans_ls)) {
          all_transitions[[key]] <- trans_ls
          all_matrices[[key]] <- create_transition_matrix(trans_ls)
          
          p <- create_sankey_plot(trans_ls,
                                  title = "Pooled - Health State Transitions",
                                  subtitle = paste0("Lifestyle: ", ls_cat, " unhealthy factors"))
          cat("  Saving:\n")
          save_plot(p, file.path(fig_dir, paste0("Sankey_", key)))
        }
      }
    } else {
      cat("WARNING: lifestyle_cat not in pooled_df!\n")
    }
  }  # End of else block for successful pooled_df
  
} else {
  cat("[WARNING] No pooled data available (pooled_data is empty)!\n")
}

# =============================================================================
# PART 6: CREATE COMPREHENSIVE EXCEL OUTPUT
# =============================================================================

cat("\n=============================================================================\n")
cat("CREATING EXCEL OUTPUT\n")
cat("=============================================================================\n\n")

# Combine all transitions
all_trans_df <- dplyr::bind_rows(all_transitions)

if (nrow(all_trans_df) > 0) {
  
  # Add metadata columns
  all_trans_df <- all_trans_df %>%
    dplyr::mutate(
      Is_Pooled = grepl("^Pooled", Group),
      Cohort = dplyr::case_when(
        grepl("^Pooled", Group) ~ "Pooled",
        grepl("^CHARLS", Group) ~ "CHARLS",
        grepl("^ELSA", Group) ~ "ELSA",
        grepl("^HRS", Group) ~ "HRS",
        grepl("^SHARE", Group) ~ "SHARE",
        grepl("^MHAS", Group) ~ "MHAS",
        TRUE ~ "Unknown"
      ),
      Lifestyle_Stratum = dplyr::case_when(
        grepl("_Overall$", Group) ~ "Overall",
        grepl("_Cat0$", Group) ~ "0",
        grepl("_Cat1$", Group) ~ "1",
        grepl("_Cat2$", Group) ~ "2",
        grepl("_Cat3plus$", Group) ~ "3+",
        TRUE ~ "Other"
      )
    )
  
  # Create summary by PPC-MM incidence
  ppcmm_states <- c("P1P2", "P1C", "P2C", "P1P2C")
  
  ppcmm_summary <- all_trans_df %>%
    dplyr::group_by(Cohort, Lifestyle_Stratum, Total_N) %>%
    dplyr::summarise(
      N_Total = sum(N),
      N_PPCMM = sum(N[followup_state %in% ppcmm_states]),
      N_P1P2 = sum(N[followup_state == "P1P2"]),
      N_P1C = sum(N[followup_state == "P1C"]),
      N_P2C = sum(N[followup_state == "P2C"]),
      N_P1P2C = sum(N[followup_state == "P1P2C"]),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Pct_PPCMM = round(N_PPCMM / N_Total * 100, 1),
      Pct_P1P2 = round(N_P1P2 / N_Total * 100, 1),
      Pct_P1C = round(N_P1C / N_Total * 100, 1),
      Pct_P2C = round(N_P2C / N_Total * 100, 1),
      Pct_P1P2C = round(N_P1P2C / N_Total * 100, 1)
    ) %>%
    dplyr::arrange(Cohort, Lifestyle_Stratum)
  
  # Create interpretation guide
  guide <- data.frame(
    Term = c("Healthy", "P1", "P2", "C", "P1P2", "P1C", "P2C", "P1P2C",
             "N", "Pct_row", "Pct_total", "N_PPCMM", "Pct_PPCMM"),
    Definition = c(
      "No physical, psychological, or cognitive condition",
      "Physical condition only (chronic disease)",
      "Psychological condition only (depression)",
      "Cognitive impairment only",
      "Physical + Psychological conditions (PPC-MM subtype)",
      "Physical + Cognitive conditions (PPC-MM subtype)",
      "Psychological + Cognitive conditions (PPC-MM subtype)",
      "All three conditions (PPC-MM subtype)",
      "Number of participants with this transition",
      "Percentage within baseline state (row percentage)",
      "Percentage of total sample",
      "Number developing any PPC-MM (2+ conditions)",
      "Percentage developing any PPC-MM"
    )
  )
  
  # Combine all transition matrices
  all_matrices_df <- dplyr::bind_rows(all_matrices)
  
  # Pooled detailed view
  pooled_detail <- all_trans_df %>%
    dplyr::filter(Cohort == "Pooled") %>%
    dplyr::select(Lifestyle_Stratum, baseline_state, followup_state, N, Pct_row, Pct_total) %>%
    dplyr::arrange(Lifestyle_Stratum, baseline_state, followup_state)
  
  # Create Excel workbook
  excel_sheets <- list(
    "Interpretation_Guide" = guide,
    "All_Transitions_Detail" = all_trans_df,
    "Transition_Matrices" = all_matrices_df,
    "PPCMM_Summary" = ppcmm_summary,
    "Pooled_Detail" = pooled_detail,
    "Pooled_Summary" = ppcmm_summary %>% dplyr::filter(Cohort == "Pooled"),
    "By_Cohort" = ppcmm_summary %>% dplyr::filter(Lifestyle_Stratum == "Overall")
  )
  
  # Save Excel
  save_to_excel(excel_sheets, "Sankey_Comprehensive_Results.xlsx")
  
  # Also save CSV files
  save_to_csv(all_trans_df, "Sankey_All_Transitions.csv")
  save_to_csv(ppcmm_summary, "Sankey_PPCMM_Summary.csv")
  
  cat("[INFO] Excel and CSV files saved.\n")
  
  # Print summary
  cat("\n--- SUMMARY ---\n")
  cat("\nPPC-MM Incidence by Lifestyle (Pooled):\n")
  pooled_sum <- ppcmm_summary %>% 
    dplyr::filter(Cohort == "Pooled") %>%
    dplyr::select(Lifestyle_Stratum, N_Total, N_PPCMM, Pct_PPCMM)
  print(as.data.frame(pooled_sum))
  
} else {
  cat("[WARNING] No transition data generated!\n")
}

# =============================================================================
# PART 7: FINAL SUMMARY
# =============================================================================

cat("\n=============================================================================\n")
cat("FINAL SUMMARY\n")
cat("=============================================================================\n\n")

# List generated files
sankey_files <- list.files(fig_dir, pattern = "Sankey_.*\\.(pdf|png)$")
cat("Generated Sankey diagrams:", length(sankey_files) / 2, "(PDF + PNG pairs)\n")

# Check for Pooled files
pooled_files <- sankey_files[grepl("Pooled", sankey_files)]
cat("Pooled Sankey diagrams:", length(pooled_files) / 2, "\n")

if (length(pooled_files) > 0) {
  cat("  Files:\n")
  for (f in unique(sub("\\.(pdf|png)$", "", pooled_files))) {
    cat("    -", f, "\n")
  }
}

cat("\nOutput directory:", fig_dir, "\n")
cat("Excel file:", file.path(OUTPUT_DIR, "Sankey_Comprehensive_Results.xlsx"), "\n")

# =============================================================================
# PART 8: GENERATE STANDALONE LEGEND FILE
# =============================================================================

cat("\n=============================================================================\n")
cat("GENERATING STANDALONE LEGEND FILE\n")
cat("=============================================================================\n\n")

# Create standalone legend
p_legend <- create_standalone_legend()

# Save legend as PDF and PNG
legend_pdf <- file.path(fig_dir, "Sankey_Legend.pdf")
legend_png <- file.path(fig_dir, "Sankey_Legend.png")

tryCatch({
  ggplot2::ggsave(legend_pdf, p_legend, width = 8, height = 3, device = "pdf")
  cat("[OK] Legend saved:", legend_pdf, "\n")
}, error = function(e) {
  cat("[ERROR] PDF:", conditionMessage(e), "\n")
})

tryCatch({
  ggplot2::ggsave(legend_png, p_legend, width = 8, height = 3, dpi = 300, device = "png")
  cat("[OK] Legend saved:", legend_png, "\n")
}, error = function(e) {
  cat("[ERROR] PNG:", conditionMessage(e), "\n")
})

cat("\n")
cat("###############################################################################\n")
cat("#                SANKEY DIAGRAM ANALYSIS COMPLETED                            #\n")
cat("###############################################################################\n")
cat("\n")
