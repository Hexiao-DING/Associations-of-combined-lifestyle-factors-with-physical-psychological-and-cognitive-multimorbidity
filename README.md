# Associations of Combined Lifestyle Factors with Physical-Psychological-Cognitive Multimorbidity

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Active-brightgreen.svg)]()

> **Full Title**: Associations of combined lifestyle factors with physical, psychological, and cognitive multimorbidity in older adults: a multi-cohort analysis of five prospective studies

---

## 📋 Table of Contents

- [Project Overview](#-project-overview)
- [Study Design](#-study-design)
- [Data Sources](#-data-sources)
- [Analysis Framework](#-analysis-framework)
- [Repository Structure](#-repository-structure)
- [Code Documentation](#-code-documentation)
- [How to Run](#-how-to-run)
- [Output Files](#-output-files)
- [Variable Definitions](#-variable-definitions)
- [Statistical Methods](#-statistical-methods)
- [Dependencies](#-dependencies)
- [Authors](#-authors)
- [Acknowledgments](#-acknowledgments)

---

## 🎯 Project Overview

This project investigates the **association between modifiable lifestyle factors and Physical-Psychological-Cognitive Multi-Morbidity (PPC-MM)** in older adults using harmonized data from five international prospective cohort studies.

### Key Research Questions

1. **Individual Effects**: What is the independent association of each lifestyle factor (drinking, smoking, physical inactivity, social isolation) with PPC-MM risk?
2. **Cumulative Effects**: Is there a dose-response relationship between the number of unhealthy lifestyle factors and PPC-MM incidence?
3. **Population Impact**: What proportion of PPC-MM cases could potentially be prevented by modifying lifestyle factors (Population Attributable Fraction)?

### Clinical Significance

- PPC-MM affects a substantial proportion of older adults globally
- Lifestyle factors are modifiable and represent potential intervention targets
- Understanding the cumulative impact informs public health strategies

---

## 📊 Study Design

### Design Type
- **Multi-cohort prospective study**
- **Harmonized individual participant data analysis**

### Inclusion Criteria
- Age ≥50 years at baseline
- Free of PPC-MM at baseline (i.e., <2 of the 3 health domains affected)
- Complete data on lifestyle exposures
- At least one follow-up wave available

### Exclusion Criteria
- Prevalent PPC-MM at baseline
- Missing data on key exposures or outcomes (complete case analysis)
- Severe neurological conditions at baseline

---

## 🗂 Data Sources

### Participating Cohorts

| Cohort | Full Name | Country/Region | Baseline Wave | Follow-up Waves | Sample Size |
|--------|-----------|----------------|---------------|-----------------|-------------|
| **CHARLS** | China Health and Retirement Longitudinal Study | China | Wave 1 (2011) | 2, 3, 4 | ~10,000 |
| **ELSA** | English Longitudinal Study of Ageing | England | Wave 7 (2014) | 8, 9 | ~5,000 |
| **HRS** | Health and Retirement Study | USA | Wave 10 (2010) | 11, 12, 13, 14 | ~8,000 |
| **SHARE** | Survey of Health, Ageing and Retirement in Europe | 20 European countries | Wave 4 (2011) | 5, 6, 7, 8 | ~20,000 |
| **MHAS** | Mexican Health and Aging Study | Mexico | Wave 3 (2012) | 4, 5 | ~8,000 |

### Data Harmonization

All cohorts were harmonized through the **Gateway to Global Aging Data (g2aging.org)** platform, ensuring:
- Consistent variable definitions across studies
- Comparable measurement instruments
- Standardized coding schemes

---

## 🔬 Analysis Framework

### Primary Outcomes

| Outcome | Definition | Code |
|---------|------------|------|
| **Overall PPC-MM** | Any ≥2 of 3 health domains affected | `event_ppcmm` |
| **P1P2** | Physical + Psychological | `event_mm_phys_psych` |
| **P1C** | Physical + Cognitive | `event_mm_phys_cog` |
| **P2C** | Psychological + Cognitive | `event_mm_psych_cog` |
| **P1P2C** | All three domains | `event_mm_all_three` |

### Exposure Variables

#### A. Individual Lifestyle Factors (Binary)

| Factor | Healthy (0) | Unhealthy (1) |
|--------|-------------|---------------|
| **Drinking** | Non-drinker or moderate | Any alcohol consumption |
| **Smoking** | Never smoker | Ever smoker |
| **Physical Activity** | Meets guidelines | Physically inactive |
| **Social Participation** | Socially engaged | Socially isolated |

#### B. Cumulative Lifestyle Score

| Analysis | Categories | Reference |
|----------|------------|-----------|
| **Main (4-level)** | 0, 1, 2, 3+ unhealthy factors | 0 (healthiest) |
| **Sensitivity (5-level)** | 0, 1, 2, 3, 4 unhealthy factors | 0 (healthiest) |

### Sensitivity Analyses

| ID | Analysis | Purpose |
|----|----------|---------|
| **S1** | 5-level lifestyle categories | Assess finer dose-response gradient |
| **S2** | Heavy drinking definition | Test robustness with stricter alcohol threshold |
| **S3** | MICE imputed data | Assess impact of missing data |
| **S4** | Drop first follow-up | Address reverse causality concerns |

### Covariates

- Age (continuous, at baseline)
- Sex (male/female)
- Education (primary/secondary/tertiary)
- Cohort (fixed effect in pooled analysis)

---

## 📁 Repository Structure

```
Project_Root/
│
├── 📂 Code/
│   ├── 00_Functions_and_Setup.R    # Environment setup, global functions
│   ├── 01_Pooled_Descriptive.R     # Descriptive statistics
│   ├── 02_Phi_ICC_Analysis.R       # Correlation and clustering analysis
│   ├── 03_Pooled_Cox_Analysis.R    # Cox proportional hazards models
│   ├── 04_Sankey_Diagram.R         # Health state transition diagrams
│   ├── 05_PAF_Analysis.R           # Population Attributable Fraction
│   ├── 06_Meta_Analysis.R          # Meta-analysis across cohorts
│   ├── main_analysis.R             # Master script (runs all)
│   └── README.md                   # This documentation
│
├── 📂 Data_Excel/
│   ├── 📂 CHARLS/
│   │   ├── CHARLS_main_analysis.csv
│   │   ├── CHARLS_mice_data.csv
│   │   ├── CHARLS_sensitivity_drop1st.csv
│   │   └── CHARLS_full_cleaned.csv
│   ├── 📂 ELSA/
│   │   └── ... (same structure)
│   ├── 📂 HRS/
│   │   ├── HRS_main_analysis.csv
│   │   ├── HRS_mice_data.csv
│   │   ├── HRS_sensitivity_drop1st.csv
│   │   └── HRS_full_cleaned.csv
│   ├── 📂 SHARE/
│   │   └── ... (same structure)
│   └── 📂 MHAS/
│       └── ... (same structure)
│
└── 📂 Output/
    ├── 📂 Figures/
    │   ├── 📂 Sankey/              # Sankey diagrams (PDF + PNG)
    │   ├── 📂 Forest/              # Forest plots
    │   ├── 📂 Funnel/              # Funnel plots
    │   └── 📂 DoseResponse/        # Dose-response curves
    ├── *.xlsx                      # Excel workbooks
    ├── *.csv                       # CSV data files
    └── *.rds                       # R data objects
```

---

## 📜 Code Documentation

### Script Descriptions

#### `00_Functions_and_Setup.R`
**Purpose**: Initialize the analysis environment

**Key Functions**:
- `load_cohort_data()`: Load data for specified cohort and type
- `get_lifestyle_vars()`: Get lifestyle variable names for a cohort
- `standardize_vars()`: Standardize variable names across cohorts
- `save_to_excel()`: Save results to Excel format
- `save_to_csv()`: Save results to CSV format
- `save_plot_safe()`: Robust plot saving with error handling
- `close_all_devices()`: Clean up graphics devices

**Global Settings**:
```r
BASELINE_WAVES <- list(CHARLS=1, ELSA=7, HRS=10, SHARE=4, MHAS=3)
DATA_PATHS <- list(...)      # Main analysis data paths
SENSITIVITY_PATHS <- list()  # Drop1st sensitivity data
MICE_PATHS <- list()         # MICE imputed data
```

---

#### `01_Pooled_Descriptive.R`
**Purpose**: Generate descriptive statistics and pooled dataset

**Outputs**:
- `Descriptive_Statistics.xlsx` - Multi-sheet Excel with:
  - Sample characteristics by cohort
  - Lifestyle factor distributions
  - Outcome event counts and rates
  - All 5 outcomes (Overall + 4 subtypes)
- `Pooled_main_data.rds` - Pooled analysis dataset
- `Table1_by_Cohort.csv` - Baseline characteristics table

**Key Steps**:
1. Load data from each cohort
2. Standardize variable names (wave-specific → common names)
3. Create pooled dataset
4. Calculate descriptive statistics
5. Generate TableOne summary

---

#### `02_Phi_ICC_Analysis.R`
**Purpose**: Assess correlations and clustering

**Analyses**:
1. **Phi Coefficients**: Correlation between binary lifestyle factors
2. **ICC (Intraclass Correlation)**: Clustering of outcomes within cohorts

**Outputs**:
- `Phi_ICC_Analysis_Results.xlsx`
- `Phi_Correlation_Heatmap.pdf/png`

**Interpretation**:
- Phi > 0.3: Strong correlation between lifestyle factors
- ICC > 0.05: Meaningful clustering, consider multilevel models

---

#### `03_Pooled_Cox_Analysis.R`
**Purpose**: Cox proportional hazards regression (pooled data)

**Analyses**:

| Model | Description | Exposure |
|-------|-------------|----------|
| **Primary Individual** | 4 lifestyle factors mutually adjusted | Binary (0/1) |
| **Primary Cumulative** | Lifestyle score categories | 0/1/2/3+ |
| **S1** | 5-level categories | 0/1/2/3/4 |
| **S2** | Heavy drinking definition | Uses `heavy_drink` |
| **S3** | MICE imputed data | Full imputation |
| **S4** | Drop first follow-up | Exclude wave 2 |

**Outputs**:
- `Pooled_Cox_Results_Comprehensive.xlsx` - All results
- `Pooled_Cox_All_Results.csv` - Long format results
- `DoseResponse_*.png` - Dose-response plots with P-trend

**Model Specification**:
```r
Surv(time_ppcmm_months, event_ppcmm) ~ 
  n_lifestyle_cat + age_baseline + sex + edu + strata(cohort)
```

---

#### `04_Sankey_Diagram.R`
**Purpose**: Visualize health state transitions

**State Definitions**:

| Baseline (4 states) | Follow-up (8 states) |
|---------------------|----------------------|
| Healthy (no condition) | Healthy |
| P1 (Physical only) | P1, P2, C (single) |
| P2 (Psychological only) | P1P2, P1C, P2C (dual) |
| C (Cognitive only) | P1P2C (triple) |

**Outputs**:
- Per-cohort Sankey diagrams (Overall + by lifestyle category)
- Pooled Sankey diagrams (Overall + stratified)
- `Sankey_Legend.pdf/png` - Standalone legend
- `Sankey_Comprehensive_Results.xlsx`:
  - Transition matrices with N (row %)
  - PPC-MM incidence by group
  - Detailed transition tables

**Features**:
- Labels show state name, N, and percentage
- Color-coded by health state
- Stratified by lifestyle category (0/1/2/3+)

---

#### `05_PAF_Analysis.R`
**Purpose**: Population Attributable Fraction calculation

**Formula (Miettinen)**:
```
PAF = P_case × (HR - 1) / HR
```

Where:
- `P_case` = Proportion of cases exposed
- `HR` = Hazard ratio for that exposure level

**Outputs**:
- `PAF_Analysis_Comprehensive.xlsx`
- `PAF_Analysis_All_Results.csv`

**Interpretation**: Proportion of cases theoretically preventable if exposure eliminated

---

#### `06_Meta_Analysis.R`
**Purpose**: Meta-analysis across cohorts

**Methods**:
- **Fixed-effects model**: Assumes common true effect
- **Random-effects model**: DerSimonian-Laird estimator
- **Heterogeneity**: Q statistic, I², Tau²
- **Publication bias**: Egger's test, funnel plots
- **Sensitivity**: Leave-one-out analysis

**Outputs**:
- `Meta_Analysis_Comprehensive.xlsx`:
  - Study characteristics
  - MA results (fixed + random)
  - Heterogeneity statistics
  - Leave-one-out results
  - Egger's test results
- Forest plots by outcome and exposure level
- Funnel plots for bias assessment

---

#### `main_analysis.R`
**Purpose**: Master script to run complete pipeline

**Execution Order**:
```
00_Functions_and_Setup.R
    ↓
01_Pooled_Descriptive.R
    ↓
02_Phi_ICC_Analysis.R
    ↓
03_Pooled_Cox_Analysis.R
    ↓
04_Sankey_Diagram.R
    ↓
05_PAF_Analysis.R
    ↓
06_Meta_Analysis.R
```

**Features**:
- Automatic error handling
- Progress tracking
- Analysis status summary
- Runtime logging

---

## 🚀 How to Run

### Prerequisites

1. **R version**: ≥ 4.0.0
2. **RStudio**: Recommended for interactive use
3. **Data files**: Place in `Data_Excel/[Cohort]/` folders

### Option 1: Complete Pipeline

```r
# Set working directory to Code folder
setwd("path/to/Project/Code")

# Run all analyses
source("main_analysis.R")

# Runtime: approximately 30-60 minutes depending on system
```

### Option 2: Individual Scripts

```r
setwd("path/to/Project/Code")

# Step 1: Setup (REQUIRED first)
source("00_Functions_and_Setup.R")

# Step 2: Choose analyses to run
source("01_Pooled_Descriptive.R")   # ~2 min
source("02_Phi_ICC_Analysis.R")     # ~1 min
source("03_Pooled_Cox_Analysis.R")  # ~10 min
source("04_Sankey_Diagram.R")       # ~5 min
source("05_PAF_Analysis.R")         # ~5 min
source("06_Meta_Analysis.R")        # ~10 min
```

### Option 3: Specific Cohort Analysis

```r
source("00_Functions_and_Setup.R")

# Load specific cohort
hrs_data <- load_cohort_data("HRS", type = "main")
```

---

## 📊 Output Files

### Summary Table

| File | Content | Format |
|------|---------|--------|
| `Descriptive_Statistics.xlsx` | Baseline characteristics, lifestyle distributions | Excel |
| `Phi_ICC_Analysis_Results.xlsx` | Phi coefficients, ICC values | Excel |
| `Pooled_Cox_Results_Comprehensive.xlsx` | All Cox regression results | Excel |
| `Sankey_Comprehensive_Results.xlsx` | Transition matrices, PPC-MM incidence | Excel |
| `PAF_Analysis_Comprehensive.xlsx` | Individual and cumulative PAF | Excel |
| `Meta_Analysis_Comprehensive.xlsx` | Meta-analysis results, heterogeneity | Excel |

### Figures

| Directory | Contents |
|-----------|----------|
| `Figures/Sankey/` | Sankey diagrams by cohort, pooled, and stratified |
| `Figures/Forest/` | Forest plots for each outcome × exposure level |
| `Figures/Funnel/` | Funnel plots for publication bias assessment |
| `Figures/DoseResponse/` | Dose-response plots with confidence intervals |

---

## 📖 Variable Definitions

### Cohort-Specific Naming Convention

Variables use wave prefix `wX_` where X = baseline wave:

| Cohort | Prefix | Examples |
|--------|--------|----------|
| CHARLS | `w1_` | `w1_unhealthy_drink`, `w1_unhealthy_score` |
| ELSA | `w7_` | `w7_unhealthy_drink`, `w7_unhealthy_score` |
| HRS | `w10_` | `w10_unhealthy_drink`, `w10_unhealthy_score` |
| SHARE | `w4_` | `w4_unhealthy_drink`, `w4_unhealthy_score` |
| MHAS | `w3_` | `w3_unhealthy_drink`, `w3_unhealthy_score` |

These are automatically standardized during analysis.

### Lifestyle Variables

| Variable | Type | Definition |
|----------|------|------------|
| `unhealthy_drink` | Binary | Any alcohol consumption (1) vs none (0) |
| `heavy_drink` | Binary | Heavy alcohol consumption (sensitivity) |
| `unhealthy_smoke` | Binary | Ever smoker (1) vs never (0) |
| `unhealthy_pa` | Binary | Physically inactive (1) vs active (0) |
| `unhealthy_soc` | Binary | Socially isolated (1) vs engaged (0) |
| `unhealthy_score` | Count (0-4) | Sum of 4 binary factors |
| `unhealthy_score_heavy` | Count (0-4) | Using heavy_drink instead |

### Health Domain Definitions

| Domain | Variable | Definition |
|--------|----------|------------|
| **Physical (P1)** | `physical_base` | ≥1 of 7 chronic diseases |
| **Psychological (P2)** | `psych_base` | CESD score above threshold |
| **Cognitive (C)** | `cog_base` | Age-standardized Z-score < -1.5 SD |

### Chronic Diseases (Physical Domain)
1. Hypertension
2. Diabetes
3. Cancer
4. Lung disease
5. Heart disease
6. Stroke
7. Arthritis

### PPC-MM Definition

**PPC-MM** = ≥2 of 3 health domains affected at follow-up

---

## 📈 Statistical Methods

### Cox Proportional Hazards Regression

```
h(t) = h₀(t) × exp(β₁X₁ + β₂X₂ + ... + βₖXₖ)
```

- **Time scale**: Months from baseline
- **Event**: First occurrence of PPC-MM (or subtype)
- **Censoring**: End of follow-up, death, or loss to follow-up
- **Stratification**: By cohort (to account for baseline hazard differences)

### Meta-Analysis

- **Fixed-effects**: Inverse variance weighting
- **Random-effects**: DerSimonian-Laird estimator
- **Heterogeneity**:
  - Q statistic (χ² test)
  - I² (% variation due to heterogeneity)
  - Tau² (between-study variance)

### Population Attributable Fraction

```
PAF = Σ P_case(i) × (HR(i) - 1) / HR(i)
```

Sum over all exposed categories

### P for Trend

- Linear trend test across ordered categories
- Exposure coded as continuous (0, 1, 2, 3)

---

## 📦 Dependencies

### Required R Packages

```r
# Data manipulation
tidyverse      # dplyr, tidyr, ggplot2, etc.
data.table     # Efficient large data handling

# Survival analysis
survival       # Cox models
survminer      # Survival visualization

# Mixed models
lme4           # Mixed effects models
performance    # ICC calculation

# Correlation
psych          # Phi coefficient
corrplot       # Correlation heatmaps

# Meta-analysis
meta           # Meta-analysis functions
metafor        # Advanced meta-analysis

# Tables
tableone       # Baseline characteristics
gtsummary      # Professional tables
flextable      # Export to Word/Excel

# I/O
haven          # SPSS/Stata files
writexl        # Excel output
broom          # Tidy model outputs

# Visualization
ggplot2        # Grammar of graphics
ggalluvial     # Sankey diagrams
scales         # Axis formatting
patchwork      # Combine plots
ggrepel        # Label placement
```

All packages are auto-installed if not present.

---

## 👥 Authors

**Hexiao Ding**  
The Hong Kong Polytechnic University (PolyU)

**Hongtao Cheng**  
Sun Yat-sen University (SYSU)

---

## 👨‍🏫 Supervisors

- **Prof. Jung Sun Yoo** - The Hong Kong Polytechnic University
- **Prof. Jung-E Zhang** - Sun Yat-sen University
- **Prof. Wei Xia** - Sun Yat-sen University

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

### Data Sources
- **Gateway to Global Aging Data** ([g2aging.org](https://g2aging.org))
- **CHARLS**: Peking University
- **ELSA**: University College London
- **HRS**: University of Michigan
- **SHARE**: Max Planck Institute
- **MHAS**: University of Texas Medical Branch

### Funding
[Add funding information if applicable]

### Software
- R Core Team (2024). R: A language and environment for statistical computing.
- RStudio Team (2024). RStudio: Integrated Development Environment for R.

---

*Last updated: December 2025*
