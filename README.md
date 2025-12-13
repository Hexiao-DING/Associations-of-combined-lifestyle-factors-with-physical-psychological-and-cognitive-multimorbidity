# Associations-of-combined-lifestyle-factors-with-physical-psychological-and-cognitive-multimorbidity
Associations of combined lifestyle factors with physical, psychological, and cognitive multimorbidity in older adults: a multi-cohort analysis of five prospective studies
# Lifestyle Factors and PPC-MM Association Study

## Project Overview

This project investigates the association between lifestyle factors (drinking, smoking, physical activity, and social participation) and Physical-Psychological-Cognitive Multi-Morbidity (PPC-MM) using data from five harmonized cohort studies:

| Cohort | Full Name | Country | Baseline Wave |
|--------|-----------|---------|---------------|
| CHARLS | China Health and Retirement Longitudinal Study | China | Wave 1 |
| ELSA | English Longitudinal Study of Ageing | England | Wave 7 |
| HRS | Health and Retirement Study | USA | Wave 10 |
| SHARE | Survey of Health, Ageing and Retirement in Europe | Europe | Wave 4 |
| MHAS | Mexican Health and Aging Study | Mexico | Wave 3 |

## Analysis Framework (Lancet Standard)

### Outcomes
- **Primary**: Overall PPC-MM (any ≥2 conditions)
- **Secondary**: 4 subtypes (P1P2, P1C, P2C, P1P2C)

### Exposures
**A. Individual Lifestyle Factors (mutually adjusted)**
- Drinking, Smoking, Physical Inactivity, Social Isolation

**B. Cumulative Effect**
- Main: 4-level (0/1/2/3+)
- Sensitivity: 5-level (0/1/2/3/4)

### Sensitivity Analyses
| Analysis | Description |
|----------|-------------|
| S1 | 5-level categories (0/1/2/3/4) |
| S2 | Heavy drinking definition |
| S3 | MICE imputed data |
| S4 | Drop first follow-up wave |

## Repository Structure

```
Code/
├── 00_Functions_and_Setup.R   # Common functions and environment setup
├── 01_Pooled_Descriptive.R    # Pooled data and descriptive statistics
├── 02_Phi_ICC_Analysis.R      # Phi coefficient and ICC analysis
├── 03_Pooled_Cox_Analysis.R   # Pooled Cox regression (comprehensive)
├── 04_Sankey_Diagram.R        # Health state transition Sankey diagrams
├── 05_PAF_Analysis.R          # Population Attributable Fraction analysis
├── 06_Meta_Analysis.R         # Meta-analysis across cohorts
├── main_analysis.R            # Master script to run all analyses
└── README.md                  # This file

Data_Excel/
├── CHARLS/
│   ├── CHARLS_main_analysis.csv
│   ├── CHARLS_mice_data.csv
│   ├── CHARLS_sensitivity_drop1st.csv
│   └── CHARLS_full_cleaned.csv
├── ELSA/
│   └── ... (same structure)
├── HRS/
│   ├── HRS_main_analysis.csv      # Main analysis (baseline=w10)
│   ├── HRS_mice_data.csv          # MICE imputed data
│   ├── HRS_sensitivity_drop1st.csv # Sensitivity (drop first follow-up)
│   └── HRS_full_cleaned.csv       # Full cleaned data
├── SHARE/
│   └── ... (same structure)
└── MHAS/
    └── ... (same structure)

Output/
├── Figures/
│   ├── Sankey/                  # Sankey diagrams
│   ├── Forest/                  # Forest plots
│   └── Funnel/                  # Funnel plots
├── *.xlsx                       # Excel results
├── *.csv                        # CSV results
└── *.rds                        # R data objects
```

## How to Run

### Option 1: Run Complete Pipeline
```r
setwd("path/to/Code")
source("main_analysis.R")
```

### Option 2: Run Individual Scripts
```r
setwd("path/to/Code")
source("00_Functions_and_Setup.R")  # Always run first
source("01_Pooled_Descriptive.R")
source("02_Phi_ICC_Analysis.R")
source("03_Pooled_Cox_Analysis.R")
source("04_Sankey_Diagram.R")
source("05_PAF_Analysis.R")
source("06_Meta_Analysis.R")
```

## Output Files

### 1. Descriptive Statistics
- `Descriptive_Statistics.xlsx`
  - Lifestyle and outcome distributions by cohort
  - All_Outcomes_Summary: Overall + 4 subtypes event counts
- `Table1_by_Cohort.csv`

### 2. Phi and ICC Analysis
- `Phi_ICC_Analysis_Results.xlsx`
  - Phi coefficients for lifestyle variable correlations
  - ICC for clustering by cohort (all 5 outcomes)
- `Phi_Correlation_Heatmap.pdf/png`

### 3. Cox Regression (Pooled)
- `Pooled_Cox_Results_Comprehensive.xlsx`
  - Primary_Individual: 4 factors mutually adjusted (5 outcomes)
  - Primary_Cumulative: 4-level (0/1/2/3+) with P-trend
  - S1_5Level, S2_HeavyDrink, S3_MICE, S4_Drop1st sheets
- `Pooled_Cox_All_Results.csv`
- `DoseResponse_*.png`: Dose-response plots

### 4. Sankey Diagrams
- `Sankey_Comprehensive_Results.xlsx`
  - Interpretation_Guide: Term definitions
  - All_Transitions_Detail: N, row %, total % for all transitions
  - Transition_Matrices: Wide format with N (%)
  - PPCMM_Summary: Incidence by cohort and lifestyle
  - Pooled_Detail/Summary, By_Cohort sheets
- `Sankey_All_Transitions.csv`, `Sankey_PPCMM_Summary.csv`
- `Figures/Sankey/`
  - `Sankey_[Cohort]_Overall.pdf/png`
  - `Sankey_[Cohort]_Cat[0/1/2/3plus].pdf/png`
  - `Sankey_Pooled_Overall.pdf/png`
  - `Sankey_Pooled_Cat[0/1/2/3plus].pdf/png`
  - `Sankey_Legend.pdf/png` (standalone legend file)

### 5. PAF Analysis
- `PAF_Analysis_Comprehensive.xlsx`
  - Primary_Individual: Each factor's PAF
  - Primary_Cumulative: Combined PAF by category
  - Sensitivity analysis sheets (S2/S3/S4)
- `PAF_Analysis_All_Results.csv`

### 6. Meta-Analysis
- `Meta_Analysis_Comprehensive.xlsx`
  - Study_Characteristics: Sample sizes, demographics
  - MA_Summary_All: Fixed + Random HR, I², Q, Tau²
  - Primary, S1-S4 sheets
  - Leave_One_Out, Eggers_Test sheets
- `Study_Specific_HR_Results.csv`
- `Meta_Analysis_Summary.csv`
- `Figures/Forest/`: Forest plots by outcome × level
- `Figures/Funnel/`: Funnel plots for k≥3 analyses

## Sankey Diagram States

### Baseline States (4 categories)
| State | Definition |
|-------|------------|
| Healthy | No physical, psychological, or cognitive condition |
| P1 | Physical condition only |
| P2 | Psychological condition only |
| C | Cognitive impairment only |

### Follow-up States (8 categories)
| State | Definition |
|-------|------------|
| Healthy | No conditions |
| P1, P2, C | Single conditions |
| P1P2 | Physical + Psychological (PPC-MM) |
| P1C | Physical + Cognitive (PPC-MM) |
| P2C | Psychological + Cognitive (PPC-MM) |
| P1P2C | All three conditions (PPC-MM) |

## Variable Definitions

### Cohort-Specific Variable Naming Convention

Variables are named with wave prefix `wX_` where X is the baseline wave number:

| Cohort | Wave Prefix | Example Variables |
|--------|-------------|-------------------|
| CHARLS | w1_ | w1_unhealthy_drink, w1_unhealthy_score |
| ELSA | w7_ | w7_unhealthy_drink, w7_unhealthy_score |
| HRS | w10_ | w10_unhealthy_drink, w10_unhealthy_score |
| SHARE | w4_ | w4_unhealthy_drink, w4_unhealthy_score |
| MHAS | w3_ | w3_unhealthy_drink, w3_unhealthy_score |

These are automatically standardized to common names (e.g., `unhealthy_drink`) during analysis.

### Lifestyle Variables (Binary: 0 = Healthy, 1 = Unhealthy)
| Variable | Definition |
|----------|------------|
| unhealthy_drink | Any alcohol consumption |
| heavy_drink | Heavy alcohol consumption (sensitivity) |
| unhealthy_smoke | Ever smoker |
| unhealthy_pa | Physical inactivity |
| unhealthy_soc | Social isolation |
| unhealthy_score | Sum of 4 factors (0-4) |
| unhealthy_score_heavy | Sum using heavy_drink instead of unhealthy_drink |

### PPC-MM Definition
≥2 of the following 3 domains:
1. **Physical disease**: Any of 7 chronic diseases
2. **Depression**: CESD score above threshold
3. **Cognitive impairment**: Age-standardized Z-score < -1.5

## Statistical Methods

- **Cox Regression**: Proportional hazards, cohort-stratified or random effect
- **Meta-Analysis**: Fixed + Random effects (DerSimonian-Laird)
- **Heterogeneity**: Q statistic, I², Tau²
- **Sensitivity**: Leave-one-out, Egger's test
- **PAF**: Miettinen formula: P_case × (HR - 1) / HR

## Dependencies

Required R packages (auto-installed):
```
tidyverse, data.table, survival, survminer, lme4, performance,
psych, corrplot, meta, metafor, tableone, gtsummary, flextable,
haven, writexl, broom, ggplot2, ggalluvial, scales, patchwork
```

## Authors
Hexiao Ding (PolyU), Hongtao Cheng (SYSU)

## Supervisors
Prof. Jung Sun Yoo (PolyU), Prof. Jung-E Zhang (SYSU), Prof. Wei Xia (SYSU)

## License

[Specify license]

## Acknowledgments

- Gateway to Global Aging Data (g2aging.org)
- Individual cohort studies (CHARLS, ELSA, HRS, SHARE, MHAS)
