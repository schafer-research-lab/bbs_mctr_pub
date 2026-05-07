# README for Code Supplement

## Quantifying Temporal Instability in Ecological Assemblages Using Modern Portfolio Theory

### Overview

This repository contains R code to reproduce all analyses and figures in the manuscript. The analysis applies modern portfolio theory metrics (portfolio volatility, marginal contribution to risk, conditional contribution to risk) to North American Breeding Bird Survey data to quantify dynamic temporal stability in ecological assemblages.

### Requirements

R version: 4.0 or higher recommended

R packages:
install.packages(c("AICcmodavg", "ggcorrplot", "ggplotify", "kableExtra", 
                   "patchwork", "psych", "tidyverse", "vcd", "viridis"))

### Repository Structure

0-master_script.R                 - Main script to reproduce all analyses
1-trait_return_data_preparation.R - Data cleaning and preparation
2-community_analysis.R            - Portfolio metrics calculation (PV, MCTR, CCTR)
3-geom_abun_hill_simp_pv_fig.R    - Time series figures for assemblage metrics
4-cctr_fig.R                      - Species CCTR contribution figures
5-return_corr_by_MCTR_fig.R       - Species correlation figures
6-biome_comparison_fig.R          - Cross-assemblage comparison figures
7-volatility_drivers.R            - Regression analyses for drivers of PV and MCTR
ewma_function.R                   - Exponentially weighted moving average functions
ewma_weights.R                    - EWMA weight calculation functions
RawData/                          - Input data files
Library/                          - Intermediate processed data objects
Figures/                          - Output figures and model summaries

### Data Sources

File                                  | Description                      | Source / DOI
--------------------------------------|----------------------------------|------------------------
Index_best_1966-2019_core_best.csv    | BBS annual abundance indices     | Pardieck et al. (2020) https://doi.org/10.5066/P9J6QUF6
BBS_1966-2019_core_best_trend.csv     | BBS trend estimates              | Pardieck et al. (2020) https://doi.org/10.5066/P9J6QUF6
rosenberg_popest.csv                  | Population size estimates        | Rosenberg et al. (2019) https://doi.org/10.1126/science.aaw1313
pigot_niche.csv                       | Foraging niche classifications   | Pigot et al. (2020) https://doi.org/10.1038/s41559-019-1070-4
wilman_bodymass.csv                   | Body mass data                   | Wilman et al. (2014) https://doi.org/10.1890/13-1917.1
Amniote_subset.csv                    | Life history traits              | Myhrvold et al. (2015) https://doi.org/10.1890/15-0846R.1
aou.csv                               | Species code                     | American Ornithological Society
SpeciesList.csv                       | Species code and common names    | https://www.sciencebase.gov/catalog/item/5ea04e9a82cefae35a129d65

### Reproducing the Analysis

1. Clone or download this repository
2. Open bbs_mctr_pub.Rproj in RStudio
3. Run 0-master_script.R

The master script will:
- Prepare trait and abundance data
- Calculate dynamic portfolio metrics using EWMA (lambda = 0.6)
- Generate all figures for the three assemblages (Aridlands, Eastern Forest, Grassland)
- Fit regression models for drivers of assemblage and species-level instability

### Key Parameters

Parameter   | Value | Description
------------|-------|-----------------------------------------------------
lambda_val  | 0.6   | EWMA smoothing parameter
weight_val  | 0.993 | Cumulative weight threshold for valid covariance matrices
min_year    | 1980  | First year with valid PV estimates after burn-in

### Outputs

Figures/ contains:
- *_volatility_cctr.png        - Time series of PV and species CCTR by assemblage
- *_correlations.png           - Distribution of pairwise species correlations by assemblage
- biome_comparison.png         - Cross-assemblage comparison of metrics
- mctr_by_reprod.png           - Species MCTR vs. reproduction regression
- pv_by_volatility_niche_abund.png - PV vs. assemblage drivers regression
- *.txt                        - Model summaries and AIC tables

Library/ contains intermediate .rda files for reproducibility without rerunning data preparation.

