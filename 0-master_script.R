library(AICcmodavg)
library(ggcorrplot)
library(ggplotify)
library(kableExtra)
library(patchwork)
library(psych)
library(tidyverse)
library(vcd)
library(viridis)

source("ewma_function.R")
source("ewma_weights.R")


# Trait data --------------------------------------------------------------
source("1-trait_return_data_preparation.R")

#---------------------------------------------------------------------------------------------------
# define parameters
# ewma smoothing parameter
lambda_val <- .6
# the amount of total exponential weight needed in order to keep cov matrix
weight_val <- .993
# minimum number of years needed for valid cov matrix
(num_years <- ewma_weights(lambda = lambda_val, mval = 50, weight_cutoff = weight_val))
# start year for valid varcov based on cumulative exponential weight
(min_year <- 1971 + num_years - 1)

#---------------------------------------------------------------------------------------------------
# source analysis script
source("2-community_analysis.R")

#---------------------------------------------------------------------------------------------------
# biome to plot
biome_vals <- c("Aridlands", "Eastern Forest", "Grassland")

for(biome_val in biome_vals){
  
  # a. Biome metric time series
  source("3-geom_abun_hill_simp_pv_fig.R")
  
  # b. Biome CCTR plots
  source("4-cctr_fig.R")
  
  # c. Mean and distribution of pairwise species correlations each year
  source("5-return_corr_by_MCTR_fig.R")

  #---------------------------------------------------------------------------------------------------
  # plot and save
  
  #1. Individual biome plots
  title_string <- paste0(biome_val,   
                         ',lambda = ', lambda_val 
  )
  
  #   a. Biome metric time series and CCTR plots
  g1 <- (p1 / p6) + 
    plot_annotation(title = title_string, theme = theme(plot.title = element_text(size = 14)))
  ggsave(paste0("Figures/", tolower(biome_val),  "_volatility_cctr.png"), width = 10, height = 10, dpi = "retina")
  
  #   b. Mean and distribution of pairwise species correlations each year
  g5 <- p17 +
  plot_annotation(title = title_string, theme = theme(plot.title = element_text(size = 14)))
  ggsave(paste0("Figures/", tolower(biome_val), "_correlations.png"), plot = g5, width = 16, height = 8, dpi = "retina")
  
}


# d. Boxplots of community metrics comparing biomes
source("6-biome_comparison_fig.R")

# e.-f. linear regression and correlation between species MCTR and traits & community PV and drivers
source("7-volatility_drivers.R")


#2. Plots comparing all biomes together
#   d. Boxplots of community metrics comparing biomes
g6 <- (p20 / p19 / p18 / p21) 
ggsave(paste0("Figures/biome_comparison.png"), plot = g6, width = 10, height = 10, dpi = "retina")

#   e. Plotting relationship between species MCTR and species traits
ggsave(paste0("Figures/mctr_by_reprod.png"), plot = p23, width = 7, height = 5, dpi = "retina")

#   e. Plotting relationship between community PV and community structure
g8 <- (p24 / p25 / p26)+ 
  plot_layout(guides = "collect")
ggsave(paste0("Figures/pv_by_volatility_niche_abund.png"), plot = g8, width = 7, height = 10, dpi = "retina")

