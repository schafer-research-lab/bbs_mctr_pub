# make needed empty objects
biomes <- sort(unique(bbs_abundance$breeding_biome))
define_it <- function() {
  x <- vector(mode = "list", length = length(biomes))
  names(x) <- biomes
  return(x)
}
r_mat <- define_it()
a_mat <- define_it()
w_mat <- define_it()
cov_mat <- define_it()
corr_mat <- define_it()
var_mat <- define_it()
mctr_mat <- define_it()
cctr_mat <- define_it()
temp <- temp2 <- temp3 <- temp4 <- define_it()

#---------------------------------------------------------------------------------------------------
# compute volatility related metrics
for (i in 1:length(biomes)) {
  # matrix of returns (intrinsic growth rates), rows are years and columns are species
  r_m <- bbs_abundance %>%
    dplyr::filter(year > 1970 & breeding_biome == biomes[i]) %>%  
    dplyr:: select(scientific_name, year, r) %>% 
    pivot_wider(id_cols = year, names_from = scientific_name, values_from = r) %>% 
    tibble::column_to_rownames('year') %>%
    data.matrix(.)
  # matrix of abundances, rows are years and columns are abundances
  a_m <-  bbs_abundance %>%
    dplyr::filter(year > 1970 & breeding_biome == biomes[i]) %>%  
    dplyr:: select(scientific_name, year, abundance) %>% 
    pivot_wider(id_cols = year, names_from = scientific_name, values_from = abundance) %>% 
    tibble::column_to_rownames('year') %>%
    data.matrix(.)
  # matrix of species weights using abundance, rows are years and columns are species
  w_m <- a_m / rowSums(a_m)
  # vector of yearly returns r weighted by species weights
  p_ret <- as.vector(rowSums(r_m * w_m))
  # list where each element t has ewma cov matrix based on species returns from 1970 to year t 
  cov_m <- ewma(rtn = r_m, lambda = lambda_val, increment = 1E-9, list_out = TRUE)
  p_vol <- NA
  mctr_m <- cctr_m <- var_m <- matrix(NA, nrow = dim(r_m)[1], ncol = dim(r_m)[2])
  rownames(mctr_m) <- rownames(var_m) <- rownames(r_m)
  colnames(mctr_m) <- colnames(var_m) <- colnames(r_m)
  for (t in 1:dim(r_m)[1]) {
  # vector of ewma portfolio volatility, where t is year
    p_vol[t] <- sqrt(t(w_m[t, ]) %*% cov_m[[t]] %*% w_m[t, ])
    # matrix of marginal contribution to risk MCTR which is sum of species weight times covariance 
    # between focal species and all other species including itself dividing by portfolio volatility
    # where rows are years columns are species
    mctr_m[t, ] <- (cov_m[[t]] %*% w_m[t, ]) / p_vol[t]
    # matrix of species variances where rows are years columns are species
    var_m[t, ] <- diag(cov_m[[t]])
  }
  # matrix of total contribution to risk CCTR by species which is weight * MCTR
  # where rows are years and columns are species
  cctr_m <- w_m * mctr_m
  # make data frame of portfolio volatility and returns
  temp[[i]] <- data.frame(
    breeding_biome = biomes[i],
    year = as.integer(rownames(r_m)), 
    portfolio_volatility = p_vol, 
    portfolio_returns = p_ret)
  # compute # of mctr values in stabilizer, low, or high risk groups for each species 
  # across all years which have valid pv and msctr values using hr_val cutoff
  # create species x year frame with w, mctr, cctr, and mctr category
  var_df <- as.data.frame(var_m) %>%
    mutate(year = as.numeric(row.names(var_m))) %>%
    dplyr::filter(year >= min_year) %>%
    pivot_longer(!year, values_to = "var_r", names_to = "scientific_name")
  w_df <- as.data.frame(w_m) %>%
    mutate(year = as.numeric(row.names(w_m))) %>%
    dplyr::filter(year >= min_year) %>%
    pivot_longer(!year, values_to = "w", names_to = "scientific_name")
  mctr_df <- as.data.frame(mctr_m) %>%
    mutate(year = as.numeric(row.names(mctr_m))) %>%
    dplyr::filter(year >= min_year) %>%
    pivot_longer(!year, values_to = "mctr", names_to = "scientific_name") %>%
    dplyr::select(scientific_name, year, mctr)  
  temp2[[i]] <- as.data.frame(cctr_m) %>%
    mutate(year = as.numeric(row.names(cctr_m))) %>%
    dplyr::filter(year >= min_year) %>%
    pivot_longer(!year, values_to = "cctr", names_to = "scientific_name") %>%
    full_join(mctr_df, by = c("scientific_name", "year")) %>%
    full_join(w_df, by = c("scientific_name", "year")) %>%
    full_join(var_df, by = c("scientific_name", "year")) %>%
    mutate(breeding_biome = biomes[i]) %>%
    dplyr::select(breeding_biome, scientific_name, year, mctr, w, var_r, cctr) 
  corr_list <- vector(mode = "list", length = (dim(r_m)[1]))
  for (t in 1:dim(r_m)[1]) {
    corr_m <- cov2cor(cov_m[[t]])
    combos <- t(combn(colnames(corr_m), 2))
    corr_list[[t]] <- data.frame(row = combos[, 1], col = combos[, 2], corr = corr_m[combos])
    corr_list[[t]]$year <- as.numeric(row.names(r_m)[t])
  }  
  temp3[[i]] <- do.call(rbind, corr_list) %>%
    mutate(breeding_biome = biomes[i]) %>%
    dplyr::select(breeding_biome, year, scientific_name1 = row, scientific_name2 = col, corr) %>%
    dplyr::filter(year >= min_year) %>%
    arrange(year, scientific_name1, scientific_name2)
  # store output in lists
  r_mat[[i]] <- r_m
  a_mat[[i]] <- a_m
  w_mat[[i]] <- w_m
  mctr_mat[[i]] <- mctr_m
  cctr_mat[[i]] <- cctr_m
  cov_mat[[i]] <- cov_m
}
#---------------------------------------------------------------------------------------------------
# create mctr species x years df
mctr_species_years <- as.data.frame(do.call(rbind, temp2)) %>%
  left_join(bbs_abundance %>% 
    dplyr::select(breeding_biome, scientific_name, common_name, AOU, species_code, year, r, abundance), 
    by = c("breeding_biome", "scientific_name", "year")) %>%
  dplyr::select(breeding_biome, scientific_name, common_name, AOU, species_code, year, r, var_r, 
    abundance, mctr, w, cctr) %>%
  arrange(breeding_biome, scientific_name, year)

#---------------------------------------------------------------------------------------------------
# create mctr corr df
mctr_corr <- as.data.frame(do.call(rbind, temp3)) %>%
  left_join(mctr_species_years %>% dplyr::select(scientific_name1 = scientific_name, year, w1 = w, 
    common_name1 = common_name, AOU1 = AOU, species_code1 = species_code), 
    by = c("scientific_name1", "year")) %>%
  left_join(mctr_species_years %>% dplyr::select(scientific_name2 = scientific_name, year, w2 = w, 
    common_name2 = common_name, AOU2 = AOU, species_code2= species_code), 
    by = c("scientific_name2", "year")) %>%
  mutate(w = w1 * w2) %>%
  dplyr::select(breeding_biome, year, scientific_name1, scientific_name2, common_name1, common_name2, 
    AOU1, AOU2, species_code1, species_code2, w1, w2, w, corr)

#---------------------------------------------------------------------------------------------------
# make with biome metrics by computing traditional community metrics
# and joining portfolio volatility and returns 
temp <- as.data.frame(do.call(rbind, temp)) 
biome_metrics <- bbs_abundance %>%
  arrange(breeding_biome, year) %>%
  group_by(breeding_biome, year) %>%
  summarise(
    geom_mean_abundance = exp(mean(log(abundance))),
    hill_simpson = 1 / sum((abundance / (sum(abundance)))^2)
  ) %>%
  ungroup() %>%
  full_join(temp, by = c("breeding_biome", "year")) %>%
  dplyr::filter(year >= min_year) %>%
  mutate(sharpe_index = portfolio_returns / portfolio_volatility) 

#---------------------------------------------------------------------------------------------------
# add additional biome metrics using community stability hypotheses from Loreau et al 2021, Gao et al 2021

# Hypoth 1. Mean volatility of dominant species drives community volatility
# Calculate weighted mean volatility of dominant species

#   a. define and select dominant species
dominant<-mctr_species_years%>%
  group_by(year,breeding_biome)%>%
  arrange(desc(w))%>%
  mutate(dominant_quant=case_when(
    #1. select species with weights in the 90th percentile in each year/biome
    w>quantile(w, probs = 0.9)~1,
    w<quantile(w, probs = 0.9)~0),
    #2. select species with top 3 weights in each year/biome
    dominant=case_when(
      w>=w[3]~1,
      w<w[3]~0))

# Save table of dominant species summaries: number of years dominant

dom_tab <- dominant |> 
  group_by(breeding_biome, species_code) |> 
  summarize(n = sum(dominant),
            avg_weight = formatC(mean(w), digits = 2, format = "e"),
            avg_cctr = formatC(mean(cctr), digits = 2, format = "e")) |> 
  arrange(breeding_biome, desc(n), species_code) 

# uncomment to produce latex table code
# dom_tab |>  
#   knitr::kable(
#   "latex",
#   booktabs = TRUE,
#   longtable = TRUE,
#   digits = 2,
#   col.names = c("Breeding Biome", "Species Code", "n", "Avg. w", "Avg. CCTR"),
#   escape = FALSE
#   ) |>
#   collapse_rows(
#     columns = 1,
#     latex_hline = "custom",
#     valign = "top",
#     custom_latex_hline = 2,
#     longtable_clean_cut = FALSE
#   ) |>
#   kable_styling(latex_options = c("hold_position", "repeat_header"))

dominant<-dominant%>%
  #use the "top 3" selection method to select dominant species
  filter(dominant==1)

wmv<-dominant%>%
  arrange(breeding_biome, year) %>%
  group_by(breeding_biome, year) %>%
  summarise(w_mean_volatility = sum(var_r*w)) %>%
  ungroup()


# Hypoth 2. Niche differentiation drives community stability (species asynchrony)
# Calculate hill-simpson evenness of foraging niches
  
niche <- group_by(traits,ForagingNiche,year,breeding_biome)%>%
#   a. sum annual abundance within each niche (use abundance, not species since species don't change)
  summarize(abundance=sum(abundance,na.rm=T))%>%
  ungroup()%>%
  arrange(breeding_biome, year) %>%
  group_by(breeding_biome, year) %>%
#   b. calculate evenness and geometric mean abundance across niches
  summarise(hill_simpson_niche = 1 / sum((abundance / (sum(abundance)))^2),
            geom_mean_niche = exp(mean(log(abundance)))) %>%
  ungroup() 

# join with biome metrics
biome_metrics <- biome_metrics %>%
  left_join(niche, by = c("breeding_biome", "year"))%>%
  left_join(wmv, by = c("breeding_biome", "year"))

#---------------------------------------------------------------------------------------------------
#add annual species dominance to data frame of annual mctr and variance

mctr_species_years <- mctr_species_years%>%
  left_join(dominant[,c("scientific_name","year","dominant")],by=c("scientific_name","year"))%>%
  mutate(status=case_when(
    dominant==1~"Dominant Species",
    is.na(dominant)~"Non-Dominant Species"
  ))%>%
  dplyr::select(-dominant)

#---------------------------------------------------------------------------------------------------
# save computed objects to library
save(biome_metrics, file = "Library/biome_metrics.rda")
save(mctr_species_years, file = "Library/mctr_species_years.rda")
save(mctr_corr, file = "Library/mctr_corr.rda")
save(r_mat, file = "Library/r_mat.rda")
save(a_mat, file = "Library/a_mat.rda")
save(w_mat, file = "Library/w_mat.rda")
save(mctr_mat, file = "Library/mctr_mat.rda")
save(cctr_mat, file = "Library/cctr_mat.rda")
save(cov_mat, file = "Library/cov_mat.rda")
