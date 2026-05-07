biome_metrics_std <- biome_metrics %>%
  group_by(breeding_biome) %>%
  mutate(gma_z = (geom_mean_abundance - mean(geom_mean_abundance)) / sd(geom_mean_abundance)) %>%
  mutate(gmn_z = (geom_mean_niche - mean(geom_mean_niche)) / sd(geom_mean_niche)) %>%
  mutate(wmv_z = (w_mean_volatility - mean(w_mean_volatility)) / sd(w_mean_volatility)) %>%
  mutate(hs_z = (hill_simpson - mean(hill_simpson)) / sd(hill_simpson)) %>%
  mutate(hsn_z = (hill_simpson_niche - mean(hill_simpson_niche)) / sd(hill_simpson_niche)) %>%
  mutate(pv_z = (portfolio_volatility - mean(portfolio_volatility, na.rm = T)) / sd(portfolio_volatility, na.rm = T)) %>%
  mutate(pr_z = (portfolio_returns - mean(portfolio_returns, na.rm = T)) / sd(portfolio_returns, na.rm = T)) %>%
  mutate(si_z = (sharpe_index - mean(sharpe_index, na.rm = T)) / sd(sharpe_index, na.rm = T)) %>%
  ungroup() %>%
  pivot_longer(cols = c(gmn_z, wmv_z, hs_z, hsn_z, pv_z, pr_z, si_z),
               names_to = "metrics", values_to = "z_score") %>%
  dplyr::select(breeding_biome, year, metrics, z_score)

#EASTERN FOREST is consistently even across species and niches
  #not impacted by evenness or dominant volatility, but correlation plays a role in PV

#GRASSLAND is consistently dominant
  #PV driven by volatility of dominant species

#ARIDLANDS has greatest variation in evenness over time
  #High vol during periods of low evenness drive PV


p1<- ggplot(
  data = biome_metrics_std %>% 
    dplyr::filter(breeding_biome == biome_val & 
                    metrics %in% c("wmv_z","hsn_z", "pv_z","gmn_z")), 
  aes(x = year, y = z_score, group = metrics, color = metrics)) +
  geom_line(linewidth = 2) +
  geom_point(size = 3) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(
    labels = c("Geometric Mean Niche Abundance", "Hill-Simpson Niche Diversity", "Portfolio Volatility","Dominant Volatility"),
    values = c("#4682B4", "#B44682", "#82B446","orange")) +
  scale_x_continuous(breaks = seq(min_year, 2020, 5)) +
  labs(x = "Year", y = "Z-Score", color = "Community Metrics:",title="A")

p2 <- ggplot(
  data = biome_metrics_std %>%
    dplyr::filter(breeding_biome == biome_val & year >= min_year & metrics %in% c("pr_z", "pv_z", "si_z")),
  aes(x = year, y = z_score, group = metrics, color = metrics)) +
  geom_line(linewidth = 2) +
  geom_point(size = 3) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_color_manual(
    labels = c("Portfolio Returns", "Portfolio Volatility", "Sharp Index"),
    values = c("#4682B4", "#82B446", "#B44682")) +
  scale_x_continuous(breaks = seq(min_year, 2020, 5)) +
  labs(x = "Year", y = "Z-Score", color = "Community Metrics:",title="B")
