#plot the mean and distribution of pairwise species correlations each year.
df <- mctr_corr %>%
  dplyr::filter(breeding_biome == biome_val)

p17<-ggplot(df, aes(x = as.character(year), y = corr)) + 
  geom_hline(yintercept = 0, color = "black", linewidth = 1, alpha = .5) +
  geom_boxplot(fill = "blue", alpha = .75, notch = TRUE) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", panel.grid.minor = element_blank()) +
  stat_summary(fun = median, colour = "red", geom = "line", group = 1, linewidth = 2) +
  scale_x_discrete(breaks = seq(min_year, 2020, 5)) +
  labs(x = "Year", y = "Correlation in Species Returns")

#pearson's r between PV and average correlation between species
corrs<-mctr_corr%>%
  group_by(breeding_biome, year)%>%
  summarise(corr=mean(corr,na.rm = T))%>%
  ungroup()%>%
  left_join(biome_metrics%>%dplyr::select(breeding_biome,year,portfolio_volatility))

corrs2<-filter(corrs, breeding_biome==biome_val)
cor(corrs2$corr,corrs2$portfolio_volatility)
