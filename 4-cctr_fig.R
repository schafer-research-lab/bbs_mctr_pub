#joined species dominance status and annual MCTR values
cctr_df <- mctr_species_years%>%
  dplyr::filter(breeding_biome == biome_val) %>%
  arrange(year, cctr) %>%
  mutate(common_name = str_remove(common_name, "'"))

# assign consistent colors across plots
myColors <- viridis(2, option="viridis")
names(myColors) <- levels(factor(cctr_df$status))

p6 <- ggplot(cctr_df, aes(y = cctr, x = year, fill = status)) +
  geom_bar(color = "white", position = "stack", stat = "identity") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  scale_fill_manual(name = "Species Dominance", values = myColors, labels = c("Dominant", "Non-Dominant")) +
  scale_x_continuous(breaks = seq(min_year, 2020, 5)) +
  geom_hline(yintercept = 0, color = "gray", linewidth = .2) +
  labs(x = "Year", y = "Total Contribution to Risk",title="B")
