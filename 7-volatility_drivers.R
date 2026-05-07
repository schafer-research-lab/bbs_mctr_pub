#Species level volatility
#---------------------------------------------------------------------------------------------
#merge species traits and MCTR values
mctr_trait_years<-left_join(select(mctr_species_years,scientific_name,year,var_r,mctr,w),
                      select(traits, scientific_name, year, breeding_biome, 
                             reprod_rate, BodyMass_Value,longevity_y))

#create a stationary dataset (average species' volatility and weight over years)
mctr_trait<-mctr_trait_years%>%
  group_by(scientific_name,breeding_biome)%>%
  select(-year)%>%
  summarise_all(mean,na.rm=T)%>%
  ungroup()

#remove species with missing trait data
mctr_trait_sub<-mctr_trait%>%
  filter(!is.na(mctr_trait$reprod_rate) &
        !is.na(mctr_trait$BodyMass_Value) &
        !is.na(mctr_trait$longevity_y))

summarise(group_by(mctr_trait_sub,breeding_biome), count=n())
#reduces n species from 30 to 29 for grassland, 56 to 55 for eastern forest, 58 to 53 for aridlands


#Relationship between species mean MCTR and species traits:
#   average annual reproduction, average longevity, average body mass,
#   mean return volatility, mean weight (dominance)

#scale all predictors
scaled_dat<-mctr_trait_sub%>%
  #log transform skewed distributions
  mutate(BodyMass_Value=log(BodyMass_Value),
         w=log(w),
         var_r=log(var_r))%>%
  mutate(across(c("var_r","w","BodyMass_Value","reprod_rate","longevity_y"),scale))

#check correlation among predictors
pairs.panels(scaled_dat%>%select(-c("scientific_name","breeding_biome")),
             scale=F, smooth=FALSE, ellipses=FALSE,pch=21, hist.col="darkgray",
             method="spearman", cex=1.2, stars=TRUE)

#linear regression
m1<-lm(mctr ~ w + var_r + BodyMass_Value + reprod_rate + longevity_y+breeding_biome, 
           data = scaled_dat)

m2<-lm(mctr ~ w + var_r + reprod_rate + longevity_y + breeding_biome, 
       data = scaled_dat)

m3<-lm(mctr ~ w + var_r + reprod_rate + breeding_biome, 
       data = scaled_dat)

m4<-lm(mctr ~ w + var_r + reprod_rate , 
       data = scaled_dat)

models <- list(m1,m2,m3,m4)

mod.names <- c('w + var_r + BodyMass_Value + reprod_rate + longevity_y + breeding_biome',
               'w + var_r + reprod_rate + longevity_y + breeding_biome', 
               'w + var_r + reprod_rate + breeding_biome', 
               'w + var_r + reprod_rate ')

#AIC model comparison. Annual reproduction is most important for predicting mean MCTR
sink(file="Figures/AIC_traits_mctr_models.txt")
print(aictab(cand.set = models, modnames = mod.names))
sink()

sink(file="Figures/lm_summary_reprod_mctr.txt")
print(summary(m1))
sink()

ci<-confint(m1)

#models using Reproductive rate and longevity, dominance and volatility have good performance. 
#More volatile spp have higher reprod rates 

#Plot relationship between reproduction and longevity with MCTR
#using regression line accounting for all predictors
p22<-ggplot(data = mctr_trait_sub, 
      mapping =aes(x=longevity_y, y = mctr)) +
  geom_point(aes(color=breeding_biome),size=2)+
  stat_smooth(method = "lm", mapping=aes(y=predict(m3,scaled_dat)), linewidth = 1, col="darkgray") + 
  scale_color_manual(values=c("#4682B4", "#B44682", "#82B446"))+
  labs(x="Longevity (Years)", y="Species MCTR",color="Breeding Biome", title="A. Species MCTR and Longevity") +
  theme_bw(base_size = 12)

p23<-ggplot(data = mctr_trait_sub, 
           mapping =aes(x=reprod_rate, y = mctr)) +
  geom_point(aes(color=breeding_biome),size=2)+
  stat_smooth(method = "lm", mapping=aes(y=predict(m3,scaled_dat)), linewidth = 1, col="darkgray") +
  scale_color_manual(values=c("#4682B4", "#B44682", "#82B446"))+
  labs(x="Annual Reproduction (clutch size * clutches per year)", y="Species MCTR",color="Breeding Biome")+#, title="B. Species MCTR and Annual Reproduction") +
  theme_bw(base_size = 12)




#Portfolio Level Volatility
#---------------------------------------------------------------------------------------------
#use z-scores to scale predictors

biome_metrics_long <- biome_metrics %>%
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
  select(c("breeding_biome","year","gmn_z","hsn_z","wmv_z","pv_z"))


#Linear regression of predictors
m1.2<-lm(pv_z ~ gmn_z+hsn_z+wmv_z+breeding_biome, 
       data = biome_metrics_long)

#try adding abundance as an interaction? 
#Higher abundances result in higher volatility and 
#changes in abundance could indicate a shift in dominance.
m2.2<-lm(pv_z ~ (hsn_z*gmn_z)+(wmv_z*gmn_z)+breeding_biome, 
       data = biome_metrics_long)

m3.2<-lm(pv_z ~ gmn_z+hsn_z+wmv_z, 
       data = biome_metrics_long)

m4.2<-lm(pv_z ~ (hsn_z*gmn_z)+(wmv_z*gmn_z), 
       data = biome_metrics_long)

models <- list(m1.2,m2.2,m3.2,m4.2)

mod.names <- c('Niche Abundance + Niche Diversity + Dominant Volatility + Biome',
               'Niche Diversity*Niche Abundance + Dominant Volatility*Niche Abundance + Biome', 
               'Niche Abundance + Niche Diversity + Dominant Volatility',
               'Niche Diversity*Niche Abundance + Dominant Volatility*Niche Abundance')

sink(file="Figures/AIC_PV_models.txt")
print(t<-aictab(cand.set = models, modnames = mod.names))
sink()
#accounting for biome decreases model performance. 
#Independent and interaction models have similar performance, but independent model is slightly better

sink(file="Figures/lm_summary_PV.txt")
print(summary(m3.2))
sink()

ci<-confint(m3.2)


p24<-ggplot(data = biome_metrics_long, 
       mapping =aes(x=wmv_z, y = pv_z)) +
  geom_point(aes(color=breeding_biome),size=2)+
  stat_smooth(method = "lm", mapping=aes(y=predict(m3.2,biome_metrics_long)), linewidth = 1, col="darkgray") +
  scale_color_manual(values=c("#4682B4", "#B44682", "#82B446"))+
  labs(x="Volatility of Dominant Species", y="Portfolio Volatility",color="Breeding Biome",title="A. Dominant Volatility and PV") +
  theme_bw(base_size = 12)+
  theme(legend.position = "none")

p25<-ggplot(data = biome_metrics_long, 
           mapping =aes(x=hsn_z, y = pv_z)) +
  geom_point(aes(color=breeding_biome),size=2)+
  stat_smooth(method = "lm",  mapping=aes(y=predict(m3.2,biome_metrics_long)), linewidth = 1, col="darkgray") +
  scale_color_manual(values=c("#4682B4", "#B44682", "#82B446"))+
  labs(x="Hill-Simpson Niche Diversity", y="Portfolio Volatility",color="Breeding Biome",title="B. Hill-Simpson Niche Diversity and PV") +
  theme_bw(base_size = 12)

p26<-ggplot(data = biome_metrics_long, 
            mapping =aes(x=gmn_z, y = pv_z)) +
  geom_point(aes(color=breeding_biome),size=2)+
  stat_smooth(method = "lm",  mapping=aes(y=predict(m3.2,biome_metrics_long)), linewidth = 1, col="darkgray") +
  scale_color_manual(values=c("#4682B4", "#B44682", "#82B446"))+
  labs(x="Geometric Mean Species Abundance", y="Portfolio Volatility",color="Breeding Biome",title="C. Niche Abundance and PV") +
  theme_bw(base_size = 12)

