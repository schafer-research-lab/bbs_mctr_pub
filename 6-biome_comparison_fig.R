#Compare dominant species volatility across biomes
p18<-ggplot(data=biome_metrics,aes(x=breeding_biome,
                                    y=w_mean_volatility,
                                    group=breeding_biome))+
  geom_violin(aes(fill=breeding_biome))+
  geom_boxplot(width=0.05)+
  scale_fill_manual(values=c("#4682B4", "#B44682", "#82B446"),guide="none")+
  theme_bw(base_size = 12)+
  theme(axis.title.x = element_blank())+
  labs(x="Breeding Biome",y="Dominant Species Volatility")+
  scale_y_sqrt()

#compare niche evenness across biomes
p19<-ggplot(data=biome_metrics,aes(x=breeding_biome,
                                    y=hill_simpson_niche,
                                    group=breeding_biome))+
  geom_violin(aes(fill=breeding_biome))+
  geom_boxplot(width=0.05)+
  scale_fill_manual(values=c("#4682B4", "#B44682", "#82B446"),guide = "none")+
  theme_bw(base_size = 12)+
  theme(axis.title.x = element_blank())+
  labs(x="Breeding Biome",y="Hill-Simpson Niche Diversity")

#compare geom mean abundance across biomes
p20<-ggplot(data=biome_metrics,aes(x=breeding_biome,
                                    y=geom_mean_niche,
                                    group=breeding_biome))+
  geom_violin(aes(fill=breeding_biome))+
  geom_boxplot(width=0.05)+
  scale_fill_manual(values=c("#4682B4", "#B44682", "#82B446"),guide="none")+
  theme_bw(base_size = 12)+
  theme(axis.title.x = element_blank())+
  labs(x="Breeding Biome",y="Geometric Mean Abundance")+
  scale_y_sqrt()

#compare portfolio volatility across biomes
p21<-ggplot(data=biome_metrics,aes(x=breeding_biome,
                                    y=portfolio_volatility,
                                    group=breeding_biome))+
  geom_violin(aes(fill=breeding_biome))+
  geom_boxplot(width=0.05)+
  scale_fill_manual(values=c("#4682B4", "#B44682", "#82B446"),guide="none")+
  theme_bw(base_size = 12)+
  labs(x="Breeding Biome",y="Portfolio Volatility")






#Old, subjectively ugly, scatter plots of biomes
#-----------------------------------------------------------------------------------------------
#compare portfolio volatility to dominant species volatility

ggplot(data=biome_metrics,aes(y=w_mean_volatility, 
                              x=portfolio_volatility, 
                              color=breeding_biome))+
  geom_point(size=2)+
  theme_bw(base_size = 12)+
  scale_color_manual(values=c("#4682B4", "#B44682", "#82B446"))+
  labs(x="Portfolio Volatility", y="Dominant Volatility", 
       color="Breeding Biome",title="A. PV x Dominant Volatility")

#compare portfolio volatility to niche evenness and geom abundance
ggplot(data=biome_metrics,aes(y=hill_simpson_niche, 
                              x=portfolio_volatility, 
                              color=breeding_biome))+
  geom_point(size=2)+
  theme_bw(base_size = 12)+
  scale_color_manual(values=c("#4682B4", "#B44682", "#82B446"))+
  labs(x="Portfolio Volatility", y="Hill-Simpson Niche Diversity", 
       color="Breeding Biome",title="B. PV x Niche Diversity")
ggplot(data=biome_metrics,aes(y=geom_mean_abundance, 
                                    x=portfolio_volatility, 
                                    color=breeding_biome))+
  geom_point(size=2)+
  theme_bw(base_size = 12)+
  scale_color_manual(values=c("#4682B4", "#B44682", "#82B446"))+
  labs(x="Portfolio Volatility", y="Geometric Mean Abundance", 
       color="Breeding Biome",title="C. PV x Abundance")


#compare abundance x richness by PV
ggplot(data=biome_metrics,aes(y=geom_mean_abundance, 
                              x=hill_simpson_niche, 
                              size=portfolio_volatility,
                              color=breeding_biome))+
  geom_point()+
  theme_bw(base_size = 12)+
  scale_color_manual(values=c("#4682B4", "#B44682", "#82B446"),guide="none")+
  labs(x="Hill-Simpson Niche Diversity", y="Geometric Mean Abundance", 
       color="Breeding Biome", size="Portfolio Volatility",title="D. Abundance x Niche Diversity")

#compare dominant volatility x richness by PV
ggplot(data=biome_metrics,aes(y=w_mean_volatility, 
                              x=hill_simpson_niche, 
                              size=portfolio_volatility,
                              color=breeding_biome))+
  geom_point()+
  theme_bw(base_size = 12)+
  scale_color_manual(values=c("#4682B4", "#B44682", "#82B446"),guide="none")+
  labs(x="Hill-Simpson Niche Diversity", y="Dominant Volatility", 
       color="Breeding Biome", size="Portfolio Volatility",title="E. Dominant Volatility x Niche Diversity")

