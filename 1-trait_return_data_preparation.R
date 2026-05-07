library(tidyverse)

# get species names and AOU from BBS
# https://www.sciencebase.gov/catalog/item/5ea04e9a82cefae35a129d65
species_list <- read.csv("RawData/SpeciesList.csv", sep = ",", header = T, stringsAsFactors = F) %>%
  mutate(scientific_name = paste(Genus, Species, sep = " ")) %>%
  dplyr::select(AOU, scientific_name, common_name = English_Common_Name)

# load BBS trend to identify trends with low credibility
bbs_bad_trends <- read.csv("RawData/BBS_1966-2019_core_best_trend.csv", sep = ",", header = T, stringsAsFactors = F) %>%
  dplyr::filter(str_trim(Credibility.Code) == "R" & Region == "SU1") %>%
  dplyr::select(AOU) %>%
  distinct() %>%
  mutate(flag = 1)

# join with BBS index data
bbs_index <- read.csv("RawData/Index_best_1966-2019_core_best.csv", sep = ",", header = T, stringsAsFactors = F) %>%
  # keep only full US non-aggregated species after 1970
  dplyr::filter(Region == "SU1" & Year >= 1970 & AOU < 30000) %>%
  left_join(species_list, by = "AOU") %>%
  dplyr::select(AOU, scientific_name, common_name, region = Region, year = Year, index = Index) %>%
  # remove species flagged as red (unreliable data) in BBS trends data
  left_join(bbs_bad_trends, by = "AOU") %>%
  dplyr::filter(is.na(flag))

# get breeding biome and population abundances from rosenberg et al 2019
popest <- read.csv("RawData/rosenberg_popest.csv", sep = ",", header = T, stringsAsFactors = F) %>%
  rename(common_name = species) %>%
  dplyr::select(common_name, breeding_biome = Breeding.Biome, first_year_popest, last_year_popest, popest)

# join with BBS index data
bbs_index <- left_join(bbs_index, popest, by = "common_name")

# Dendrocygna autumnalis Black-bellied Whistling-Duck is not in popest and is therefore dropped
# this species is a wetland species anyway which is not a breeding biome we are using
bbs_index %>%
  dplyr::filter(is.na(breeding_biome)) %>%
  dplyr::select(AOU, scientific_name, common_name, breeding_biome) %>%
  distinct()

# by keeping only species that are in popest and the eastern forest, grassland, and aridland biomes
bbs_index <- bbs_index %>%
  dplyr::filter(breeding_biome %in% c("Eastern Forest", "Aridlands", "Grassland"))

# get species codes from AOU csv file
# https://www.pwrc.usgs.gov/bbl/manual/speclist.cfm
species_codes <- read.csv("RawData/aou.csv", sep = ",", header = T, stringsAsFactors = F) %>%
  dplyr::select(AOU = Species.Number, species_code = Alpha.Code)

bbs_index <- left_join(bbs_index, species_codes, by = "AOU") %>%
  # manually add missing species codes from https://www.birdpop.org/pages/birdSpeciesCodes.php 
  mutate(species_code = replace(species_code, AOU == 2890, "NOBO")) %>%
  mutate(species_code = replace(species_code, AOU == 2930, "SCQU")) %>%
  mutate(species_code = replace(species_code, AOU == 2940, "CAQU")) %>%
  mutate(species_code = replace(species_code, AOU == 2950, "GAQU")) %>%
  mutate(species_code = replace(species_code, AOU == 3080, "STGR")) %>%
  mutate(species_code = replace(species_code, AOU == 3090, "GRSG")) %>%
  mutate(species_code = replace(species_code, AOU == 3050, "GRPC")) %>%
  # make index numeric as all values kept at this point are numeric
  mutate(index = as.numeric(index)) %>%
  # re-scale indices relative to a common base-year (1970) ("indstd")
  dplyr::select(AOU, species_code, scientific_name, common_name, breeding_biome, year, index,
    first_year_popest, last_year_popest, popest) %>%
  arrange(AOU, year) %>%
  group_by(AOU) %>%
  mutate(indstd = (index) / (first(index))) %>%
  ungroup()
  
# calculate scaling factor based on population size estimate (popest) from Rosenberg et al 2019
species_weights <- bbs_index %>%
  # for each species keep years between first and last years in popest file 
  dplyr::filter(year >= first_year_popest & year <= last_year_popest) %>%
  mutate(time_length = last_year_popest - first_year_popest + 1) %>%
  # divide popest by indstd from that year
  # popest represents the species' mean population size between first and last years
  mutate(species_weight = popest / indstd) %>%
  group_by(AOU, time_length) %>%
  # sum across years
  summarise(species_weight = sum(species_weight)) %>%
  ungroup() %>%
  # and divide by the number of years
  mutate(species_weight = species_weight / time_length) %>%
  dplyr::select(AOU, species_weight)

bbs_abundance <- left_join(bbs_index, species_weights, by = "AOU") %>%
  mutate(abundance = indstd * species_weight) %>%
  dplyr::select(AOU, species_code, scientific_name, common_name, breeding_biome, year, index,
    indstd, species_weight, abundance) %>%
  arrange(AOU, year) %>%
  group_by(AOU) %>%
  mutate(r = log(abundance / lag(abundance))) %>%
  ungroup()

save(bbs_abundance, file = "Library/bbs_abundance.rda")



#Species Trait Data:
#-------------------------------------------------------------------------------------------
##1. population estimates, Rosenberg et al. 2019
#    "popest"

##2. Niche classifications from Pigot et al. 2020,
#    "original.scientific" is the name used in the original dataset that was updated to match BBS
niche <- read.csv("RawData/pigot_niche.csv", sep=",", header=T, stringsAsFactors = F) 

##3. species biomass from Wilman et al. 2014
biomass <- read.csv("RawData/wilman_bodymass.csv", sep=",", header=T, stringsAsFactors = F)%>%
  select(Scientific.aou,BodyMass_Value)

##4. species life history from Myhrvold et al. 2015
#    "original_scientific" is the name used in the original dataset that was updated to match BBS
#    subset only reproductive/longevity traits from original Amniote dataset 
lhist<- read.csv("RawData/Amniote_subset.csv")



# combine with species abundance data
# include year replicates to calculate evenness and abundance across traits
traits<-left_join(bbs_abundance,
                  niche, by = c("scientific_name" = "scientific.name")) %>% 
  left_join(lhist, by = c("scientific_name" = "scientific.name")) %>% 
  mutate(original.scientific=case_when(
    original.scientific==""~scientific_name,
    original.scientific!=""~original.scientific))%>%
  left_join(biomass, by=c("original.scientific"="Scientific.aou")) 


#Classify species with missing Foraging Niches as Omnivores (Based on Trophic Niche)
traits$ForagingNiche[is.na(traits$ForagingNiche)]<-"Omnivore"

#Manually add missing trait records based on other species in the same genus and Cornell info
traits$longevity_y[traits$common_name=="Scaled Quail"]<-5.0 
traits$longevity_y[traits$common_name=="Lesser Nighthawk"]<-4.5
traits$longevity_y[traits$common_name=="Scissor-tailed Flycatcher"]<-12.0
traits$litters_or_clutches_per_y[traits$common_name=="Scissor-tailed Flycatcher"]<-1.5
traits$longevity_y[traits$common_name=="Couch's Kingbird"]<-9.0
traits$litters_or_clutches_per_y[traits$common_name=="Couch's Kingbird"]<-1
traits$longevity_y[traits$common_name=="Lawrence's Goldfinch"]<-10.0
traits$litters_or_clutches_per_y[traits$common_name=="Lawrence's Goldfinch"]<-1.5
traits$longevity_y[traits$common_name=="McCown's Longspur"]<-5.0
traits$longevity_y[traits$common_name=="Black-chinned Sparrow"]<-8.0
traits$litters_or_clutches_per_y[traits$common_name=="lack-chinned Sparrow"]<-1
traits$longevity_y[traits$common_name=="Cassin's Sparrow"]<-4.0
traits$litters_or_clutches_per_y[traits$common_name=="assin's Sparrow"]<-1.5
traits$longevity_y[traits$common_name=="Sprague's Pipit"]<-6.0
traits$litters_or_clutches_per_y[traits$common_name=="Sprague's Pipit"]<-2
traits$longevity_y[traits$common_name=="Rock Wren"]<-5.0
traits$longevity_y[traits$common_name=="Canyon Wren"]<-4.0
traits$litters_or_clutches_per_y[traits$common_name=="Canyon Wren"]<- 1.5
traits$longevity_y[traits$common_name=="Sedge Wren"]<-5.0
traits$longevity_y[traits$common_name=="Black-tailed Gnatcatcher"]<-6.0
traits$litters_or_clutches_per_y[traits$common_name=="Chestnut-collared Longspur"]<- 2.5
traits$litters_or_clutches_per_y[traits$common_name== "LeConte's Sparrow"]<-  1
traits$litters_or_clutches_per_y[traits$common_name=="Ladder-backed Woodpecker"]<-  1
traits$litters_or_clutches_per_y[traits$common_name=="Green Jay"]<- 1
traits$litters_or_clutches_per_y[traits$common_name=="Gilded Flicker"]<- 1
traits$litters_or_clutches_per_y[traits$common_name=="Pyrrhuloxia"]<- 1.5
traits$litters_or_clutches_per_y[traits$common_name=="Eastern Wood-Pewee"]<- 1

#replace additional missing data values with NA
traits<-naniar::replace_with_na(traits,replace = list(longevity_y = c(-999),
                                                            litter_or_clutch_size_n = c(-999),
                                                            litters_or_clutches_per_y = c(-999)))

#create annual reproduction trait (clutch size * clutches per year)
traits<-mutate(traits, reprod_rate=litter_or_clutch_size_n*litters_or_clutches_per_y)


save(traits, file = "Library/traits.rda")


