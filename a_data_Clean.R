## Prepares the data for constructing matrix population models
## The following packages are needed. If not installed, they need to be installed.


#devtools::install_github("james-thorson/FishLife")
#install.packages("rfishbase")
#install.package("tidyverse")

library(tidyverse)

rm(list=ls()) 

load("species.Rdata") # List of species [[1]] Common Names, [[2]] Scientific Names

NSP <- lengths(SP)[1] # The number of species

## Separate genus and species to prepare for FishLife Packaage
genus<-unlist(map(strsplit(SP$scientificName," "), ~.x[1]))
species<-unlist(map(strsplit(SP$scientificName," "), ~.x[2]))

species[species =="spp."] <- NA  # Eliminate the ones without species

INDEX<-which(!is.na(species))

#########################################
#### FishLife Package (Thorson 2019) ####

## Run just once and save. 
## Recommend checking the scientific names (many have multiple names)

# PREDICTION<-map( INDEX , ~append(list(SciName=SP$scientificName[.x], commonName = SP$commonName[.x], genus=genus[.x],species=species[.x]), FishLife::Plot_taxa( FishLife::Search_species(Genus=genus[.x],Species=species[.x])$match_taxonomy, mfrow=c(2,2) )[[1]]))
# graphics.off()
# save(PREDICTION, file="FishLifeTraits.Rdata")

load("FishLifeTraits.Rdata")

#########################################
#### FishBase  W = a_lw * L ^b_lw  (length is in cm and weight in gram of wet weight)

## list of some of functions in FishBase (these were looked at for potential uses)
## c("species","estimate","fecundity","ecology","maturity","fooditems","reproduction","length_weight")

## Obtain allometric data from Fishbase and Save -- This should be run only once. 

# TRAITS <- rfishbase::length_weight(SP$scientificName) 
# LL<- rfishbase::length_length(SP$scientificName) 
# save(LL,TRAITS, file="FishBaseTraits.Rdata")

load("FishBaseTraits.Rdata") # Load the information rather than obtaining the data from fishBase each time

# Calculate the mean coefficients for length-weight relationships 
LW <- TRAITS %>%
  group_by (Species,Type) %>%
  summarise (a=mean(a, rm.na=TRUE),at=mean(aTL,rm.na=TRUE),b=mean(b, rm.na=TRUE))

ALLOM<-data.frame(Species=unique(TRAITS$Species))

## Species with the allometric relationship based on the total length
LW_T<-LW %>% filter(Type == "TL") %>% select(Species, Type, a, b)

## Species with the allometric relationship based on the standard length but not total length
LW_S<-LW %>% filter(Type == "SL") %>% 
  mutate(a=if_else(is.na(at),a,at),Type=if_else(is.na(at),"SL","TL")) %>% 
  select(Species, Type, a, b) %>% 
  filter(Species %in% setdiff(Species,LW_T$Species))

## Species with the allometric relationship based on the fork length but not total or standard length
LW_F<-LW %>% filter(Type == "FL") %>% 
  mutate(a=if_else(is.na(at),a,at),Type=if_else(is.na(at),"FL","TL")) %>% 
  select(Species, Type, a, b) %>%
  filter(Species %in% setdiff(Species,c(LW_T$Species,LW_S$Species)))

## Combine all of the allometric relationship data
LW <- rbind(rbind(LW_T,LW_S),LW_F) %>%
  right_join(ALLOM)

LW_T <- LW %>% filter(Type == "TL") %>% mutate(alpha = 0, beta= 1)
LW_S <- LW %>% filter(Type == "SL") %>% mutate(alpha =NA, beta =NA)
LW_F <- LW %>% filter(Type == "FL") %>% mutate(alpha =NA, beta =NA)
LW_N <- LW %>% filter(is.na(Type))  %>% mutate(alpha =NA, beta =NA)

## Get the length-length parameters for the L-W model with STANDARD lengths
for (k in 1:nrow(LW_S)){
  PARAM <- LL %>% 
    filter(Species == as.character(LW_S[k,"Species"])) %>%
    filter (Length2 %in% c("TL","SL")) %>%
    filter (Length1 %in% c("TL","SL"))
  if (nrow(PARAM)>0){
  if (PARAM$Length2 == "TL"){
    LW_S[k,"alpha"] <- -as.numeric(PARAM$a)/as.numeric(PARAM$b)
    LW_S[k,"beta"] <- 1/as.numeric(PARAM$b)
  } else {
    LW_S[k,"alpha"] <- as.numeric(PARAM$a)
    LW_S[k,"beta"] <- as.numeric(PARAM$b)
  }} else
  {LW_N <- rbind(LW_N,LW_S[k,])
  LW_S<- LW_S[-k,]}
}

## Get the length-length parameters for the L-W model with FORK lengths
for (k in 1:nrow(LW_F)){
  PARAM <- LL %>% 
    filter(Species == as.character(LW_F[k,"Species"])) %>%
    filter (Length2 %in% c("TL","FL")) %>%
    filter (Length1 %in% c("TL","FL"))
  if (nrow(PARAM)>0){
    if (PARAM$Length2 == "TL"){
      LW_F[k,"alpha"] <- -as.numeric(PARAM$a)/as.numeric(PARAM$b)
      LW_F[k,"beta"] <- 1/as.numeric(PARAM$b)
    } else {
      LW_F[k,"alpha"] <- as.numeric(PARAM$a)
      LW_F[k,"beta"] <- as.numeric(PARAM$b)
    }} else
    {LW_N <- rbind(LW_N,LW_F[k,])
    LW_F<- LW_F[-k,]}
}

## Eliminated the ones without allometric relationship info
LW <- rbind(LW_T,LW_S)
LW <- rbind(LW,LW_F)


#########################################
### Lorenzen (1996): 90% CI is given in the paper M <- M_u * W ^ b_L (W is wet weight in gram)
M_u   <- 3.00     # Mortality at 1 g
B_L   <- -0.288


#########################################
### Scientific names in FishLife package
sci_name_fishlife<- as.data.frame(unlist(map(PREDICTION,~.x$SciName)))
names(sci_name_fishlife)<-"SciName"

## Save for the analysis
save(sci_name_fishlife, M_u, B_L,LW,PREDICTION, file="MPM_DATA.Rdata")