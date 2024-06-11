## This code is to plot the figures
## After running "c_variabilityPAR.R", please restart R. 
## There appears to be a conflict among packages. 

library(tidyverse)
library(scales)

rm(list=ls())

load(file="Demo_Results_Variable_V2.Rdata")
load("Demo_Results.Rdata")

FINAL_RESULTS_DF <- bind_rows(FINAL_RESULTS) 

MEDIANS<- FINAL_RESULTS_DF %>%
  filter(!is.na(GEN))%>%
  group_by(commonName) %>%
  summarize(median_value = median (GEN))


figure1 <- FINAL_RESULTS_DF %>%
  left_join(MEDIANS,by = "commonName") %>%
  mutate(commonName = factor(commonName, levels = MEDIANS$commonName[order(log10(MEDIANS$median_value))])) %>%
  filter (GEN>10^-2) %>%
  ggplot(aes(y=GEN, x= commonName)) +
  geom_boxplot() +
  theme_bw(base_size = 16)+
  scale_y_continuous(trans="log10") +
  coord_flip() + xlab("Common Name") + ylab("Generation Time")
print(figure1)

MEDIANS<- FINAL_RESULTS_DF %>%
  filter(!is.na(DPR))%>%
  group_by(commonName) %>%
  summarize(median_value = median (DPR))

figure2 <- FINAL_RESULTS_DF %>%
  left_join(MEDIANS,by = "commonName") %>%
  mutate(commonName = factor(commonName, levels = MEDIANS$commonName[order(log10(MEDIANS$median_value))])) %>%
  filter (DPR>0.5) %>%
  ggplot(aes(y=DPR, x= commonName)) +
  geom_boxplot() +
  theme_bw(base_size = 16)+
  scale_y_continuous(trans="log10") +
  coord_flip() + xlab("Common Name") + ylab("Damping Ratio")
print(figure2)


MEDIANS<- FINAL_RESULTS_DF %>%
  filter(!is.na(EXP_RES))%>%
  group_by(commonName) %>%
  summarize(median_value = median (EXP_RES))

figure3 <- FINAL_RESULTS_DF %>%
  left_join(MEDIANS,by = "commonName") %>%
  mutate(commonName = factor(commonName, levels = MEDIANS$commonName[order(log10(MEDIANS$median_value))])) %>%
  #filter (EXP_RES>0.01) %>%
  ggplot(aes(y=EXP_RES, x= commonName)) +
  geom_boxplot() +
  theme_bw(base_size = 16)+
  scale_y_continuous(trans="log10") +
  coord_flip() + xlab("Common Name") + ylab("Resilience")
print(figure3)

PAIRS_DATA <- FINAL_RESULTS_DF %>%
  filter(!is.na(GEN))%>%
  group_by(commonName)%>%
  summarize (Generation_TIME = median(GEN), Damping_Ratio = median(DPR), Resilience = median(EXP_RES))%>%
  select(-commonName)%>%
  pairs(col="blue")

INDICES <- c("(a)", "(b)","(c)","(d)","(e)", "(f)","(g)","(h)","(i)")
commonName<-sapply(dem_results, function(x)x[["commonName"]] )
SP <- data.frame(commonName=commonName[c(4,13,14,18,19,24,21,29,23)]) %>%
  mutate(SP=paste(INDICES,commonName,sep=" "))  


SENS<-function (LI){
  INDICES=which(LI$S==LI$S,arr.ind = TRUE)
  LI$S <- LI$S * (LI$M>0)
  MATRIX=data.frame(sen=as.vector(LI$S),commonName=LI$commonName, row=INDICES[,1],col=INDICES[,2])
}

figure4<-bind_rows(lapply(dem_results,SENS)) %>%
  filter(commonName %in% SP$commonName) %>%
  left_join(SP) %>%
  ggplot(aes(x=col,y=row,fill=sen))+
  geom_tile()+
  theme_classic(base_size = 10)+
  theme (strip.background = element_blank (),
         strip.text = element_text(size =10, hjust = 0))+
    scale_fill_gradientn(colours = c("white", "green", "red"), values = c(0,0.1,1),name="Sensitivity")+
  scale_x_continuous(breaks=scales::pretty_breaks())+scale_y_reverse(breaks=scales::pretty_breaks())+
  ylab("Age class to transition into")+xlab("Age class to transition from")+
  facet_wrap(~SP,scales ="free")
print(figure4)



ELAS<-function (LI){
  INDICES=which(LI$E==LI$E,arr.ind = TRUE)
  MATRIX=data.frame(elas=as.vector(LI$E),commonName=LI$commonName, row=INDICES[,1],col=INDICES[,2])
}

figure5<-bind_rows(lapply(dem_results,ELAS)) %>%
  filter(commonName %in% SP$commonName) %>%
  left_join(SP) %>%
  ggplot(aes(x=col,y=row,fill=elas))+
  geom_tile()+
  theme_classic(base_size = 10)+
  theme (strip.background = element_blank (),
         strip.text = element_text(size =10, hjust = 0))+
  scale_fill_gradientn(colours = c("white", "green", "red"), values = c(0,0.1,1),name="Elasticity")+
  scale_x_continuous(breaks=scales::pretty_breaks())+scale_y_reverse(breaks=scales::pretty_breaks())+
  ylab("Age class to transition into")+xlab("Age class to transition from")+
   facet_wrap(~SP,scales ="free")
print(figure5)

SSD <-function (LI){
  INDICES=which(LI$w==LI$w)
  MATRIX=data.frame(ssd=as.vector(LI$w),commonName=LI$commonName, row=INDICES)
}

figure6<-bind_rows(lapply(dem_results,SSD)) %>%
  filter(commonName %in% SP$commonName) %>%
  left_join(SP) %>%
  ggplot(aes(x=row,y=ssd)) +
  geom_col()+
  theme_classic(base_size = 10)+
  theme (strip.background = element_blank (),
         strip.text = element_text(size =10, hjust = 0))+
  ylab("Stable Stage Distribution")+ 
  scale_x_continuous(breaks=scales::pretty_breaks())+
  xlab ("Age Class")+
  facet_wrap(~SP,scales ="free")
print(figure6)


REP <-function (LI){
  INDICES=which(LI$v==LI$v)
  MATRIX=data.frame(rep=as.vector(LI$v),commonName=LI$commonName, row=INDICES)
}

figure7<-bind_rows(lapply(dem_results,REP)) %>%
  filter(commonName %in% SP$commonName) %>%
  left_join(SP) %>%
  ggplot(aes(x=row,y=rep)) +
  geom_col()+
  theme_classic(base_size = 10)+
  theme (strip.background = element_blank (),
         strip.text = element_text(size =10, hjust = 0))+
  ylab("Reproductive Value")+ 
  scale_x_continuous(breaks=pretty_breaks())+
  xlab ("Age Class")+
  facet_wrap(~SP,scales ="free")
print(figure7)

## Export figures. Exporting one at a time manually is recommended (does not seem to work if all are exported in the same code)

# tiff("GTime.tiff", units="in", width=10, height=7, res=600, compression="jpeg")
# figure1
# dev.off()

# tiff("DPRatio.tiff", units="in", width=10, height=7, res=600, compression="jpeg")
# figure2
# dev.off()

# tiff("Resilience.tiff", units="in", width=10, height=7, res=600, compression="jpeg")
# figure3
# dev.off()

# tiff("SENS.tiff", units="in", width=10, height=7, res=600, compression="jpeg")
# figure4
# dev.off()

# tiff("ELAS.tiff", units="in", width=10, height=7, res=600, compression="jpeg")
# figure5
# dev.off()

# tiff("SSD.tiff", units="in", width=10, height=7, res=600, compression="jpeg")
# figure6
# dev.off()

# tiff("REPRO.tiff", units="in", width=10, height=7, res=600, compression="jpeg")
# figure7
# dev.off()