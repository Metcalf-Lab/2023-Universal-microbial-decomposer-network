library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(cowplot)
library(qiime2R)
library(reshape2)
library(nortest)
library(FSA)
library(ggalluvial)
library(readr)
library(forcats)
library(MuMIn)
library(lmerTest)

setwd("/Users/zacharyburcham/Dropbox/CSU/PMI_3_analyses/multi-omics_data/shotgun/Soil-modeling-Metcalf-lab-main/Bins/Jan2022_mapping/")

# soil data import
table = read.csv('ATP.csv', header=TRUE)
CMU_table = subset(table, facility=="CMU")
UTK_table = subset(table, facility=="UTK")
SHSU_table = subset(table, facility=="SHSU")

# clean up data
table_melt = melt(table, id.vars=c("SampleID","facility","Conv_time","host_subject_ID"), measure.vars = c("Carbohydrates","Lipids","AminoAcids"))
CMU_melt = melt(CMU_table, id.vars=c("SampleID","facility","Conv_time","host_subject_ID"), measure.vars = c("Carbohydrates","Lipids","AminoAcids"))
UTK_melt = melt(UTK_table, id.vars=c("SampleID","facility","Conv_time","host_subject_ID"), measure.vars = c("Carbohydrates","Lipids","AminoAcids"))
SHSU_melt = melt(SHSU_table, id.vars=c("SampleID","facility","Conv_time","host_subject_ID"), measure.vars = c("Carbohydrates","Lipids","AminoAcids"))

# Plot each separately
ggplot(subset(table_melt,variable=="AminoAcids"), aes(x=Conv_time, y=value, fill = facility, color = facility)) + 
  geom_point(size=1)+
  geom_smooth(method=lm)+
  theme_bw() +
  xlim(0,575) +
  ylim(0,0.67) +
  labs(title = "Amino Acid Metabolism Efficiency", y = "ATP per C-mol Substrate", x = "ADD")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5), 
        legend.title = element_blank(),
        legend.position = c(0.1, 0.91),
        legend.key.size = unit(0.4, 'cm'))+
  annotate("text", x = c(442,559), y = c(0.45,0.41), label = "***", size=8) +
  scale_fill_manual(values=c('ARF'='forestgreen', 'STAFS'='chartreuse2', 'FIRS'='mediumorchid4'))+
  scale_color_manual(values=c('ARF'='forestgreen', 'STAFS'='chartreuse2', 'FIRS'='mediumorchid4'))
ggsave("AA_ATPm.png", height=5, width=5, dpi=600, units="in")

ggplot(subset(table_melt,variable=="Carbohydrates"), aes(x=Conv_time, y=value, fill = facility, color = facility)) + 
  geom_point(size=1)+
  geom_smooth(method=lm)+
  theme_bw() +
  xlim(0,575) +
  ylim(0,0.67) +
  labs(title = "Carbohydrate Metabolism Efficiency", y = "ATP per C-mol Substrate", x = "ADD")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5), 
        legend.title = element_blank(),
        legend.position = c(0.1, 0.91),
        legend.key.size = unit(0.4, 'cm'))+
  annotate("text", x = c(559), y = c(0.31), label = "***", size=8) +
  scale_fill_manual(values=c('ARF'='forestgreen', 'STAFS'='chartreuse2', 'FIRS'='mediumorchid4'))+
  scale_color_manual(values=c('ARF'='forestgreen', 'STAFS'='chartreuse2', 'FIRS'='mediumorchid4'))
ggsave("Carb_ATPm.png", height=5, width=5, dpi=600, units="in")

ggplot(subset(table_melt,variable=="Lipids"), aes(x=Conv_time, y=value, fill = facility, color = facility)) + 
  geom_point(size=1)+
  geom_smooth(method=lm)+
  theme_bw() +
  xlim(0,575) +
  ylim(0,0.67) +
  labs(title = "Lipid Metabolism Efficiency", y = "ATP per C-mol Substrate", x = "ADD")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.title = element_text(hjust = 0.5), 
        legend.title = element_blank(),
        legend.position = c(0.1, 0.91),
        legend.key.size = unit(0.4, 'cm'))+
  annotate("text", x = c(559), y = c(0.2), label = "*", size=8) +
  scale_fill_manual(values=c('ARF'='forestgreen', 'STAFS'='chartreuse2', 'FIRS'='mediumorchid4'))+
  scale_color_manual(values=c('ARF'='forestgreen', 'STAFS'='chartreuse2', 'FIRS'='mediumorchid4'))
ggsave("Lipid_ATPm.png", height=5, width=5, dpi=600, units="in")


# AA linear models
FIRSlinearMod<-lmer((value)~Conv_time + (1|host_subject_ID), 
                    data=subset(CMU_melt, variable=='AminoAcids'), 
                    REML=TRUE)
FIRSsum = anova(FIRSlinearMod)
FIRSsum
FIRSR2m = r.squaredGLMM(FIRSlinearMod)[1]
FIRSR2m
ARFlinearMod<-lmer((value)~Conv_time + (1|host_subject_ID), 
                   data=subset(UTK_melt, variable=='AminoAcids'), 
                   REML=TRUE)
ARFsum = anova(ARFlinearMod)
ARFsum
ARFR2m = r.squaredGLMM(ARFlinearMod)[1]
ARFR2m
STAFSlinearMod<-lmer((value)~Conv_time + (1|host_subject_ID), 
                     data=subset(SHSU_melt, variable=='AminoAcids'), 
                     REML=TRUE)
STAFSsum = anova(STAFSlinearMod)
STAFSsum
STAFSR2m = r.squaredGLMM(STAFSlinearMod)[1]
STAFSR2m

# Carb linear models
FIRSlinearMod<-lmer((value)~Conv_time + (1|host_subject_ID), 
                    data=subset(CMU_melt, variable=='Carbohydrates'), 
                    REML=TRUE)
FIRSsum = anova(FIRSlinearMod)
FIRSsum
FIRSR2m = r.squaredGLMM(FIRSlinearMod)[1]
FIRSR2m
ARFlinearMod<-lmer((value)~Conv_time + (1|host_subject_ID), 
                   data=subset(UTK_melt, variable=='Carbohydrates'), 
                   REML=TRUE)
ARFsum = anova(ARFlinearMod)
ARFsum
ARFR2m = r.squaredGLMM(ARFlinearMod)[1]
ARFR2m
STAFSlinearMod<-lmer((value)~Conv_time + (1|host_subject_ID), 
                     data=subset(SHSU_melt, variable=='Carbohydrates'), 
                     REML=TRUE)
STAFSsum = anova(STAFSlinearMod)
STAFSsum
STAFSR2m = r.squaredGLMM(STAFSlinearMod)[1]
STAFSR2m

# Lipids linear models
FIRSlinearMod<-lmer((value)~Conv_time + (1|host_subject_ID), 
                    data=subset(CMU_melt, variable=='Lipids'), 
                    REML=TRUE)
FIRSsum = anova(FIRSlinearMod)
FIRSsum
FIRSR2m = r.squaredGLMM(FIRSlinearMod)[1]
FIRSR2m
ARFlinearMod<-lmer((value)~Conv_time + (1|host_subject_ID), 
                   data=subset(UTK_melt, variable=='Lipids'), 
                   REML=TRUE)
ARFsum = anova(ARFlinearMod)
ARFsum
ARFR2m = r.squaredGLMM(ARFlinearMod)[1]
ARFR2m
STAFSlinearMod<-lmer((value)~Conv_time + (1|host_subject_ID), 
                     data=subset(SHSU_melt, variable=='Lipids'), 
                     REML=TRUE)
STAFSsum = anova(STAFSlinearMod)
STAFSsum
STAFSR2m = r.squaredGLMM(STAFSlinearMod)[1]
STAFSR2m
