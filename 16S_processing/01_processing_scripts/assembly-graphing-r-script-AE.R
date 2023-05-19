##### Libraries #####
library(tidyverse)
library(dplyr)
library(readr)
library(readxl)
library(ggplot2)
library(ggthemes)
library(viridis)
library(forcats)
library(ggcorrplot)
library(stringr)
library(tidyr)
library(RColorBrewer)
library(patchwork)
library(qiime2R)
library(phyloseq)
library(grid)
library(gridExtra)
library(ggtext)
library(vegan)
library(reshape2)
library(ggpubr)
library(vegan)
library(FSA)
library(nortest)

setwd('~/Dropbox/PMI_3_analyses/multi-omics_data/amplicon/16S/')

## Load Data
cmu_control_early <- read.csv('deterministic_outputs_comp-to-control/hip_only/ControlvEarly_CMU_hip_S3_Community_bNTI_TotalCounts.csv',
                        row.names = 1)
cmu_control_active <- read.csv('deterministic_outputs_comp-to-control/hip_only/ControlvActive_CMU_hip_S3_Community_bNTI_TotalCounts.csv',
                               row.names = 1)
cmu_control_advanced <- read.csv('deterministic_outputs_comp-to-control/hip_only/ControlvAdvanced_CMU_hip_S3_Community_bNTI_TotalCounts.csv',
                               row.names = 1)

shsu_control_early <- read.csv('deterministic_outputs_comp-to-control/hip_only/ControlvEarly_SHSU_hip_S3_Community_bNTI_TotalCounts.csv',
                              row.names = 1)
shsu_control_active <- read.csv('deterministic_outputs_comp-to-control/hip_only/ControlvActive_SHSU_hip_S3_Community_bNTI_TotalCounts.csv',
                               row.names = 1)
shsu_control_advanced <- read.csv('deterministic_outputs_comp-to-control/hip_only/ControlvAdvanced_SHSU_hip_S3_Community_bNTI_TotalCounts.csv',
                                 row.names = 1)

utk_control_early <- read.csv('deterministic_outputs_comp-to-control/hip_only/ControlvEarly_UTK_hip_S3_Community_bNTI_TotalCounts.csv',
                               row.names = 1)
utk_control_active <- read.csv('deterministic_outputs_comp-to-control/hip_only/ControlvActive_UTK_hip_S3_Community_bNTI_TotalCounts.csv',
                                row.names = 1)
utk_control_advanced <- read.csv('deterministic_outputs_comp-to-control/hip_only/ControlvAdvanced_UTK_hip_S3_Community_bNTI_TotalCounts.csv',
                                  row.names = 1)

# Make matrices and melt
cmu_control_early <- melt(as.matrix(cmu_control_early))
cmu_control_active <- melt(as.matrix(cmu_control_active))
cmu_control_advanced <- melt(as.matrix(cmu_control_advanced))

shsu_control_early <- melt(as.matrix(shsu_control_early))
shsu_control_active <- melt(as.matrix(shsu_control_active))
shsu_control_advanced <- melt(as.matrix(shsu_control_advanced))

utk_control_early <- melt(as.matrix(utk_control_early))
utk_control_active <- melt(as.matrix(utk_control_active))
utk_control_advanced <- melt(as.matrix(utk_control_advanced))

## Make metadata
metadata = read_tsv('03_metadata/combined-metadata-simple-nov2020_R.txt') %>% 
  dplyr::rename('Samples' = sample) %>% mutate_if(is.factor,as.character)

# remove self-comparisons and NAs
m_cmu_early = cmu_control_early %>%
  filter(as.character(Var1) != as.character(Var2)) %>% filter(!is.na(value))%>%
  mutate_if(is.factor,as.character)
m_cmu_active = cmu_control_active %>%
  filter(as.character(Var1) != as.character(Var2)) %>% filter(!is.na(value))%>%
  mutate_if(is.factor,as.character) 
m_cmu_advanced = cmu_control_advanced %>% 
  filter(as.character(Var1) != as.character(Var2)) %>%  filter(!is.na(value))%>%
  mutate_if(is.factor,as.character)

m_shsu_early = shsu_control_early %>%
  filter(as.character(Var1) != as.character(Var2)) %>% filter(!is.na(value))%>%
  mutate_if(is.factor,as.character)
m_shsu_active = shsu_control_active %>%
  filter(as.character(Var1) != as.character(Var2)) %>% filter(!is.na(value))%>%
  mutate_if(is.factor,as.character) 
m_shsu_advanced = shsu_control_advanced %>% 
  filter(as.character(Var1) != as.character(Var2)) %>%  filter(!is.na(value))%>%
  mutate_if(is.factor,as.character)

m_utk_early = utk_control_early %>%
  filter(as.character(Var1) != as.character(Var2)) %>% filter(!is.na(value))%>%
  mutate_if(is.factor,as.character)
m_utk_active = utk_control_active %>%
  filter(as.character(Var1) != as.character(Var2)) %>% filter(!is.na(value))%>%
  mutate_if(is.factor,as.character)
m_utk_advanced = utk_control_advanced %>% 
  filter(as.character(Var1) != as.character(Var2)) %>%  filter(!is.na(value))%>%
  mutate_if(is.factor,as.character)


#bind dfs
combined_cmu<-rbind(m_cmu_early,m_cmu_active,m_cmu_advanced)

combined_shsu<-rbind(m_shsu_early,m_shsu_active,m_shsu_advanced)

combined_utk<-rbind(m_utk_early,m_utk_active,m_utk_advanced)

#new metadata for filtering
sd = metadata %>%
  select(Samples,soil_control,decomp_stage,facility,season) %>%
  mutate_if(is.factor,as.character)

# combined distances with sample data
colnames(sd) = c("Var1", "soil_control","decomp_stage","facility","season")
cmu_c_sd = left_join(combined_cmu, sd, by = "Var1")
colnames(sd) = c("Var2", "soil_control2","decomp_stage2","facility2","season2")
cmu_c_sd= left_join(cmu_c_sd, sd, by = "Var2")

colnames(sd) = c("Var1", "soil_control","decomp_stage","facility","season")
shsu_c_sd = left_join(combined_shsu, sd, by = "Var1")
colnames(sd) = c("Var2", "soil_control2","decomp_stage2","facility2","season2")
shsu_c_sd= left_join(shsu_c_sd, sd, by = "Var2")

colnames(sd) = c("Var1", "soil_control","decomp_stage","facility","season")
utk_c_sd = left_join(combined_utk, sd, by = "Var1")
colnames(sd) = c("Var2", "soil_control2","decomp_stage2","facility2","season2")
utk_c_sd= left_join(utk_c_sd, sd, by = "Var2")

combined_all_orig<-rbind(cmu_c_sd,shsu_c_sd,utk_c_sd)

#filter so that only soil samples and control pairwise comparisons remain 
combined_all<-combined_all_orig %>% filter(soil_control != soil_control2) 
#labels need to be the correct decomp stage for plotting regardless of which column (Var1 or Var2) the control is in
combined_all$decomp_stage<-ifelse(combined_all$decomp_stage=="soil_control",combined_all$decomp_stage2,combined_all$decomp_stage)
combined_all$decomp_stage2<-ifelse(combined_all$decomp_stage2=="soil_control",combined_all$decomp_stage,combined_all$decomp_stage2)

#filter so that we are only looking at within stage comparisons
combined_all_within<-combined_all_orig %>% 
  filter(decomp_stage == decomp_stage2) %>%
  filter(soil_control == soil_control2)
#labels need to be the correct decomp stage for plotting regardless of which column (Var1 or Var2) the control is in
combined_all_within$decomp_stage<-ifelse(combined_all_within$decomp_stage=="soil_control",combined_all_within$decomp_stage2,combined_all_within$decomp_stage)
combined_all_within$decomp_stage2<-ifelse(combined_all_within$decomp_stage2=="soil_control",combined_all_within$decomp_stage,combined_all_within$decomp_stage2)


# Plot of absolute value BNTI
#between stages (comparison of control vs stage)
combined_all %>% 
  ggplot(aes(x = fct_relevel(decomp_stage, 'early', 'active', 'advanced'), y = value)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = Inf, 
           fill = "#FCAE12", alpha = .3, color = NA) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -2, 
           fill = "cornflowerblue", alpha = .3, color = NA) +
  geom_violin(aes(fill=facility)) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  stat_summary(fun=mean, geom="point", shape=5, size=1, color="black", position = position_dodge(0.5)) +
  facet_grid(. ~ facility, scales = 'free') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.text.x = element_text(size=7, color = "black", angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(size=7, color = "black", angle = 0, hjust = 0.5, vjust = 0.5),
        axis.title.x = element_text(size=8, color = "black"),
        axis.title.y = element_text( size=8, color = "black"),
        plot.title = element_text(hjust = 0.5, color = "black", size = 10),
        legend.position = "none", 
        plot.background = element_rect(fill='white')) +
  labs(x = 'Decomposition Stage',
       y = 'bNTI Value') +
  geom_hline(yintercept = 2, linetype = 2, alpha=0.75) +
  geom_hline(yintercept = -2, linetype = 2, alpha=0.75) +
  scale_fill_manual(name='Facility', values=c('CMU'='mediumorchid4', 'SHSU'='royalblue3', 'UTK'='mediumseagreen'))
ggsave("../../figures/abs_bnti_facility_controlvdecomp.png",dpi=600, width=5, height=4, units="in")

# Statistical comparison
# CMU
CMUdata_all=subset(combined_all, facility=="CMU")
ad.test(CMUdata_all$value) # not normal
kruskal.test((CMUdata_all$value) ~ CMUdata_all$decomp_stage)
dunnTest((CMUdata_all$value) ~ CMUdata_all$decomp_stage, method='bh')

# SHSU
SHSUdata_all=subset(combined_all, facility=="SHSU")
ad.test(SHSUdata_all$value) # not normal
kruskal.test((SHSUdata_all$value) ~ SHSUdata_all$decomp_stage)
dunnTest((SHSUdata_all$value) ~ SHSUdata_all$decomp_stage, method='bh')

# UTK
UTKdata_all=subset(combined_all, facility=="UTK")
ad.test(UTKdata_all$value) # not normal
kruskal.test((UTKdata_all$value) ~ UTKdata_all$decomp_stage)
dunnTest((UTKdata_all$value) ~ UTKdata_all$decomp_stage, method='bh')

#####within stages
data.segm1<-data.frame(x=0.8,y=3,xend=1.2,yend=3, facility="CMU")
data.segm2<-data.frame(x=1.8,y=5,xend=2.2,yend=5, facility="CMU")
data.segm3<-data.frame(x=2.8,y=5.7,xend=3.2,yend=5.7, facility="CMU")
data.segm4<-data.frame(x=0.8,y=5.2,xend=1.2,yend=5.2, facility="SHSU")
data.segm5<-data.frame(x=1.8,y=6.6,xend=2.2,yend=6.6, facility="SHSU")
data.segm6<-data.frame(x=2.8,y=6,xend=3.2,yend=6, facility="SHSU")
data.segm7<-data.frame(x=0.8,y=3,xend=1.2,yend=3, facility="UTK")
data.segm8<-data.frame(x=1.8,y=6.6,xend=2.2,yend=6.6, facility="UTK")
data.segm9<-data.frame(x=2.8,y=4.75,xend=3.2,yend=4.75, facility="UTK")

ann_text1 <- data.frame(x=1,y=3,facility="CMU",label="*")
ann_text2 <- data.frame(x=2,y=5.5,facility="CMU",label="ns")
ann_text3 <- data.frame(x=3,y=5.7,facility="CMU",label="**")
ann_text4 <- data.frame(x=1,y=5.2,facility="SHSU",label="***")
ann_text5 <- data.frame(x=2,y=6.6,facility="SHSU",label="***")
ann_text6 <- data.frame(x=3,y=6,facility="SHSU",label="***")
ann_text7 <- data.frame(x=1,y=3,facility="UTK",label="***")
ann_text8 <- data.frame(x=2,y=6.6,facility="UTK",label="***")
ann_text9 <- data.frame(x=3,y=4.75,facility="UTK",label="***")


combined_all_within %>% 
  ggplot(aes(x = fct_relevel(decomp_stage, 'early', 'active', 'advanced'), y = value, fill=soil_control)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = Inf, 
           fill = "#FCAE12", alpha = .3, color = NA) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -2, 
           fill = "cornflowerblue", alpha = .3, color = NA) +
  #geom_boxplot(width = 0.5, outlier.shape =  21, outlier.size = 0.5, outlier.color = NA, position = position_dodge(width=0.55)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, position = position_dodge(width=0.55)) +
  stat_summary(fun=mean, geom="point", shape=5, size=1.75, color="black", position = position_dodge(0.55)) +
  facet_grid(. ~ facility, scales = 'free') +
  ylim(-6,10) +
  geom_segment(data=data.segm1, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE)+
  geom_segment(data=data.segm2, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE)+
  geom_segment(data=data.segm3, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE)+
  geom_segment(data=data.segm4, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE)+
  geom_segment(data=data.segm5, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE)+
  geom_segment(data=data.segm6, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE)+
  geom_segment(data=data.segm7, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE)+
  geom_segment(data=data.segm8, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE)+
  geom_segment(data=data.segm9, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE)+
  geom_text(data=ann_text1,aes(x=x,y=y,label=label,size=1),inherit.aes=FALSE,show.legend = FALSE)+
  geom_text(data=ann_text2,aes(x=x,y=y,label=label,size=1),inherit.aes=FALSE,show.legend = FALSE)+
  geom_text(data=ann_text3,aes(x=x,y=y,label=label,size=1),inherit.aes=FALSE,show.legend = FALSE)+
  geom_text(data=ann_text4,aes(x=x,y=y,label=label,size=1),inherit.aes=FALSE,show.legend = FALSE)+
  geom_text(data=ann_text5,aes(x=x,y=y,label=label,size=1),inherit.aes=FALSE,show.legend = FALSE)+
  geom_text(data=ann_text6,aes(x=x,y=y,label=label,size=1),inherit.aes=FALSE,show.legend = FALSE)+
  geom_text(data=ann_text7,aes(x=x,y=y,label=label,size=1),inherit.aes=FALSE,show.legend = FALSE)+
  geom_text(data=ann_text8,aes(x=x,y=y,label=label,size=1),inherit.aes=FALSE,show.legend = FALSE)+
  geom_text(data=ann_text9,aes(x=x,y=y,label=label,size=1),inherit.aes=FALSE,show.legend = FALSE)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.text.x = element_text(size=7, color = "black", angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(size=7, color = "black", angle = 0, hjust = 0.5, vjust = 0.5),
        axis.title.x = element_text(size=8, color = "black"),
        axis.title.y = element_text( size=8, color = "black"),
        plot.title = element_text(hjust = 0.5, color = "black", size = 10),
        legend.text = element_text(size=6),
        legend.title = element_blank(),
        legend.key.width = unit(0.15, "in"),
        legend.key.height = unit(0.15, "in"),
        legend.key = element_blank(),
        legend.position = c(0.11,0.96),
        legend.direction = "vertical",
        legend.background = element_rect(fill=NA,color=NA, size=0.25), 
        plot.background = element_rect(fill='white')) +
  labs(x = 'Decomposition Stage',
       y = 'bNTI Value') +
  geom_hline(yintercept = 2, linetype = 2, alpha=0.75) +
  geom_hline(yintercept = -2, linetype = 2, alpha=0.75) +
  scale_fill_manual(labels = c("Decomposition Soil","Control Soil"), values=c('y'='gray80', 'n'='gray40'))
ggsave("../../figures/bnti_within_facility.png",dpi=600, width=5, height=4, units="in")

# Statistical comparison
# CMU
CMUdata=subset(combined_all_within, facility=="CMU")
res.aov2 = aov(value ~ decomp_stage * soil_control, data=CMUdata)
summary(res.aov2)
TukeyHSD(res.aov2)

CMUearlydata=subset(CMUdata, decomp_stage=="early")
wilcox.test(value ~ soil_control, data=CMUearlydata)

CMUactivedata=subset(CMUdata, decomp_stage=="active")
wilcox.test(value ~ soil_control, data=CMUactivedata)

CMUadvanceddata=subset(CMUdata, decomp_stage=="advanced")
wilcox.test(value ~ soil_control, data=CMUadvanceddata)

# SHSU
SHSUdata=subset(combined_all_within, facility=="SHSU")
res.aov2 = aov(value ~ decomp_stage * soil_control, data=SHSUdata)
summary(res.aov2)
TukeyHSD(res.aov2)

SHSUearlydata=subset(SHSUdata, decomp_stage=="early")
wilcox.test(value ~ soil_control, data=SHSUearlydata)

SHSUactivedata=subset(SHSUdata, decomp_stage=="active")
wilcox.test(value ~ soil_control, data=SHSUactivedata)

SHSUadvanceddata=subset(SHSUdata, decomp_stage=="advanced")
wilcox.test(value ~ soil_control, data=SHSUadvanceddata)

# UTK
UTKdata=subset(combined_all_within, facility=="UTK")
res.aov2 = aov(value ~ decomp_stage * soil_control, data=UTKdata)
summary(res.aov2)
TukeyHSD(res.aov2)

UTKearlydata=subset(UTKdata, decomp_stage=="early")
wilcox.test(value ~ soil_control, data=UTKearlydata)

UTKactivedata=subset(UTKdata, decomp_stage=="active")
wilcox.test(value ~ soil_control, data=UTKactivedata)

UTKadvanceddata=subset(UTKdata, decomp_stage=="advanced")
wilcox.test(value ~ soil_control, data=UTKadvanceddata)







# Just deterministic
# Barplot of percenteages relative abundance
categorical <- combined_all %>% 
  mutate(deterministic = case_when(
           value >= 2 ~ 'Heterogenous\nSelection',
           value <= -2 ~ 'Homogenous\nSelection',
           value < 2 & value > -2 ~ 'Stochastic')) %>% 
  filter(deterministic != 'NA')
 
categorical %>% 
  group_by(decomp_stage, facility, deterministic) %>% 
  dplyr::summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  ungroup() %>% 
  ggplot(aes(x = fct_relevel(decomp_stage, 'early', 'active', 'advanced'), y = freq, fill = deterministic)) +
  geom_bar(stat = 'identity',  alpha = 0.8) +
  facet_grid(. ~ facility, scales = 'free') +
  scale_fill_manual(values = c("#FCAE12", "#380F63", '#B93556')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.text.x = element_text(size=7, color = "black", angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(size=7, color = "black", angle = 0, hjust = 0.5, vjust = 0.5),
        axis.title.x = element_text(size=8, color = "black"),
        axis.title.y = element_text( size=8, color = "black"),
        plot.title = element_text(hjust = 0.5, color = "black", size = 10),
        legend.title.align = 0, 
        plot.background = element_rect(fill='white'),
        legend.text=element_text(size=7, color = "black"),
        legend.title =element_text(size=8, color = "black"),
        legend.background = element_rect(fill = "white")) +
  labs(x = 'Decomposition Stage',
       y = 'Relative Abundance',
       fill = 'Assembly Force') +
  scale_y_continuous(labels=scales::percent)
ggsave("../../figures/assembly_force_facility.png",dpi=600, width=5, height=4, units="in")

# Barplot of percenteages relative abundance within stage
categorical <- combined_all_within %>% 
  mutate(deterministic = case_when(
    value >= 2 ~ 'Heterogenous\nSelection',
    value <= -2 ~ 'Homogenous\nSelection',
    value < 2 & value > -2 ~ 'Stochastic')) %>% 
  filter(deterministic != 'NA') %>%
  filter(decomp_stage != 'control')

categorical %>% 
  group_by(decomp_stage, facility, deterministic) %>% 
  dplyr::summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  ungroup() %>% 
  ggplot(aes(x = fct_relevel(decomp_stage, 'early', 'active', 'advanced'), y = freq, fill = deterministic)) +
  geom_bar(stat = 'identity',  alpha = 0.8) +
  facet_grid(. ~ facility, scales = 'free') +
  scale_fill_manual(values = c("#FCAE12", "#380F63", '#B93556')) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.text.x = element_text(size=7, color = "black", angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(size=7, color = "black", angle = 0, hjust = 0.5, vjust = 0.5),
        axis.title.x = element_text(size=8, color = "black"),
        axis.title.y = element_text( size=8, color = "black"),
        plot.title = element_text(hjust = 0.5, color = "black", size = 10),
        legend.title.align = 0, 
        plot.background = element_rect(fill='white'),
        legend.text=element_text(size=7, color = "black"),
        legend.title =element_text(size=8, color = "black"),
        legend.background = element_rect(fill = "white")) +
  labs(x = 'Decomposition Stage',
       y = 'Relative Abundance',
       fill = 'Assembly Force') +
  scale_y_continuous(labels=scales::percent)
ggsave("../../figures/assembly_force_within_facility.png",dpi=600, width=5, height=4, units="in")
