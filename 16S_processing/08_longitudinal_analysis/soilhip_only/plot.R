
library(ggplot2)
library(tidyverse)
library(scales) 
library(qiime2R)
library(cowplot)
library(Rmisc)
library(dplyr)
library(viridis)
library(ggpmisc)
library(ggsignif)
library(FSA)
library(hrbrthemes)
library(car)
library(lsmeans)
library(reshape2)
library(multcomp)
library(lme4)
library(lmerTest)
library(merTools)
library(redres)
library(ggrepel)
library(forcats)
library(tidytext)

setwd("~/Dropbox/PMI_3_analyses/multi-omics_data/amplicon/16S/longitudinal/soilhip")

SHSU_vol = read.table("pc_SHSU_vol.tsv", header = TRUE, sep = "\t")

p1=ggplot(SHSU_vol, aes(x=add_0c, y= Axis1, color=sample_site)) +
  geom_rect(aes(xmin = 50, xmax = 200, ymin = -Inf, ymax = Inf),
            alpha = 1, linetype = 0,
            fill = "gray85") +
  #scale_y_continuous(breaks=pretty_breaks()) +
  #stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=0.4) +
  #stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "ADD", y = "Diversity Distance", title="STAFS Phylogenetic Volatility") +
  guides(color = guide_legend(override.aes = list(size=2))) +
  theme_bw() + 
  xlim(0,650) +
  ylim(-0.3,0.75) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=9),
        plot.title = element_text(size=12, hjust = 0.5),
        strip.text = element_text(size=6.5),
        legend.text = element_text(size=7),
        legend.title = element_blank(),
        legend.position = c(0.81, 0.86),
        legend.direction = "horizontal",
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.01, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.background = element_rect(fill="white")) +
  scale_colour_manual(values=c('Decomposition Soil'='royalblue3', 'Control Soil'='gray20')) +
  geom_vline(xintercept=50, linetype="dashed") +
  geom_vline(xintercept=200, linetype="dashed")
p1

UTK_vol = read.table("pc_UTK_vol.tsv", header = TRUE, sep = "\t")

p2=ggplot(UTK_vol, aes(x=add_0c, y= Axis1, color=sample_site)) +
  geom_rect(aes(xmin = 50, xmax = 200, ymin = -Inf, ymax = Inf),
            alpha = 1, linetype = 0,
            fill = "gray85") +
  #scale_y_continuous(breaks=pretty_breaks()) +
  #stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=0.4) +
  #stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "ADD", y = "Diversity Distance", title="ARF Phylogenetic Volatility") +
  guides(color = guide_legend(override.aes = list(size=2))) +
  theme_bw() + 
  xlim(0,525) +
  ylim(-0.3,0.75) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=9),
        plot.title = element_text(size=12, hjust = 0.5),
        strip.text = element_text(size=6.5),
        legend.text = element_text(size=7),
        legend.title = element_blank(),
        legend.position = c(0.81, 0.86),
        legend.direction = "horizontal",
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.01, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.background = element_rect(fill="white")) +
  scale_colour_manual(values=c('Decomposition Soil'='mediumseagreen', 'Control Soil'='gray20')) +
  geom_vline(xintercept=50, linetype="dashed") +
  geom_vline(xintercept=200, linetype="dashed")
p2

CMU_vol = read.table("pc_CMU_vol.tsv", header = TRUE, sep = "\t")

p3=ggplot(CMU_vol, aes(x=add_0c, y= Axis1, color=sample_site)) +
  geom_rect(aes(xmin = 50, xmax = 200, ymin = -Inf, ymax = Inf),
            alpha = 1, linetype = 0,
            fill = "gray85") +
  #scale_y_continuous(breaks=pretty_breaks()) +
  #stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  stat_summary(geom="line", fun.data=mean_se, orientation="x", size=0.4) +
  #stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "ADD", y = "Diversity Distance", title="FIRS Phylogenetic Volatility") +
  guides(color = guide_legend(override.aes = list(size=2))) +
  theme_bw() + 
  xlim(0,600) +
  ylim(-0.3,0.75) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=9),
        plot.title = element_text(size=12, hjust = 0.5),
        strip.text = element_text(size=6.5),
        legend.text = element_text(size=7),
        legend.title = element_blank(),
        legend.position = c(0.81, 0.86),
        legend.direction = "horizontal",
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.01, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.background = element_rect(fill="white")) +
  scale_colour_manual(values=c('Decomposition Soil'='mediumorchid4', 'Control Soil'='gray20')) +
  geom_vline(xintercept=50, linetype="dashed") +
  geom_vline(xintercept=200, linetype="dashed")
p3

p4 = plot_grid(p2, p3, p1, ncol = 1, nrow = 3)
p4
save_plot("vol_multiplot.png", p4, base_height = 5, base_width = 6, dpi=600)





### richness LME
SHSU_richness = read.table("lme/richness-SHSU-data.tsv", header = TRUE, sep = "\t")

p1=ggplot(SHSU_richness, aes(x=add_0c, y= observed_otus, color=sample_site)) +
  geom_rect(aes(xmin = 50, xmax = 200, ymin = -Inf, ymax = Inf),
            alpha = 1, linetype = 0,
            fill = "gray85") +
  #scale_y_continuous(breaks=pretty_breaks()) +
  #stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  #stat_summary(geom="line", fun.data=mean_se, orientation="x", size=0.4) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "ADD", y = "Observed ASVs", title="STAFS Bacterial Richness") +
  guides(color = guide_legend(override.aes = list(size=2))) +
  theme_bw() + 
  xlim(0,650) +
  ylim(0,1600) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=9),
        plot.title = element_text(size=12, hjust = 0.5),
        strip.text = element_text(size=6.5),
        legend.text = element_text(size=7),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.87),
        legend.direction = "vertical",
        legend.key.height = unit(0.001, "cm"),
        legend.key.width = unit(0.001, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.background = element_rect(fill=NA)) +
  scale_colour_manual(values=c('Decomposition Soil'='royalblue3', 'Control Soil'='gray20')) +
  geom_vline(xintercept=50, linetype="dashed") +
  geom_vline(xintercept=200, linetype="dashed")
p1

UTK_richness = read.table("lme/richness-UTK-data.tsv", header = TRUE, sep = "\t")

p2=ggplot(UTK_richness, aes(x=add_0c, y= observed_otus, color=sample_site)) +
  geom_rect(aes(xmin = 50, xmax = 200, ymin = -Inf, ymax = Inf),
            alpha = 1, linetype = 0,
            fill = "gray85") +
  #scale_y_continuous(breaks=pretty_breaks()) +
  #stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  #stat_summary(geom="line", fun.data=mean_se, orientation="x", size=0.4) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "ADD", y = "Observed ASVs", title="ARF Bacterial Richness") +
  guides(color = guide_legend(override.aes = list(size=2))) +
  theme_bw() + 
  xlim(0,650) +
  ylim(0,1600) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=9),
        plot.title = element_text(size=12, hjust = 0.5),
        strip.text = element_text(size=6.5),
        legend.text = element_text(size=7),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.87),
        legend.direction = "vertical",
        legend.key.height = unit(0.001, "cm"),
        legend.key.width = unit(0.001, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.background = element_rect(fill=NA)) +
  scale_colour_manual(values=c('Decomposition Soil'='mediumseagreen', 'Control Soil'='gray20')) +
  geom_vline(xintercept=50, linetype="dashed") +
  geom_vline(xintercept=200, linetype="dashed")
p2

CMU_richness = read.table("lme/richness-CMU-data.tsv", header = TRUE, sep = "\t")

p3=ggplot(CMU_richness, aes(x=add_0c, y= observed_otus, color=sample_site)) +
  geom_rect(aes(xmin = 50, xmax = 200, ymin = -Inf, ymax = Inf),
            alpha = 1, linetype = 0,
            fill = "gray85") +
  #scale_y_continuous(breaks=pretty_breaks()) +
  #stat_summary(geom="errorbar", fun.data=mean_se, width=0, orientation="x") +
  #stat_summary(geom="line", fun.data=mean_se, orientation="x", size=0.4) +
  stat_summary(geom="point", fun.data=mean_se, orientation="x", size=1) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "ADD", y = "Observed ASVs", title="FIRS Bacterial Richness") +
  guides(color = guide_legend(override.aes = list(size=2))) +
  theme_bw() + 
  xlim(0,650) +
  ylim(0,1600) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=9),
        plot.title = element_text(size=12, hjust = 0.5),
        strip.text = element_text(size=6.5),
        legend.text = element_text(size=7),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.87),
        legend.direction = "vertical",
        legend.key.height = unit(0.001, "cm"),
        legend.key.width = unit(0.001, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        legend.background = element_rect(fill=NA)) +
  scale_colour_manual(values=c('Decomposition Soil'='mediumorchid4', 'Control Soil'='gray20')) +
  geom_vline(xintercept=50, linetype="dashed") +
  geom_vline(xintercept=200, linetype="dashed")
p3

p4 = plot_grid(p2, p3, p1, ncol = 1, nrow = 3)
p4
save_plot("richness_multiplot.png", p4, base_height = 5, base_width = 6, dpi=600)


