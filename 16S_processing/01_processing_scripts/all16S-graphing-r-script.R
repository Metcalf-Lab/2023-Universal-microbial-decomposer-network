###############################################################################
# Project: PMI Spring 18S
# Name: Aeriel
# Generating Figures
# 02-13-2019
###############################################################################
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


setwd('~/projects/Dropbox/PMI_3_analyses/multi-omics_data/amplicon/16S/')

##
##### Richness Evaluation #####
richness <- read_tsv('03_metadata/combined-metadata-simple-oct2020_R.txt') %>% 
  filter(control == 'n') %>% 
  filter(facility != 'control') %>% 
  filter(host_subject_id != 'not_applicable') %>% 
  filter(richness_16S != '#N/A') %>% 
  mutate(richness_16S = as.numeric(richness_16S)) %>% 
  mutate(decomp_stage = fct_relevel(decomp_stage, 'early', 'active'))

richness_18S <- read_tsv('03_metadata/combined-metadata-simple-oct2020_R.txt') %>% 
  filter(control == 'n') %>% 
  filter(richness_18S != '#N/A') %>% 
  mutate(richness_18S = as.numeric(richness_18S))

richness %>% 
  group_by(season_facility) %>% 
  summarize(median = median(richness_16S))

## 16S
richness %>% 
  ggplot(aes(x = season_facility, y = richness_16S, fill = season)) +
  geom_boxplot(notch = TRUE) +
  geom_jitter(alpha = 0.3) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95))

richness %>% 
  filter(days_since_placement == '0' & sample_type == 'soil') %>% 
  ggplot(aes(x = season_facility, y = richness_16S, fill = season)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.3) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95))

richness %>% 
  ggplot(aes(x = season, y = richness_16S, fill = season)) +
  geom_boxplot(notch = TRUE) +
  geom_jitter(alpha = 0.3) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95))

richness %>% 
  ggplot(aes(x = facility, y = richness_16S, fill = facility)) +
  geom_boxplot(notch = TRUE) +
  geom_jitter(alpha = 0.3) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95))
  
richness %>% 
  filter(richness_16S <= 2500 | sample_type == 'soil')  %>% 
  ggplot(aes(x = host_subject_id, y = richness_16S, fill = season)) +
  geom_boxplot(notch = TRUE) +
  geom_jitter(alpha = 0.3) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95)) +
  facet_wrap(. ~ sample_type, scales = 'free', ncol = 1)

richness %>% 
  filter(richness_16S <= 3000 | sample_type == 'soil')  %>% 
  filter(extraction_plate != '#N/A') %>% 
  mutate(extraction_plate = as.character(extraction_plate)) %>% 
  ggplot(aes(x = extraction_plate, y = richness_16S, fill = sample_type)) +
  geom_boxplot(notch = TRUE) +
  geom_jitter(alpha = 0.3) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95)) +
  facet_wrap(. ~ sample_type, scales = 'free', ncol = 1)

richness %>% 
  filter(richness_16S <= 1000 | sample_type == 'soil')  %>% 
  ggplot(aes(x = season_facility, y = richness_16S, fill = season)) +
  geom_boxplot(notch = TRUE) +
  geom_jitter(alpha = 0.3) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95)) +
  facet_wrap(. ~ sample_type, scales = 'free', ncol = 1)

richness %>% 
  mutate(add_0c = as.numeric(add_0c)) %>% 
  ggplot(aes(x = add_0c, y = richness_16S, color = decomp_stage)) +
  geom_point() +
  theme_classic() +
  facet_grid(facility ~ sample_type, scales = 'free') +
  geom_smooth(method = "loess", se = FALSE, color = 'black') +
  scale_color_viridis(discrete =  TRUE, end = 0.8) +
  labs(x = 'Accumulated Degree Days',
       y = '16S Richness',
       color = 'Decomposition \nStage')

####FIGURES FOR PUB ###
richness <- read_tsv('03_metadata/combined-metadata-simple-oct2020_R.txt') %>% 
  filter(control == 'n') %>% 
  filter(facility != 'control') %>% 
  filter(host_subject_id != 'not_applicable') %>% 
  filter(richness_16S != '#N/A') %>% 
  filter(faiths_16S != '#N/A') %>% 
  mutate(richness_16S = as.numeric(richness_16S)) %>% 
  mutate(decomp_stage = fct_relevel(decomp_stage, 'early', 'active'))


richness %>% 
  mutate(add_0c = as.numeric(add_0c)) %>% 
  ggplot(aes(x = add_0c, y = faiths_16S, color = decomp_stage)) +
  geom_point() +
  theme_classic() +
  facet_grid(facility ~ sample_type, scales = 'free') +
  geom_smooth(method = "loess", se = FALSE, color = 'black') +
  scale_color_viridis(discrete =  TRUE, end = 0.8) +
  labs(x = 'Accumulated Degree Days',
       y = "16S Faith's Phylogenetic Diversity",
       color = 'Decomposition \nStage') -> faiths_16S_plot

ggsave('16S-faiths-add.png', plot = faiths_16S_plot, path = '08_results/figures/',
      device = 'png', width = 8, height = 4, units = "in")



richness %>% 
  mutate(add_0c = as.numeric(add_0c)) %>% 
  ggplot(aes(x = add_0c, y = faiths_18S, color = decomp_stage)) +
  geom_point() +
  theme_classic() +
  facet_grid(facility ~ sample_type, scales = 'free') +
  geom_smooth(method = "loess", se = FALSE, color = 'black') +
  scale_color_viridis(discrete =  TRUE, end = 0.8) +
  labs(x = 'Accumulated Degree Days',
       y = "18S Faith's Phylogenetic Diversity",
       color = 'Decomposition \nStage') -> faiths_18S_plot

ggsave('18S-faiths-add.png', plot = faiths_18S_plot, path = '../18S/08_results/figures/',
       device = 'png', width = 8, height = 4, units = "in")


####Day 0 ###
alpha <- read_tsv('03_metadata/combined-metadata-simple-oct2020_R.txt') %>% 
  filter(control == 'n') %>% 
  filter(facility != 'control') %>% 
  filter(host_subject_id != 'not_applicable') %>% 
  filter(richness_16S != '#N/A') %>% 
  filter(faiths_16S != '#N/A') %>% 
  mutate(richness_16S = as.numeric(richness_16S))

alpha %>% 
  filter(days_since_placement == '0') %>% 
  
  gather(alpha_metric, alpha_diversity, 
         "richness_16S", "richness_18S", "faiths_16S", "faiths_18S") %>%
  mutate(alpha_metric = fct_relevel(alpha_metric, "richness_16S", "richness_18S", "faiths_16S", "faiths_18S"),
         alpha_metric = fct_recode(alpha_metric, "16S Richness" = "richness_16S"),
         alpha_metric = fct_recode(alpha_metric, "18S Richness" = "richness_18S"),
         alpha_metric = fct_recode(alpha_metric, "16S Faiths" = "faiths_16S"),
         alpha_metric = fct_recode(alpha_metric, "18S Faiths" = "faiths_18S")) %>% 
  ggplot(aes(x = facility, y = alpha_diversity, color = facility)) +
  geom_boxplot() +
  geom_jitter() +
  theme_classic() +
  theme(legend.position = 'none') +
  facet_grid(alpha_metric ~ ., scales = 'free') +
  scale_color_viridis(discrete =  TRUE, end = 0.8) +
  labs(x = 'Facility',
       y = "Alpha Diversity") -> day0alpha

ggsave('Day0Alpha', day0alpha, 'png', '08_results/figures/', width = 5, height = 8)




#### ### ###


read_tsv('~/Downloads/faiths.txt') %>% 
  group_by(decomp_stage) %>% 
  filter(sample_type == 'skin') %>% 
  summarize(faith = mean(faith_pd),
            rich = mean(richness_16S))

read_tsv('~/Downloads/faiths.txt') %>% 
  group_by(facility) %>% 
  filter(sample_type == 'soil' & days_since_placement == '0') %>% 
  summarize(faith = mean(faith_pd),
            rich = mean(richness_16S)) 
  
  
pairwise.wilcox.test(stat$faith_pd, stat$facility)

pairwise.wilcox.test(stat$richness_16S, stat$facility)



read_tsv('~/Downloads/alpha_18.txt') %>% 
  group_by(decomp_stage) %>% 
  filter(sample_type == 'skin') %>% 
  summarize(faith = mean(faith_pd),
            rich = mean(richness_18S))

read_tsv('~/Downloads/alpha_18.txt') %>% 
  filter(sample_type == 'soil' & days_since_placement == '0') %>% 
  group_by(facility) %>% 
  summarize(faith = mean(faith_pd),
            rich = mean(richness_18S))

pairwise.wilcox.test(stat2$faith_pd, stat2$facility)
  
pairwise.wilcox.test(stat2$richness_18S, stat2$facility)


richness_18S %>% 
  ggplot(aes(x = season_facility, y = richness_18S, fill = season)) +
  geom_boxplot(notch = TRUE) +
  geom_jitter(alpha = 0.3) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95))

richness_18S %>% 
  filter(days_since_placement == '0' & sample_type == 'soil') %>% 
  ggplot(aes(x = season_facility, y = richness_18S, fill = season)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.3) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95))

richness_18S %>% 
  ggplot(aes(x = season_facility, y = richness_18S, fill = season)) +
  geom_boxplot(notch = TRUE) +
  geom_jitter(alpha = 0.3) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95)) +
  facet_wrap(. ~ sample_type, scales = 'free', ncol = 1)

richness_18S %>% 
  ggplot(aes(x = season, y = richness_18S, fill = season)) +
  geom_boxplot(notch = TRUE) +
  geom_jitter(alpha = 0.3) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95))

richness_18S %>% 
  ggplot(aes(x = facility, y = richness_18S, fill = facility)) +
  geom_boxplot(notch = TRUE) +
  geom_jitter(alpha = 0.3) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95))

richness_18S %>% 
  ggplot(aes(x = host_subject_id, y = richness_18S, fill = season)) +
  geom_boxplot(notch = TRUE) +
  geom_jitter(alpha = 0.3) +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95))

richness_18S %>% 
  filter(extraction_plate != '#N/A') %>% 
  mutate(extraction_plate = as.character(extraction_plate)) %>% 
  ggplot(aes(x = extraction_plate, y = richness_18S)) +
  geom_boxplot(notch = TRUE) +
  geom_jitter(alpha = 0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95))


richness %>% 
  mutate(add_0c = as.numeric(add_0c)) %>% 
  ggplot(aes(x = add_0c, y = faiths_18S, color = decomp_stage)) +
  geom_point() +
  theme_classic() +
  facet_grid(facility ~ sample_type, scales = 'free') +
  geom_smooth(method = "loess", se = FALSE, color = 'black') +
  scale_color_viridis(discrete =  TRUE, end = 0.8) +
  labs(x = 'Accumulated Degree Days',
       y = "18S Faith's Phylogenetic Diversity",
       color = 'Decomposition \nStage') -> faiths_18S_plot

ggsave('18S-faiths-add.png', plot = faiths_18S_plot, path = '../18S/08_results/figures/',
       device = 'png', width = 8, height = 4, units = "in")




##
##### Taxonomy barplots#####
## Phylum
phylum_tax <- read_tsv('08_results/rel-freq-tables/phylum_relfreq_table.tsv')
metadata <- read_tsv('03_metadata/combined-metadata-simple-oct2020_R.txt')

phylum_tax %>% 
  gather(key = 'taxon', value = 'abundance', -sample) %>% 
  group_by(taxon) %>% 
  mutate_at(.vars = vars(abundance), .funs = list(sum = sum)) %>% 
  ungroup() %>% 
  mutate(taxon = replace(taxon, sum <= 10, 'other'),
         taxon = gsub('D_0__', '', taxon),
         taxon = gsub('D_1__', ' ', taxon)) %>% 
  left_join(metadata, by = 'sample') %>% 
  mutate(decomp_stage = fct_relevel(decomp_stage, 'early', 'active'),) %>% 
  ggplot(aes(x = reorder(sample, add_0c), y = abundance, fill = taxon)) +
  geom_bar(stat = 'identity', width = 1) +
  facet_wrap(. ~ sample_type + facility, scales = 'free') +
  theme_few() +
  theme(legend.position = 'bottom',
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(fill = 'Bacterial Phyla', y = 'Relative Abundance')

## Genus
genus_tax <- read_tsv('08_results/rel-freq-tables/genus_relfreq_table.tsv')
metadata <- read_tsv('03_metadata/combined-metadata-simple-oct2020_R.txt')

#get a count of genera per phylum
genus_tax %>% 
  gather(key = 'taxon', value = 'abundance', -sample) %>% 
  #group_by(taxon) %>% 
  #mutate_at(.vars = vars(abundance), .funs = list(sum = sum)) %>% 
  #ungroup() %>% 
  group_by(sample) %>% 
  mutate_at(.vars = vars(abundance), .funs = list(total_abundance = sum)) %>% 
  ungroup() %>% 
  mutate(percent_abundance = (abundance/total_abundance)*100) %>%  
  left_join(metadata, by = 'sample') %>% 
  mutate(decomp_stage = fct_relevel(decomp_stage, 'early', 'active')) %>% 
  filter(sample_type == 'soil') %>% 
  mutate(taxon = gsub('D_0__Bacteria;', '', taxon),
         taxon = gsub('D_1__', ' ', taxon),
         taxon = gsub('D_2__', ' ', taxon),
         taxon = gsub('D_3__', ' ', taxon),
         taxon = gsub('D_4__', ' ', taxon),
         taxon = gsub('D_5__', ' ', taxon)) %>% 
  mutate(taxon = replace(taxon, percent_abundance <= 10, 'other'),
         taxon = replace(taxon, n_distinct(taxon) <= 2, 'other')) %>%  
  separate(taxon, 
           into = c('kingdom', 'phylum', 'class', 'order', 'family', 'genus'), 
           sep = ';', remove = FALSE) %>%
  mutate(taxon = replace(taxon, genus == '__', 'other'),
         taxon = replace(taxon, family == '__', 'other'),
         taxon = replace(taxon, order == '__', 'other'),
         taxon = replace(taxon, class == '__', 'other')) %>% 
  distinct(taxon) -> distinct_taxa




#make palettes
actino_pal <- colorRampPalette(rev(brewer.pal(3, "YlOrRd")))
bactero_pal <- colorRampPalette(rev(brewer.pal(3, "YlGnBu")))
cyano_pal <- colorRampPalette(rev(brewer.pal(1, "YlOrBr")))
firmi_pal <- colorRampPalette(rev(brewer.pal(6, "Purples")))
proteo_pal <- colorRampPalette(rev(brewer.pal(11, "YlGn")))
verruco_pal <- colorRampPalette(rev(brewer.pal(1, "OrRd")))
other_pal <- colorRampPalette(rev(brewer.pal(1, "Greys")))


#final plot - soil
genus_tax %>% 
  gather(key = 'taxon', value = 'abundance', -sample) %>% 
  #group_by(taxon) %>% 
  #mutate_at(.vars = vars(abundance), .funs = list(sum = sum)) %>% 
  #ungroup() %>% 
  group_by(sample) %>% 
  mutate_at(.vars = vars(abundance), .funs = list(total_abundance = sum)) %>% 
  ungroup() %>% 
  mutate(percent_abundance = (abundance/total_abundance)*100) %>%  
  left_join(metadata, by = 'sample') %>% 
  mutate(decomp_stage = fct_relevel(decomp_stage, 'early', 'active')) %>% 
  filter(sample_type == 'soil') %>% 
  mutate(taxon = gsub('D_0__Bacteria;', '', taxon),
         taxon = gsub('D_1__', ' ', taxon),
         taxon = gsub('D_2__', ' ', taxon),
         taxon = gsub('D_3__', ' ', taxon),
         taxon = gsub('D_4__', ' ', taxon),
         taxon = gsub('D_5__', ' ', taxon)) %>% 
  mutate(taxon = replace(taxon, percent_abundance <= 1, 'other'),
         taxon = replace(taxon, n_distinct(taxon) <= 2, 'other')) %>%  
  separate(taxon, 
           into = c('kingdom', 'phylum', 'class', 'order', 'family', 'genus'), 
           sep = ';', remove = FALSE) %>%
  mutate(taxon = replace(taxon, genus == '__', 'other'),
         taxon = replace(taxon, family == '__', 'other'),
         taxon = replace(taxon, order == '__', 'other'),
         taxon = replace(taxon, class == '__', 'other')) %>% 
  ggplot(aes(x = reorder(sample, add_0c), y = abundance, fill = taxon)) +
  geom_bar(stat = 'identity', width = 1, 
           position = position_fill(reverse = TRUE)) +
  facet_wrap(. ~ facility, scales = 'free') +
  theme_few() +
  #scale_fill_manual(values = c(actino_pal(3), bactero_pal(3), cyano_pal(1),
  #                             firmi_pal(6), proteo_pal(11), verruco_pal(1),
  #                             other_pal(1))) +
  theme(legend.position = 'right',
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(fill = 'Bacterial Genera', y = 'Relative Abundance') +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(labels = scales::percent)

#final plot - skin
genus_tax %>% 
  gather(key = 'taxon', value = 'abundance', -sample) %>% 
  group_by(taxon) %>% 
  mutate_at(.vars = vars(abundance), .funs = list(sum = sum)) %>% 
  ungroup() %>% 
  left_join(metadata, by = 'sample') %>% 
  mutate(decomp_stage = fct_relevel(decomp_stage, 'early', 'active')) %>% 
  filter(sample_type == 'skin') %>% 
  mutate(taxon = gsub('D_0__Bacteria;', '', taxon),
         taxon = gsub('D_1__', ' ', taxon),
         taxon = gsub('D_2__', ' ', taxon),
         taxon = gsub('D_3__', ' ', taxon),
         taxon = gsub('D_4__', ' ', taxon),
         taxon = gsub('D_5__', ' ', taxon)) %>% 
  mutate(taxon = replace(taxon, sum <= 20, 'other'),
         taxon = replace(taxon, n_distinct(taxon) <= 2, 'other')) %>%  
  separate(taxon, 
           into = c('kingdom', 'phylum', 'class', 'order', 'family', 'genus'), 
           sep = ';', remove = FALSE) %>%
  mutate(taxon = replace(taxon, genus == '__', 'other'),
         taxon = replace(taxon, family == '__', 'other'),
         taxon = replace(taxon, order == '__', 'other'),
         taxon = replace(taxon, class == '__', 'other')) %>% 
  ggplot(aes(x = reorder(sample, add_0c), y = abundance, fill = taxon)) +
  geom_bar(stat = 'identity', width = 1, 
           position = position_fill(reverse = TRUE)) +
  facet_wrap(. ~ facility, scales = 'free') +
  theme_few() +
 # scale_fill_manual(values = c(actino_pal(3), bactero_pal(3), cyano_pal(1),
  #                             firmi_pal(6), proteo_pal(11), verruco_pal(1),
   #                            other_pal(1))) +
  theme(legend.position = 'right',
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(fill = 'Bacterial Genera', y = 'Relative Abundance') +
  guides(fill=guide_legend(ncol=1))

##### Taxonomy heatmaps #####
phyloseq <- qza_to_phyloseq(
  features="04_artifacts-and-visualizations/PMI-16S-nochlomito-filtered-table.qza",
  tree="05_trees/16S_sepp_tree_filtered.qza",
  taxonomy="07_taxonomy/pmi3-16S-nochlomito-silva-taxonomy.qza",
  metadata="03_metadata/combined-metadata-simple-nov2020.txt"
)


#### Try just ggplot
## 16S
genus_tax <- read_tsv('08_results/rel-freq-tables/genus_relfreq_table.tsv')
metadata <- read_tsv('03_metadata/combined-metadata-simple-oct2020_R.txt')

genus_tax %>% 
  gather(key = 'taxon', value = 'abundance', -sample) %>% 
  left_join(metadata, by = 'sample') %>% 
  mutate(taxon = gsub('D_0__Bacteria;', '', taxon),
         taxon = gsub('D_1__', ' ', taxon),
         taxon = gsub('D_2__', ' ', taxon),
         taxon = gsub('D_3__', ' ', taxon),
         taxon = gsub('D_4__', ' ', taxon),
         taxon = gsub('D_5__', ' ', taxon),
         taxon = replace(taxon, abundance <= 0.1, 'other;other;other;other;other;other;other'),
         decomp_stage = fct_relevel(decomp_stage, 'early', 'active')) %>%
  filter(sample_type != 'control',
         facility != 'control') %>% 
  group_by(decomp_stage, facility, taxon, sample_type) %>% 
  summarise(rel_abundance = mean(abundance, na.rm = TRUE)) -> partial

## 18S
genus_tax_18S <- read_tsv('../18S/08_results/rel-freq-tables/genus_relfreq_table.tsv')

genus_tax_18S %>% 
  gather(key = 'taxon', value = 'abundance', -sample) %>% 
  left_join(metadata, by = 'sample') %>% 
  mutate(taxon = gsub('d__Eukaryota;', '', taxon),
         taxon = gsub('k__', ' ', taxon),
         taxon = gsub('p__', ' ', taxon),
         taxon = gsub('ps__', ' ', taxon),
         taxon = gsub('c__', ' ', taxon),
         taxon = gsub('cs__', ' ', taxon),
         taxon = gsub('o__', ' ', taxon),
         taxon = gsub('os__', ' ', taxon),
         taxon = gsub('f__', ' ', taxon),
         taxon = gsub('fs__', ' ', taxon),
         taxon = replace(taxon, abundance <= 0.2, 'other;other;other;other;other;other;other'),
         decomp_stage = fct_relevel(decomp_stage, 'early', 'active')) %>%
  filter(sample_type != 'control',
         facility != 'control') %>% 
  group_by(decomp_stage, facility, taxon, sample_type) %>% 
  summarise(rel_abundance = mean(abundance, na.rm = TRUE)) %>% 
  mutate(euk = 'euk') %>% 
  unite(col = taxon, euk, taxon, sep = ';') %>% 
  complete(facility, decomp_stage, sample_type, taxon, fill = list(rel_abundance = 0)) %>% 
  separate(taxon, 
           into = c('domain', 'phylum', 'phylum2', 'subphylum', 'class', 
                    'subclass','order', 'suborder', 'family', 'subfamily'), 
           sep = ';', remove = FALSE) -> partial_18


## Combining 16S and 18S
partial %>% 
  mutate(bac = 'bac') %>% 
  unite(col = taxon, bac, taxon, sep = ';') %>% 
  complete(facility, decomp_stage, sample_type, taxon, fill = list(rel_abundance = 0)) %>% 
  separate(taxon, 
           into = c('domain', 'phylum', 'class', 'order', 'family', 'genus'), 
           sep = ';', remove = FALSE) %>% 
  bind_rows(partial_18) %>% 
  filter(facility != 'control',
         sample_type != 'NA',
         decomp_stage != 'not_applicable',
         phylum != ' Chloroflexi',
         phylum != ' Deinococcus-Thermus',
         phylum != ' Cyanobacteria',
         phylum != ' Chloroplastida',
         order != '__',
         order != ' uncultured',
         order != ' Unknown Family',
         order != ' uncultured Acidobacteriales bacterium',
         order != ' 11-24',
         order != ' Longimicrobiales',
         order != ' Acidobacteria bacterium SCN 69-37',
         order != 'other',
         order != ' LKM11',
         order != ' LKM15',
         order != ' SG8-5',
         order != ' Mammalia',
         order != ' Incertae_Sedis',
         order != ' Cryptomonadales',
         order !=  ' metagenome') %>% 
  unite(col = phylum_order, phylum, order, remove = FALSE) %>% 
  group_by(phylum_order) %>% 
  mutate(phylum_order = cur_group_id()) %>% 
  ungroup() %>% 
  mutate(phylum_order = as.numeric(phylum_order), 
         order = fct_reorder(order, phylum_order, .desc = TRUE)) %>% 
  ggplot(aes(x = decomp_stage, y = order, fill = rel_abundance)) +
  geom_tile() +
  facet_grid(domain ~ sample_type + facility, scales = 'free_y') +
  theme_classic() +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.spacing.x = unit(0, 'lines'),
        strip.text.y = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  labs(fill = 'Relative \nAbundance') +
  scale_fill_viridis(begin = 0.15, end = 0.8) -> combined_taxa_plot

ggsave('all_taxa.png', plot = combined_taxa_plot, path = '../figures/',
       device = 'png', width = 10, height = 20, units = "in")

##ORDER
order_tax <- read_tsv('08_results/rel-freq-tables/order_relfreq_table.tsv')
metadata <- read_tsv('03_metadata/combined-metadata-simple-oct2020_R.txt')

order_tax %>% 
  gather(key = 'taxon', value = 'abundance', -sample) %>% 
  left_join(metadata, by = 'sample') %>% 
  mutate(taxon = gsub('D_0__Bacteria;', '', taxon),
         taxon = gsub('D_1__', ' ', taxon),
         taxon = gsub('D_2__', ' ', taxon),
         taxon = gsub('D_3__', ' ', taxon),
         taxon = replace(taxon, abundance <= 0.15, 'other;other;other;other'),
         decomp_stage = fct_relevel(decomp_stage, 'early', 'active')) %>%
  filter(sample_type != 'control',
         facility != 'control') %>% 
  group_by(decomp_stage, facility, taxon, sample_type) %>% 
  summarise(rel_abundance = mean(abundance, na.rm = TRUE)) -> partial

## 18S
genus_tax_18S <- read_tsv('../18S/08_results/rel-freq-tables/genus_relfreq_table.tsv')

genus_tax_18S %>% 
  gather(key = 'taxon', value = 'abundance', -sample) %>% 
  left_join(metadata, by = 'sample') %>% 
  mutate(taxon = gsub('d__Eukaryota;', '', taxon),
         taxon = gsub('k__', ' ', taxon),
         taxon = gsub('p__', ' ', taxon),
         taxon = gsub('ps__', ' ', taxon),
         taxon = gsub('c__', ' ', taxon),
         taxon = gsub('cs__', ' ', taxon),
         taxon = gsub('o__', ' ', taxon),
         taxon = gsub('os__', ' ', taxon),
         taxon = gsub('f__', ' ', taxon),
         taxon = gsub('fs__', ' ', taxon),
         taxon = replace(taxon, abundance <= 0.3, 'other;other;other;other;other;other;other'),
         decomp_stage = fct_relevel(decomp_stage, 'early', 'active')) %>%
  filter(sample_type != 'control',
         facility != 'control') %>% 
  group_by(decomp_stage, facility, taxon, sample_type) %>% 
  summarise(rel_abundance = mean(abundance, na.rm = TRUE)) %>% 
  mutate(euk = 'euk') %>% 
  unite(col = taxon, euk, taxon, sep = ';') %>% 
  complete(facility, decomp_stage, sample_type, taxon, fill = list(rel_abundance = 0)) %>% 
  separate(taxon, 
           into = c('domain', 'phylum', 'phylum2', 'subphylum', 'class', 
                    'subclass','order', 'suborder', 'family', 'subfamily'), 
           sep = ';', remove = FALSE) -> partial_18


## Combining 16S and 18S
partial %>% 
  mutate(bac = 'bac') %>% 
  unite(col = taxon, bac, taxon, sep = ';') %>% 
  complete(facility, decomp_stage, sample_type, taxon, fill = list(rel_abundance = 0)) %>% 
  separate(taxon, 
           into = c('domain', 'phylum', 'class', 'order'), 
           sep = ';', remove = FALSE) %>% 
  bind_rows(partial_18) %>% 
  filter(facility != 'control',
         sample_type != 'NA',
         decomp_stage != 'not_applicable',
         phylum != ' Chloroflexi',
         phylum != ' Deinococcus-Thermus',
         phylum != ' Cyanobacteria',
         phylum != ' Chloroplastida',
         order != '__',
         order != ' uncultured',
         order != ' Unknown Family',
         order != ' uncultured Acidobacteriales bacterium',
         order != ' 11-24',
         order != ' Longimicrobiales',
         order != ' Acidobacteria bacterium SCN 69-37',
         order != 'other',
         order != ' LKM11',
         order != ' LKM15',
         order != ' SG8-5',
         order != ' Mammalia',
         order != ' Incertae_Sedis',
         order != ' Cryptomonadales',
         order !=  ' metagenome') %>% 
  unite(col = phylum_order, phylum, order, remove = FALSE) %>% 
  group_by(phylum_order) %>% 
  mutate(phylum_order = cur_group_id()) %>% 
  ungroup() %>% 
  mutate(phylum_order = as.numeric(phylum_order), 
         order = fct_reorder(order, phylum_order, .desc = TRUE)) %>% 
  ggplot(aes(x = decomp_stage, y = order, fill = rel_abundance)) +
  geom_tile() +
  facet_grid(domain ~ sample_type + facility, scales = 'free_y') +
  theme_classic() +
  theme(legend.position = 'right',
        text = element_text(size = 18, family = 'Times New Roman'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.spacing.x = unit(0, 'lines'),
        strip.text.y = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        panel.background = element_rect(fill = "#453781")) +
  labs(fill = 'Relative \nAbundance') +
  scale_fill_viridis(begin = 0.15, end = 0.8) -> combined_taxa_plot

ggsave('all_taxa.png', plot = combined_taxa_plot, path = '../figures/',
       device = 'png', width = 10, height = 20, units = "in")

## Wanna see if there are patterns by body
#genus

genus_tax %>% 
  summarize_if(is.numeric, mean, na.rm = TRUE) -> genus_totals

genus_tax %>% 
  bind_rows(totals) %>%
  gather(key = 'taxon', value = 'abundance', -sample) %>% 
  left_join(metadata, by = 'sample') %>% 
  mutate(taxon = gsub('D_0__Bacteria;', '', taxon),
         taxon = gsub('D_1__', ' ', taxon),
         taxon = gsub('D_2__', ' ', taxon),
         taxon = gsub('D_3__', ' ', taxon),
         taxon = gsub('D_4__', ' ', taxon),
         taxon = gsub('D_5__', ' ', taxon),
         taxon = fct_reorder(taxon, sample == 'NA')) %>%
  filter(sample_type != 'control',
         facility != 'control',
         abundance <= 0.15) %>% 
  complete(host_subject_id, days_since_placement, 
           taxon, fill = list(abundance = 0)) -> bacteria_genus

bacteria_genus %>% 
  mutate(days_since_placement = as.numeric(days_since_placement)) %>% 
  filter(taxon != 'other', 
         taxon != 'D_0__Archaea; Thaumarchaeota; Nitrososphaeria') %>% 
  ggplot(aes(x = days_since_placement, y = taxon, fill = abundance)) +
  geom_tile() +
  facet_wrap(. ~ host_subject_id, nrow = 3) +
  theme_classic() +
  theme(legend.position = 'right',
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.spacing.x = unit(0, 'lines'),
        strip.text.y = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  labs(fill = 'Relative \nAbundance') +
  scale_fill_viridis(begin = 0.15, end = 0.8) 


# Class
class_tax <- read_tsv('08_results/rel-freq-tables/class_relfreq_table.tsv')

class_tax %>% 
  summarize_if(is.numeric, mean, na.rm = TRUE) -> totals

class_tax %>% 
  bind_rows(totals) %>% 
  gather(key = 'taxon', value = 'abundance', -sample) %>% 
  left_join(metadata, by = 'sample') %>% 
  mutate(taxon = gsub('D_0__Bacteria;', '', taxon),
         taxon = gsub('D_1__', ' ', taxon),
         taxon = gsub('D_2__', ' ', taxon),
         taxon = gsub('D_3__', ' ', taxon),
         taxon = gsub('D_4__', ' ', taxon),
         taxon = gsub('D_5__', ' ', taxon),
         taxon = replace(taxon, abundance <= 0.1, 'other'),
         taxon = fct_reorder(taxon, sample == 'NA')) %>%
  filter(sample_type != 'control',
         facility != 'control') %>% 
  complete(host_subject_id, days_since_placement, 
           taxon, fill = list(abundance = 0))  -> bacteria

bacteria %>% 
  mutate(days_since_placement = as.numeric(days_since_placement)) %>% 
  filter(taxon != 'other', 
         taxon != 'D_0__Archaea; Thaumarchaeota; Nitrososphaeria') %>% 
  ggplot(aes(x = days_since_placement, y = taxon, fill = abundance)) +
  geom_tile() +
  facet_wrap(. ~ host_subject_id, nrow = 3) +
  theme_classic() +
  theme(legend.position = 'right',
        axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.spacing.x = unit(0, 'lines'),
        strip.text.y = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  labs(fill = 'Relative \nAbundance') +
  scale_fill_viridis(begin = 0.15, end = 0.8) 
  

##











  

##
##### Longitudinal #####
volatility <- read_tsv('06_core-metrics/filtered-core-metrics-5000/longitudinal/beta_lonitudinal.tsv')

volatility %>% 
  filter(id != '#q2:types') %>% 
  mutate(add_0c = as.numeric(add_0c),
         `Axis 1` = as.numeric(`Axis 1`)) %>% 
  ggplot(aes(x = add_0c, y = `Axis 1`, color = sample_site)) +
  geom_line(aes(group = indiv_id), alpha = 0.3) +
  geom_smooth(se = FALSE) +
  facet_grid(facility ~ .) +
  theme_clean() +
  labs(x = 'Sampling Event',
       y = 'Distance From Previous Event',
       color = 'Room Function')

##
##### Decay Test #####
## Load Data
cmu_control <- read.csv('ecology_tables/outputs/CMU_control_soil_S3_Community_bNTI_TotalCounts.csv',
                        row.names = 1)
cmu_early <- read.csv('ecology_tables/outputs/CMU_early_soil_S3_Community_bNTI_TotalCounts.csv',
                      row.names = 1)
cmu_active <- read.csv("ecology_tables/outputs/CMU_active_soil_S3_Community_bNTI_TotalCounts.csv",
                      row.names = 1)
cmu_advanced <- read.csv('ecology_tables/outputs/CMU_advanced_soil_S3_Community_bNTI_TotalCounts.csv',
                         row.names = 1)

shsu_control <- read.csv('ecology_tables/outputs/SHSU_control_soil_S3_Community_bNTI_TotalCounts.csv',
                         row.names = 1)
shsu_early <- read.csv('ecology_tables/outputs/SHSU_early_soil_S3_Community_bNTI_TotalCounts.csv',
                       row.names = 1)
shsu_active <- read.csv('ecology_tables/outputs/SHSU_active_soil_S3_Community_bNTI_TotalCounts.csv',
                        row.names = 1)
shsu_advanced <- read.csv('ecology_tables/outputs/SHSU_advanced_soil_S3_Community_bNTI_TotalCounts.csv',
                        row.names = 1)

utk_control <- read.csv('ecology_tables/outputs/',
                        row.names = 1)
utk_early <- read.csv('ecology_tables/outputs/UTK_early_soil_S3_Community_bNTI_TotalCounts.csv',
                      row.names = 1)
utk_active <- read.csv('ecology_tables/outputs/UTK_active_soil_S3_Community_bNTI_TotalCounts.csv',
                        row.names = 1)
utk_advanced <- read.csv('ecology_tables/outputs/UTK_advanced_soil_S3_Community_bNTI_TotalCounts.csv',
                        row.names = 1)
    

# Turning bNTI into two-sided matrix
cmu_control = as.matrix(as.dist(cmu_control))
diag(cmu_control) = NA
cmu_early = as.matrix(as.dist(cmu_early))
diag(cmu_early) = NA
cmu_active = as.matrix(as.dist(cmu_active))
diag(cmu_active) = NA
cmu_advanced = as.matrix(as.dist(cmu_advanced))
diag(cmu_advanced) = NA

shsu_control = as.matrix(as.dist(shsu_control))
diag(shsu_control) = NA
shsu_early = as.matrix(as.dist(shsu_early))
diag(shsu_early) = NA
shsu_active = as.matrix(as.dist(shsu_active))
diag(shsu_active) = NA
shsu_advanced = as.matrix(as.dist(shsu_advanced))
diag(shsu_advanced) = NA

utk_control = as.matrix(as.dist(utk_control))
diag(utk_control) = NA
utk_early = as.matrix(as.dist(utk_early))
diag(utk_early) = NA
utk_active = as.matrix(as.dist(utk_active))
diag(utk_active) = NA
utk_advanced = as.matrix(as.dist(utk_advanced))
diag(utk_advanced) = NA


## Make metadata
metadata = read_tsv('03_metadata/combined-metadata-simple-oct2020_R.txt') %>% 
  dplyr::rename('Samples' = sample)

cmu_control_factors = data.frame(Samples = colnames(cmu_control))
cmu_early_factors = data.frame(Samples = colnames(cmu_early))
cmu_active_factors = data.frame(Samples = colnames(cmu_active))
cmu_advanced_factors = data.frame(Samples = colnames(cmu_advanced))

shsu_control_factors = data.frame(Samples = colnames(shsu_control))
shsu_early_factors = data.frame(Samples = colnames(shsu_early))
shsu_active_factors = data.frame(Samples = colnames(shsu_active))
shsu_advanced_factors = data.frame(Samples = colnames(shsu_advanced))
  
utk_control_factors = data.frame(Samples = colnames(utk_control))           
utk_early_factors = data.frame(Samples = colnames(utk_early))
utk_active_factors = data.frame(Samples = colnames(utk_active))
utk_advanced_factors = data.frame(Samples = colnames(utk_advanced))


# Adding in factors
cmu_control= data.frame(cmu_control, Samples = cmu_control_factors$Samples) %>% 
  melt(id.vars = c("Samples"))
cmu_early = data.frame(cmu_early, Samples = cmu_early_factors$Samples) %>% 
  melt(id.vars = c("Samples"))
cmu_active = data.frame(cmu_active, Samples = cmu_active_factors$Samples) %>% 
  melt(id.vars = c("Samples"))
cmu_advanced = data.frame(cmu_advanced, Samples = cmu_advanced_factors$Samples) %>% 
  melt(id.vars = c("Samples"))

shsu_control = data.frame(shsu_control, Samples = shsu_control_factors$Samples) %>% 
  melt(id.vars = c("Samples"))
shsu_early = data.frame(shsu_early, Samples = shsu_early_factors$Samples) %>% 
  melt(id.vars = c("Samples"))
shsu_active = data.frame(shsu_active, Samples = shsu_active_factors$Samples) %>% 
  melt(id.vars = c("Samples"))
shsu_advanced = data.frame(shsu_advanced, Samples = shsu_advanced_factors$Samples) %>% 
  melt(id.vars = c("Samples"))

utk_control = data.frame(utk_control, Samples = utk_control_factors$Samples) %>% 
  melt(id.vars = c("Samples"))
utk_early = data.frame(utk_early, Samples = utk_early_factors$Samples) %>% 
  melt(id.vars = c("Samples"))
utk_active = data.frame(utk_active, Samples = utk_active_factors$Samples) %>% 
  melt(id.vars = c("Samples"))
utk_advanced = data.frame(utk_advanced, Samples = utk_advanced_factors$Samples) %>% 
  melt(id.vars = c("Samples"))

# All data boxplots
combined <- cmu_control %>% 
  bind_rows(cmu_early, cmu_active, cmu_advanced) %>% 
  bind_rows(shsu_control, shsu_early, shsu_active, shsu_advanced) %>% 
  bind_rows(utk_early, utk_active, utk_advanced) %>% 
  left_join(metadata, by = 'Samples') %>% 
  select('Samples', 'days_since_placement', 'value', 'facility', 'decomp_stage',
         'host_number', 'season') %>% 
  mutate(decomp_stage = fct_relevel(decomp_stage, 'control', 'early', 'active'),
         color = ifelse(value < -2, '#380F63', ifelse(value > 2, '#FCAE12', '#000000'))) %>% 
  filter(facility != 'NA', decomp_stage != 'NA')

# Plot of all BNTI
combined %>% 
  ggplot(aes(x = decomp_stage, y = value)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -2, 
           fill = "#380F63", alpha = .3, color = NA) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = Inf, 
           fill = "#FCAE12", alpha = .3, color = NA) +
  geom_hline(yintercept = 2) +
  geom_hline(yintercept = -2) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  facet_grid(. ~ facility) +
  theme_classic() +
  theme(panel.border = element_rect(color="black", fill = NA),
        text = element_text(family = 'Times New Roman', size = 18)) +
  labs(x = 'Decomposition Stage',
       y = 'bNTI Value')

# Plot of absolute value BNTI
combined %>% 
  ggplot(aes(x = decomp_stage, y = abs(value))) +
  #annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -2, 
  #         fill = "#380F63", alpha = .3, color = NA) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = Inf, 
           fill = "#FCAE12", alpha = .3, color = NA) +
  geom_hline(yintercept = 2) +
  #geom_hline(yintercept = -2) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  facet_grid(. ~ facility, scales = 'free') +
  theme_classic() +
  theme(panel.border = element_rect(color="black", fill = NA),
        text = element_text(family = 'Times New Roman', size = 18)) +
  labs(x = 'Decomposition Stage',
       y = '|bNTI| Value')

# Just deterministic
combined 
  

# Barplot of percenteages relative abundance
categorical <- combined %>% 
  mutate(deterministic = case_when(
           value >= 2 ~ 'Heterogenous Selection',
           value <= -2 ~ 'Homogenous Selection',
           value < 2 & value > -2 ~ 'Stochastic Assembly')) %>% 
  filter(deterministic != 'NA')
 
categorical %>% 
  group_by(decomp_stage, facility, deterministic) %>% 
  dplyr::summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>% 
  ungroup() %>% 
  ggplot(aes(x = decomp_stage, y = freq, fill = deterministic)) +
  geom_bar(stat = 'identity',  alpha = 0.3) +
  facet_grid(. ~ facility, scales = 'free') +
  scale_fill_manual(values = c("#FCAE12", "#380F63", '#B93556')) +
  theme_classic() +
  theme(panel.border = element_rect(color="black", fill = NA),
        text = element_text(family = 'Times New Roman', size = 18),
        legend.position = 'bottom') +
  labs(x = 'Decomposition Stage',
       y = 'Relative Abundance',
       fill = 'Community Assembly\nMethod') +
  scale_y_continuous(labels=scales::percent)
           





  
  

  
  
  
  
  
  
  
  
  
##
##### Comparing 16S and Metagenomics #####
### Alpha ###
alpha <- read_tsv('03_metadata/combined-metadata-simple-oct2020_R.txt')

#faiths
cor.test(alpha$faiths_16S, alpha$faiths_meta,
         method = 'pearson', conf.level = 0.95)

alpha %>% 
  filter(decomp_stage != 'not_applicable') %>% 
  mutate(decomp_stage = fct_relevel(decomp_stage, 'early', 'active')) %>% 
  ggplot(aes(x = faiths_16S, y = faiths_meta)) +
  geom_point(aes(color = decomp_stage, shape = facility), size = 2) +
  geom_text(label = 'r^2 == 0.78', x = 250, y = 125, 
            color = 'black', parse = TRUE) +
  geom_smooth(method = lm, se = FALSE, color = 'black') +
  theme_classic() +
  labs(x = "16S Faith's Diversity",
       y = "Metagenomics Faith's Diversity",
       color = "Decomposition Stage",
       shape = "Facility") +
  scale_color_viridis(discrete =  TRUE, end = 0.8)

#richness
cor.test(alpha$richness_16S, alpha$richness_meta,
         method = 'pearson', conf.level = 0.95)

alpha %>% 
  filter(decomp_stage != 'not_applicable') %>% 
  mutate(decomp_stage = fct_relevel(decomp_stage, 'early', 'active')) %>% 
  ggplot(aes(x = richness_16S, y = richness_meta)) +
  geom_point(aes(color = decomp_stage, shape = facility), size = 2) +
  geom_text(label = 'r^2 == 0.6989', x = 4000, y = 1000, 
            color = 'black', parse = TRUE) +
  geom_smooth(method = lm, se = FALSE, color = 'black') +
  theme_classic() +
  labs(x = "16S Richness",
       y = "Metagenomics Richness",
       color = "Decomposition Stage",
       shape = "Facility") +
  scale_color_viridis(discrete =  TRUE, end = 0.8)

### Taxa Rel Abundances - Class###
# Import Data
amplicon_class <- read_tsv('08_results/rel-freq-tables/class_relfreq_table.tsv')
meta_class <- read_tsv('../../shotgun/qiime2-2020-8/WoL/feature_tables/final-class-relfreq-table.tsv')

metadata <- read_tsv('03_metadata/combined-metadata-simple-oct2020_R.txt') %>% 
  mutate(decomp_stage = as.factor(decomp_stage),
         facility = as.factor(facility))

# Manage amplicon class data so it can be merged
amplicon_class %>% 
  gather(key = 'taxon', value = 'amplicon_abundance', -sample) %>% 
  mutate(taxon = gsub('D_0__Bacteria;', '', taxon),
         taxon = gsub('D_0__Archaea;', '', taxon),
         taxon = gsub('D_1__', ' ', taxon),
         taxon = gsub('D_2__', ' ', taxon)) %>% 
  unite(sample_tax, sample, taxon, sep = ' --')  -> amplicon_class2

# Metge data, create data frame of both 16S and shotgun
meta_class %>% 
  gather(key = 'taxon', value = 'meta_abundance', -sample) %>% 
  mutate(taxon = gsub('k__Bacteria;', '', taxon),
         taxon = gsub('k__Archaea;', '', taxon),
         taxon = gsub('p__', ' ', taxon),
         taxon = gsub('c__', ' ', taxon)) %>% 
  unite(sample_tax, sample, taxon, sep = ' --', remove = FALSE) %>% 
  inner_join(amplicon_class2, by = 'sample_tax') -> full_class

# Run grouped correlation test, filter non-significant and low correlation
corr_test <- full_class %>% 
  group_by(taxon) %>% 
  summarize(corr = cor.test(amplicon_abundance, meta_abundance)$estimate, 
            p = cor.test(amplicon_abundance, meta_abundance)$p.value) 
  filter(p <= 0.05 & corr >= 0.4)

# Final clean-up of table, make plot
full_class %>% 
  left_join(corr_test, by = 'taxon') %>% 
  left_join(metadata, by = 'sample') %>%
  #filter(p <= 0.05 & corr >= 0.4) %>% 
  mutate(decomp_stage = fct_relevel(decomp_stage, 'early', 'active')) %>% 
  ggplot(aes(x = amplicon_abundance, y = meta_abundance)) +
  geom_point(aes(color = decomp_stage, shape = facility)) +
  geom_smooth(method = lm, se = FALSE, color = 'black') +
  facet_wrap(.~ taxon, scales = 'free') +
  theme_classic() +
  labs(x = "16S Relative Abundance",
       y = "Metagenomics Relative Abundance",
       color = "Decomposition Stage",
       shape = "Facility") +
  scale_color_viridis(discrete =  TRUE, end = 0.8) +
  geom_text(data = corr_test, aes(label = round(corr, 4), 
                                  x = -Inf, y = Inf - 1,
                                  hjust = 0, vjust = 1))

### Taxa Rel Abundances - Phylum ###
# Import Data
amplicon_phylum <- read_tsv('08_results/rel-freq-tables/phylum_relfreq_table.tsv')
meta_phylum <- read_tsv('../../shotgun/qiime2-2020-8/WoL/feature_tables/final-phylum-relfreq-table.tsv')

metadata <- read_tsv('03_metadata/combined-metadata-simple-oct2020_R.txt') %>% 
  mutate(decomp_stage = as.factor(decomp_stage),
         facility = as.factor(facility))

# Manage amplicon class data so it can be merged
amplicon_phylum %>% 
  gather(key = 'taxon', value = 'amplicon_abundance', -sample) %>% 
  mutate(taxon = gsub('D_0__Bacteria;', '', taxon),
         taxon = gsub('D_0__Archaea;', '', taxon),
         taxon = gsub('D_1__', ' ', taxon),
         taxon = gsub('D_2__', ' ', taxon)) %>% 
  unite(sample_tax, sample, taxon, sep = ' --')  -> amplicon_phylum2

# Metge data, create data frame of both 16S and shotgun
meta_phylum %>% 
  gather(key = 'taxon', value = 'meta_abundance', -sample) %>% 
  mutate(taxon = gsub('k__Bacteria;', '', taxon),
         taxon = gsub('k__Archaea;', '', taxon),
         taxon = gsub('p__', ' ', taxon),
         taxon = gsub('c__', ' ', taxon)) %>% 
  unite(sample_tax, sample, taxon, sep = ' --', remove = FALSE) %>% 
  inner_join(amplicon_phylum2, by = 'sample_tax') %>% 
  filter(taxon == ' Acidobacteria' | taxon == ' Actinobacteria' |
         taxon == ' Bacteroidetes' | taxon == ' Chloroflexi' |
         taxon == ' Cyanobacteria' | taxon == ' Firmicutes' |
         taxon == ' Proteobacteria' | taxon == ' Verrucomicrobia') -> full_phylum

# Run grouped correlation test, filter non-significant and low correlation
corr_test_phylum <- full_phylum %>% 
  group_by(taxon) %>% 
  summarize(corr = cor.test(amplicon_abundance, meta_abundance)$estimate, 
            p = cor.test(amplicon_abundance, meta_abundance)$p.value)
#  filter(p <= 0.05 & corr >= 0.4)

# Final clean-up of table, make plot
full_phylum %>% 
  left_join(corr_test_phylum, by = 'taxon') %>% 
  left_join(metadata, by = 'sample') %>%
  #filter(p <= 0.05 & corr >= 0.4) %>% 
  mutate(decomp_stage = fct_relevel(decomp_stage, 'early', 'active')) %>% 
  ggplot(aes(x = amplicon_abundance, y = meta_abundance)) +
  geom_point(aes(color = decomp_stage, shape = facility)) +
  geom_smooth(method = lm, se = FALSE, color = 'black') +
  facet_wrap(.~ taxon, scales = 'free', ncol = 2) +
  theme_classic() +
  labs(x = "16S Relative Abundance",
       y = "Metagenomics Relative Abundance",
       color = "Decomposition Stage",
       shape = "Facility") +
  scale_color_viridis(discrete =  TRUE, end = 0.8) +
  geom_text(data = corr_test_phylum, aes(label = round(corr, 4), 
                                  x = -Inf, y = Inf - 1,
                                  hjust = 0, vjust = 1)) +
  theme(legend.position = 'bottom') +
  guides(colour = guide_legend(title.position = 'top', ncol = 1),
         shape = guide_legend(title.position = 'top', ncol = 1)) -> taxa_comp_plot

ggsave('taxa_comp_plot.png', plot = taxa_comp_plot, path = '../figures/',
       device = 'png', width = 4, height = 8, units = "in")

### Comparing Beta Div ###
#https://rdrr.io/rforge/vegan/man/mantel.html


amplicon_beta <- read_qza('06_core-metrics/filtered-core-metrics-5000/compare-with-shotgun/generalized-unifrac05-shotgun-samples-distance-matrix.qza')
shotgun_beta <- read_qza('../../shotgun/qiime2-2020-8/WoL/core-metrics-results-all-9039/beta_diversity/matrices/generalized-unifrac05-distance-matrix.qza')

mantel(amplicon_beta$data, shotgun_beta$data)

##
##### Model errors #####
errors <- read_excel('08_results/MAE-tables-facility-3-27-19.xlsx', sheet = 'ASV')
colnames(errors)

errors %>% 
  mutate(type = replace(type, which(model == 'General'), 'general'),
         type = fct_relevel(type, levels = c('within', 'cross-facility', 'general'))) %>% 
  ggplot(aes(x = model, y = error, fill = type, group = factor(model_facility))) +
  geom_bar(stat = "identity", color = 'black') +
  geom_text(aes(label = error), color = "white", size = 3, vjust = -0.1) +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.8) +
  scale_x_discrete(limits = c("CMU", "SHSU", "UTK",
                              "CMU on SHSU","CMU on UTK", "SHSU on CMU",
                              "SHSU on UTK", "UTK on CMU", "UTK on SHSU", 
                              "General"))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95,
                                   color = 'white'),
        axis.text.y = element_text(color = 'white'),
        axis.title = element_text(color = 'white'),
        plot.background = element_rect(fill = 'black'),
        legend.text = element_text(color = 'white'),
        legend.title = element_text(color = 'white'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'right') +
  labs(x = "Model Facility" , y = "Mean Absolute Error (Days)", 
       fill = "Model Type") +
  scale_y_continuous(limits = c(0,15), expand = c(0, 0))


## White for proposal
errors %>% 
  mutate(type = replace(type, which(model == 'General'), 'general'),
         type = fct_relevel(type, levels = c('within', 'cross-facility', 'general'))) %>% 
  ggplot(aes(x = model, y = error, fill = type, group = factor(model_facility))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = error), color = "black", size = 3.5, vjust = -0.1) +
  scale_x_discrete(limits = c("CMU", "SHSU", "UTK",
                              "CMU on SHSU","CMU on UTK", "SHSU on CMU",
                              "SHSU on UTK", "UTK on CMU", "UTK on SHSU", 
                              "General"))+
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.8) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95,
                                   color = 'black', size = 10),
        axis.text.y = element_text(color = 'black', size = 10),
        axis.title = element_text(color = 'black', size = 10),
        plot.background = element_rect(fill = 'white'),
        legend.text = element_text(color = 'black', size = 10),
        legend.title = element_text(color = 'black', size = 10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'right') +
  labs(x = "Model Facility" , y = "Mean Absolute Error (Days)", 
       fill = "Model Type") +
  scale_y_continuous(limits = c(0,15), expand = c(0, 0))

ggsave('08_results/allseasons-ASV-errors2.png', plot = white_errors, height = 4, width = 6)

## Just SHSU for proposal
SHSU_errors <- read_excel('08_results/MAE-tables-SHSU-seasons-3-28-19.xlsx', sheet = 'ASV')

SHSU_errors %>% 
  mutate(type = replace(type, which(model == 'General'), 'general'),
         type = fct_relevel(type, levels = c('within', 'cross-facility', 'general'))) %>% 
  ggplot(aes(x = model, y = error, fill = type, group = factor(model_season))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = error), color = "black", size = 3, vjust = -0.1) +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.8) +
  scale_x_discrete(limits = c("Spring", "Summer", "Fall", "Winter",
                              "Spring on Summer", "Spring on Fall", "Spring on Winter",
                              "Summer on Spring", "Summer on Fall", "Summer on Winter",
                              "Fall on Spring", "Fall on Summer", "Fall on Winter",
                              "Winter on Spring", "Winter on Summer", "Winter on Fall", 
                              "General"))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95,
                                   color = 'black'),
        axis.text.y = element_text(color = 'black'),
        axis.title = element_text(color = 'black'),
        plot.background = element_rect(fill = 'white'),
        legend.text = element_text(color = 'black'),
        legend.title = element_text(color = 'black'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'right') +
  labs(x = "Model Facility" , y = "Mean Absolute Error (Days)", 
       fill = "Model Type") +
  scale_y_continuous(limits = c(0,16.5), expand = c(0, 0))

SHSUboth <- SHSU_errors %>% 
  filter(type == 'within') %>% 
  mutate(type = replace(type, which(model == 'General'), 'general'),
         type = fct_relevel(type, levels = c('within', 'cross-facility', 'general'))) %>% 
  ggplot(aes(x = model, y = error, fill = type, group = factor(model_season))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = error), color = "black", size = 3.5, vjust = -0.1) +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.8) +
  scale_x_discrete(limits = c("Spring","Summer","Fall","Winter","General")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95,
                                   color = 'black', size = 15),
        axis.text.y = element_text(color = 'black', size = 15),
        axis.title = element_text(color = 'black', size = 15),
        plot.background = element_rect(fill = 'white'),
        legend.text = element_text(color = 'black', size = 10),
        legend.title = element_text(color = 'black', size = 10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'right') +
  labs(x = "Model Facility" , y = "Mean Absolute Error (Days)", 
       fill = "Model Type") +
  scale_y_continuous(limits = c(0,10), expand = c(0, 0))

SHSUboth
ggsave('08_results/allseasonsSHSU-ASV-errors.png', plot = SHSUboth, height = 4, width = 6)

## Within facility by season for zach AAFS
max_add <- read_excel('08_results/MAE-tables-forAAFS2020.xlsx', sheet = 'max_add')

read_excel('08_results/MAE-tables-forAAFS2020.xlsx', sheet = 'combined_errors') %>% 
  mutate(type = replace(type, which(model == 'General'), 'general'),
         type = fct_relevel(type, levels = c('within', 'cross-facility', 'general')),
         error = round(error, 2),
         season = fct_relevel(season, 'Spring', 'Summer', 'Fall')) %>% 
  ggplot(aes(x = model, y = error, fill = type, group = factor(model_facility))) +
  geom_bar(stat = "identity", color = 'black') +
  geom_text(aes(label = error), color = "white", size = 3, vjust = -0.1) +
  #geom_point(aes(x = model, y = max_add), 
  #           color = 'darkgray', shape = 6) +
  facet_wrap(. ~ season) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95,
                                   color = 'white'),
        axis.text.y = element_text(color = 'white'),
        axis.title = element_text(color = 'white'),
        plot.background = element_rect(fill = 'black'),
        legend.text = element_text(color = 'white'),
        legend.title = element_text(color = 'white'),
        #strip.background = element_rect(color = 'white'),
        strip.text = element_text(color = 'white'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'right') +
  labs(x = "Model Facility" , y = "Mean Absolute Error (ADD)", 
       fill = "Model Type") +
  guides(shape = FALSE) +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.8) +
  scale_x_discrete(limits = c("CMU", "SHSU", "UTK",
                              "CMU on SHSU","CMU on UTK", "SHSU on CMU",
                              "SHSU on UTK", "UTK on CMU", "UTK on SHSU", 
                              "General")) +
  scale_y_continuous(limits = c(0,300), expand = c(0, 0))


##### Model errors combined 16S and 18S#####
errorsboth <- read_excel('08_results/16S18S-MAE_tables_facility_2-13-19.xlsx', sheet = 'ASV')
colnames(errors)

errorsboth %>% 
  mutate(type = replace(type, which(model == 'General'), 'general'),
         type = fct_relevel(type, levels = c('within', 'cross-facility', 'general'))) %>%
  ggplot(aes(x = model, y = error, fill = type, group = factor(model_facility))) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.8) +
  scale_x_discrete(limits = c("CMU", "SHSU", "UTK",
                              "CMU on SHSU","CMU on UTK", "SHSU on CMU",
                              "SHSU on UTK", "UTK on CMU", "UTK on SHSU", 
                              "General"))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95,
                                   color = 'white'),
        axis.text.y = element_text(color = 'white'),
        axis.title = element_text(color = 'white'),
        plot.background = element_rect(fill = 'black'),
        legend.text = element_text(color = 'white'),
        legend.title = element_text(color = 'white'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'right') +
  labs(x = "Model Facility" , y = "Mean Absolute Error (Days)", 
       fill = "Model Type") +
  scale_y_continuous(limits = c(0,7.5), expand = c(0, 0))

##### Model errors to compare amplicons #####
## Comparing just 16S, just 18S, both

error16S <- read_excel('/Users/Aeriel/Dropbox/projects/PMI_3_analyses/deblur_16S/deblur_16S_spring/qiime2_analysis/results/MAE_tables_facility_4-21-18.xlsx', sheet = 'Sheet1') %>% 
  mutate(amplicon = '16S',
         model_amplicon = model_facility) %>% 
  unite(model_amplicon,amplicon, model_amplicon, sep = " ", remove = FALSE)

error18S <- read_excel('08_results/MAE_tables_facility_2-12-19.xlsx', sheet = 'ASV') %>% 
  mutate(amplicon = '18S',
         model_amplicon = model_facility) %>% 
  unite(model_amplicon,amplicon, model_amplicon, sep = " ", remove = FALSE)

errorboth <-  read_excel('08_results/16S18S-MAE_tables_facility_2-13-19.xlsx', sheet = 'ASV') %>% 
  mutate(amplicon = 'Both',
         model_amplicon = model_facility) %>% 
  unite(model_amplicon,amplicon, model_amplicon, sep = " ", remove = FALSE)

fullerror <- full_join(error16S, error18S) %>% 
  full_join(errorboth)
  

fullerror %>% 
  filter(type == "within") %>% 
  ggplot(aes(x = model_amplicon, y = error, fill = amplicon, group = factor(model_facility))) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.8) +
  scale_x_discrete(limits = c("16S General", "18S General", "Both General",
                              "16S CMU", "18S CMU", "Both CMU",
                              "16S SHSU", "18S SHSU", "Both SHSU",
                              "16S UTK", "18S UTK", "Both UTK")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95,
                                   color = 'white'),
        axis.text.y = element_text(color = 'white'),
        axis.title = element_text(color = 'white'),
        plot.background = element_rect(fill = 'black'),
        legend.text = element_text(color = 'white'),
        legend.title = element_text(color = 'white'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
        #legend.position = 'bottom') +
  labs(x = "Model Facility" , y = "Mean Absolute Error (Days)", 
       fill = "Model Type") +
  scale_y_continuous(limits = c(0,5), expand = c(0, 0))


##### Model errors on just most important features #####
highimperrors <- read_excel('08_results/16S18S-highimport-MAE_tables_facility_2-13-19.xlsx', sheet = 'ASV')
colnames(errors)

highimperrors %>% 
  mutate(type = replace(type, which(model == 'General'), 'general'),
         type = fct_relevel(type, levels = c('within', 'cross-facility', 'general'))) %>% 
  ggplot(aes(x = model, y = error, fill = type, group = factor(model_facility))) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +
  scale_fill_viridis(discrete = TRUE, begin = 0.15, end = 0.8) +
  scale_x_discrete(limits = c("CMU", "SHSU", "UTK",
                              "CMU on SHSU","CMU on UTK", "SHSU on CMU",
                              "SHSU on UTK", "UTK on CMU", "UTK on SHSU", 
                              "General"))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95,
                                   color = 'white'),
        axis.text.y = element_text(color = 'white'),
        axis.title = element_text(color = 'white'),
        plot.background = element_rect(fill = 'black'),
        legend.text = element_text(color = 'white'),
        legend.title = element_text(color = 'white'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'right') +
  labs(x = "Model Facility" , y = "Mean Absolute Error (Days)", 
       fill = "Model Type") +
  scale_y_continuous(limits = c(0,7), expand = c(0, 0))


##### Model Errors FINAL ######
errors <- read_tsv('../machine_learning/CV_errors.txt')

errors %>% 
  separate(model, 
           into = c('amplicon', 'sample_type', 'sample_site', 'taxon_level', 
                    'model_lab', 'results'), 
           sep = '_', remove = FALSE) %>% 
  mutate(taxon_level = fct_recode(taxon_level, 'species' = 'model'),
         taxon_level = fct_relevel(taxon_level, 'phylum', 'class', 'order', 
                                   'family', 'genus')) %>% 
  select(-model_lab, -results) %>% 
  unite(model_type, sample_type, sample_site, sep = '.', remove = FALSE) %>% 
  ggplot(aes(x = taxon_level, y = error, fill = sample_type)) +
  geom_boxplot() +
  facet_grid(amplicon ~ model_type) +
  theme_classic() +
  scale_fill_viridis(discrete =  TRUE, end = 0.8) +
  labs(x = 'Taxonomic Level',
       y = "Nested CV MAE Scores (ADD)",
       fill = 'Sample Type')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95),
        panel.border = element_rect(fill = NA, color = "black", size = 1)) -> errors_plot

ggsave('errors_plot.png', plot = errors_plot, path = '../figures/',
       device = 'png', width = 8, height = 6, units = "in")


best_errors <- read_tsv('../machine_learning/gen_errors.txt')

best_errors %>% 
  separate(model, 
           into = c('amplicon', 'sample_type', 'sample_site', 'taxon_level', 
                    'model_lab', 'results'), 
           sep = '_', remove = FALSE) %>% 
  mutate(taxon_level = fct_recode(taxon_level, 'species' = 'model'),
         taxon_level = fct_relevel(taxon_level, 'phylum', 'class', 'order', 
                                   'family', 'genus')) %>% 
  select(-model_lab, -results) %>% 
  unite(model_type, sample_type, sample_site, sep = '.', remove = FALSE) -> best_errors

write_tsv(best_errors,"../machine_learning/best_errors.txt")

best_errors %>% 
  group_by(amplicon, sample_type) %>% 
  summarize(mean = mean(best_error))












