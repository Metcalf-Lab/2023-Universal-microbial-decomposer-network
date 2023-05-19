#!/bin/sh
#SBATCH --job-name=richness
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=shas
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aeriel.belk@gmail.com

#Activate qiime
source activate qiime2-2020.2

#Command
qiime feature-table rarefy \
  --i-table 04_artifacts-and-visualizations/PMI-16S-table.qza \
  --p-sampling-depth 4548 \
  --o-rarefied-table 04_artifacts-and-visualizations/PMI-16S-4548-unfiltered-table.qza

qiime diversity alpha \
  --i-table 04_artifacts-and-visualizations/PMI-16S-4548-unfiltered-table.qza \
  --p-metric 'observed_otus' \
  --o-alpha-diversity PMI-unfiltered-full-4548-richness.qza

qiime diversity alpha-group-significance \
  --i-alpha-diversity PMI-unfiltered-full-4548-richness.qza \
  --m-metadata-file 03_metadata/combined-metadata-simple-oct2020.txt \
  --o-visualization PMI-unfiltered-full-4548-richness.qzv
