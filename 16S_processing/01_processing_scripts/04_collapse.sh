#!/bin/sh
#SBATCH --job-name=16S_collapse_species
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=shas
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aeriel.belk@gmail.com

#Activate qiime
source activate qiime2-2020.2

#Command
qiime taxa collapse \
  --i-table 06_core-metrics/filtered-core-metrics-5000/rarefied_table.qza \
  --i-taxonomy 07_taxonomy/pmi3-16S-nochlomito-silva-taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table 04_artifacts-and-visualizations/rarefied_species_table.qza

qiime tools export \
  --input-path 04_artifacts-and-visualizations/rarefied_species_table.qza \
  --output-path 04_artifacts-and-visualizations/rarefied_species_exported-feature-table

mv 04_artifacts-and-visualizations/rarefied_species_exported-feature-table/feature-table.biom 04_artifacts-and-visualizations/species-normalized-table.biom

rm -r 04_artifacts-and-visualizations/rarefied_species_exported-feature-table
