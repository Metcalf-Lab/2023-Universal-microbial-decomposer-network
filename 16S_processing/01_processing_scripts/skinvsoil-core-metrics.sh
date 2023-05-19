#!/bin/sh
#SBATCH --job-name=skinvsoil
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=shas
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aeriel.belk@gmail.com

#Activate qiime
source activate qiime2-2020.2

#Command
qiime feature-table filter-samples \
  --i-table 04_artifacts-and-visualizations/combined-16S-table-filtered-noChloMito.qza \
  --m-metadata-file 03_metadata/combined-metadata-simple-Jan2020.txt \
  --p-where "sample_type='skin'" \
  --o-filtered-table 04_artifacts-and-visualizations/skin-combined-16S-table-filtered-noChloMito.qza

qiime feature-table filter-samples \
  --i-table 04_artifacts-and-visualizations/combined-16S-table-filtered-noChloMito.qza \
  --m-metadata-file 03_metadata/combined-metadata-simple-Jan2020.txt \
  --p-where "sample_type='soil'" \
  --o-filtered-table 04_artifacts-and-visualizations/soil-combined-16S-table-filtered-noChloMito.qza

qiime diversity core-metrics-phylogenetic \
  --i-table 04_artifacts-and-visualizations/skin-combined-16S-table-filtered-noChloMito.qza \
  --i-phylogeny 05_trees/pmi3-16S-sepp-tree.qza \
  --m-metadata-file 03_metadata/combined-metadata-simple-Jan2020.txt \
  --p-sampling-depth 5000 \
  --output-dir 06_core-metrics/skin-5000-sepp-results

qiime diversity core-metrics-phylogenetic \
  --i-table 04_artifacts-and-visualizations/soil-combined-16S-table-filtered-noChloMito.qza \
  --i-phylogeny 05_trees/pmi3-16S-sepp-tree.qza \
  --m-metadata-file 03_metadata/combined-metadata-simple-Jan2020.txt \
  --p-sampling-depth 5000 \
  --output-dir 06_core-metrics/soil-5000-sepp-results

qiime diversity alpha-group-significance \
  --i-alpha-diversity 06_core-metrics/skin-5000-sepp-results/shannon_vector.qza \
  --m-metadata-file 03_metadata/combined-metadata-simple-Jan2020.txt \
  --o-visualization 06_core-metrics/skin-5000-sepp-results/shannon_group_significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity 06_core-metrics/skin-5000-sepp-results/observed_otus_vector.qza \
  --m-metadata-file 03_metadata/combined-metadata-simple-Jan2020.txt \
  --o-visualization 06_core-metrics/skin-5000-sepp-results/observed_otus_group_significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity 06_core-metrics/soil-5000-sepp-results/shannon_vector.qza \
  --m-metadata-file 03_metadata/combined-metadata-simple-Jan2020.txt \
  --o-visualization 06_core-metrics/soil-5000-sepp-results/shannon_group_significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity 06_core-metrics/soil-5000-sepp-results/observed_otus_vector.qza \
  --m-metadata-file 03_metadata/combined-metadata-simple-Jan2020.txt \
  --o-visualization 06_core-metrics/soil-5000-sepp-results/observed_otus_group_significance.qzv
