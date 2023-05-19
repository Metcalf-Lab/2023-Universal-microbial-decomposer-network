#!/bin/sh
#SBATCH --job-name=pmi3_import
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=shas
#SBATCH --qos=long
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aeriel.belk@gmail.com


#Activate qiime
source activate qiime2-2020.2

#Command
cd testing-spring-issue

qiime tools import \
  --input-path 58343_all.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path 58343_table.qza

qiime tools import \
  --input-path 66413_all.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path 66413_table.qza

qiime tools import \
  --input-path 66411_all.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path 66411_table.qza

qiime tools import \
  --input-path 66371_all.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path 66371_table.qza

qiime tools import \
  --input-path 66405_all.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path 66405_table.qza

qiime tools import \
  --input-path 66405_all.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path 66405_table.qza
  
qiime tools import \
  --input-path 66420_all.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path 66420_table.qza

qiime tools import \
  --input-path 66415_all.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path 66415_table.qza

qiime tools import \
  --input-path 66360_all.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path 66360_table.qza
  
# merge
qiime feature-table merge \
  --i-tables 58343_table.qza 66413_table.qza 66411_table.qza 66371_table.qza 66405_table.qza 66405_table.qza 66420_table.qza 66415_table.qza 66360_table.qza \
  --p-overlap-method sum \
  --o-merged-table qiita-all-biom-test.qza
