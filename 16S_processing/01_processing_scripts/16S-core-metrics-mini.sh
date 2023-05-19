#!/bin/sh
#SBATCH --job-name=16S-diversity
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --qos=long
#SBATCH --partition=shas
#SBATCH --time=100:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aeriel.belk@gmail.com

#Activate qiime
source activate qiime2-2020.2

#Command
qiime fragment-insertion sepp \
  --i-representative-sequences 16S-deblur-mini/combined-16S-rep-seqs-mini.qza \
  --i-reference-database sepp-refs-gg-13-8.qza \
  --p-threads 5 \
  --o-tree 05_trees/16S_sepp_tree_mini.qza \
  --o-placements 05_trees/16S_sepp_placements_mini.qza 
  
qiime diversity core-metrics-phylogenetic \
  --i-table 16S-deblur-mini/combined-16S-table-mini.qza \
  --i-phylogeny 05_trees/16S_sepp_tree_mini.qza \
  --m-metadata-file 03_metadata/combined-metadata-simple-Jan2020.txt \
  --p-sampling-depth 5000 \
  --output-dir core-metrics-5000-mini
