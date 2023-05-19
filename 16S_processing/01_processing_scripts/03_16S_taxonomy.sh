#!/bin/sh
#SBATCH --job-name=16S_taxonomy
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=smem
#SBATCH --qos=long
#SBATCH --time=100:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aeriel.belk@gmail.com

#Activate qiime
source activate qiime2-2020.2

#Command
qiime feature-classifier classify-sklearn \
  --i-reads 04_artifacts-and-visualizations/PMI-16S-rep-seqs.qza \
  --i-classifier silva-132-99-515-806-nb-classifier.qza \
  --p-n-jobs 4 \
  --o-classification 07_taxonomy/pmi3-16S-silva-taxonomy.qza

qiime metadata tabulate \
  --m-input-file 07_taxonomy/pmi3-16S-silva-taxonomy.qza \
  --o-visualization 07_taxonomy/pmi3-16S-silva-taxonomy.qzv
