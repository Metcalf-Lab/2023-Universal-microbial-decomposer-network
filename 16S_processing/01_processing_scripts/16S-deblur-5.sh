#!/bin/sh
#SBATCH --job-name=16S-deblur-5
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --partition=shas
#SBATCH --qos=long
#SBATCH --time=148:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aeriel.belk@gmail.com

#Activate qiime
source activate qiime2-2020.2

#Command
qiime deblur denoise-16S \
  --i-demultiplexed-seqs 04_artifacts-and-visualizations/demux-files/16S-9_demux.qza \
  --p-trim-length 150 \
  --p-jobs-to-start 5 \
  --p-sample-stats \
  --p-min-reads 1 \
  --o-table 04_artifacts-and-visualizations/deblur-files/16S-9_table.qza \
  --o-representative-sequences 04_artifacts-and-visualizations/deblur-files/16S-9_rep-seqs.qza \
  --o-stats 04_artifacts-and-visualizations/deblur-files/16S-9_deblur-stats.qza

qiime deblur denoise-16S \
  --i-demultiplexed-seqs 04_artifacts-and-visualizations/demux-files/16S-10_demux.qza \
  --p-trim-length 150 \
  --p-jobs-to-start 5 \
  --p-sample-stats \
  --p-min-reads 1 \
  --o-table 04_artifacts-and-visualizations/deblur-files/16S-10_table.qza \
  --o-representative-sequences 04_artifacts-and-visualizations/deblur-files/16S-10_rep-seqs.qza \
  --o-stats 04_artifacts-and-visualizations/deblur-files/16S-10_deblur-stats.qza
