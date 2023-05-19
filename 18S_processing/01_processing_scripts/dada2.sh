#!/bin/sh
#SBATCH --job-name=dada2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=shas
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL


#Activate qiime
source activate qiime2-2020.8


#Command
#1-6
qiime dada2 denoise-single \
  --i-demultiplexed-seqs cutadapt/4733_NIJ_1-6_18S_noprimer_seqs.qza \
  --p-trunc-len 0 \
  --o-table dada2/4733_NIJ_1-6_18S_noprimer_table.qza \
  --o-representative-sequences dada2/4733_NIJ_1-6_18S_noprimer_rep-seqs.qza  \
  --o-denoising-stats dada2/4733_NIJ_1-6_18S_noprimer_stats.qza \
  --p-n-threads 24

qiime metadata tabulate \
  --m-input-file dada2/4733_NIJ_1-6_18S_noprimer_stats.qza \
  --o-visualization dada2/4733_NIJ_1-6_18S_noprimer_stats.qzv

#7-12
qiime dada2 denoise-single \
  --i-demultiplexed-seqs cutadapt/4733_NIJ_7-12_18S_noprimer_seqs.qza \
  --p-trunc-len 0 \
  --o-table dada2/4733_NIJ_7-12_18S_noprimer_table.qza \
  --o-representative-sequences dada2/4733_NIJ_7-12_18S_noprimer_rep-seqs.qza  \
  --o-denoising-stats dada2/4733_NIJ_7-12_18S_noprimer_stats.qza \
  --p-n-threads 24

qiime metadata tabulate \
  --m-input-file dada2/4733_NIJ_7-12_18S_noprimer_stats.qza \
  --o-visualization dada2/4733_NIJ_7-12_18S_noprimer_stats.qzv

#13-17_52
qiime dada2 denoise-single \
  --i-demultiplexed-seqs cutadapt/5465_NIJ_13-17_52_18S_noprimer_seqs.qza \
  --p-trunc-len 0 \
  --o-table dada2/5465_NIJ_13-17_52_18S_noprimer_table.qza \
  --o-representative-sequences dada2/5465_NIJ_13-17_52_18S_noprimer_rep-seqs.qza  \
  --o-denoising-stats dada2/5465_NIJ_13-17_52_18S_noprimer_stats.qza \
  --p-n-threads 5 \
  --verbose
 
qiime metadata tabulate \
  --m-input-file dada2/5465_NIJ_13-17_52_18S_noprimer_stats.qza \
  --o-visualization dada2/5465_NIJ_13-17_52_18S_noprimer_stats.qzv

#18-21_24-25
qiime dada2 denoise-single \
  --i-demultiplexed-seqs cutadapt/6398_NIJ_18-21_24_25_18S_noprimer_seqs.qza \
  --p-trunc-len 0 \
  --o-table dada2/6398_NIJ_18-21_24_25_18S_noprimer_table.qza \
  --o-representative-sequences dada2/6398_NIJ_18-21_24_25_18S_noprimer_rep-seqs.qza  \
  --o-denoising-stats dada2/6398_NIJ_18-21_24_25_18S_noprimer_stats.qza \
  --p-n-threads 5 \
  --verbose

qiime metadata tabulate \
  --m-input-file dada2/6398_NIJ_18-21_24_25_18S_noprimer_stats.qza \
  --o-visualization dada2/6398_NIJ_18-21_24_25_18S_noprimer_stats.qzv

#6353_NIJ_1
qiime dada2 denoise-single \
  --i-demultiplexed-seqs cutadapt/6353_NIJ_1_18S_noprimer_seqs.qza \
  --p-trunc-len 0 \
  --o-table dada2/6353_NIJ_1_18S_noprimer_table.qza \
  --o-representative-sequences dada2/6353_NIJ_1_18S_noprimer_rep-seqs.qza  \
  --o-denoising-stats dada2/6353_NIJ_1_18S_noprimer_stats.qza \
  --p-n-threads 24

qiime metadata tabulate \
  --m-input-file dada2/6353_NIJ_1_18S_noprimer_stats.qza \
  --o-visualization dada2/6353_NIJ_1_18S_noprimer_stats.qzv

#6351_NIJ_2
qiime dada2 denoise-single \
  --i-demultiplexed-seqs cutadapt/6351_NIJ_2_18S_noprimer_seqs.qza \
  --p-trunc-len 0 \
  --o-table dada2/6351_NIJ_2_18S_noprimer_table.qza \
  --o-representative-sequences dada2/6351_NIJ_2_18S_noprimer_rep-seqs.qza  \
  --o-denoising-stats dada2/6351_NIJ_2_18S_noprimer_stats.qza \
  --p-n-threads 24

qiime metadata tabulate \
  --m-input-file dada2/6351_NIJ_2_18S_noprimer_stats.qza \
  --o-visualization dada2/6351_NIJ_2_18S_noprimer_stats.qzv

#32_37-41
qiime dada2 denoise-single \
  --i-demultiplexed-seqs cutadapt/5403_NIJ_32_37-41_18S_noprimer_seqs.qza \
  --p-trunc-len 0 \
  --o-table dada2/5403_NIJ_32_37-41_18S_noprimer_table.qza \
  --o-representative-sequences dada2/5403_NIJ_32_37-41_18S_noprimer_rep-seqs.qza  \
  --o-denoising-stats dada2/5403_NIJ_32_37-41_18S_noprimer_stats.qza \
  --p-n-threads 5 \
  --verbose

qiime metadata tabulate \
  --m-input-file dada2/5403_NIJ_32_37-41_18S_noprimer_stats.qza \
  --o-visualization dada2/5403_NIJ_32_37-41_18S_noprimer_stats.qzv

#42-47
qiime dada2 denoise-single \
  --i-demultiplexed-seqs cutadapt/5402_Metcalf_NIJ_42-47_18S_noprimer_seqs.qza \
  --p-trunc-len 0 \
  --o-table dada2/5402_Metcalf_NIJ_42-47_18S_noprimer_table.qza \
  --o-representative-sequences dada2/5402_Metcalf_NIJ_42-47_18S_noprimer_rep-seqs.qza  \
  --o-denoising-stats dada2/5402_Metcalf_NIJ_42-47_18S_noprimer_stats.qza \
  --p-n-threads 24

qiime metadata tabulate \
  --m-input-file dada2/5402_Metcalf_NIJ_42-47_18S_noprimer_stats.qza \
  --o-visualization dada2/5402_Metcalf_NIJ_42-47_18S_noprimer_stats.qzv

# Merge 18S tables
cd dada2/

qiime feature-table merge \
  --i-tables 5402_Metcalf_NIJ_42-47_18S_noprimer_table.qza 5403_NIJ_32_37-41_18S_noprimer_table.qza 6351_NIJ_2_18S_noprimer_table.qza 6353_NIJ_1_18S_noprimer_table.qza 4733_NIJ_1-6_18S_noprimer_table.qza 4733_NIJ_7-12_18S_noprimer_table.qza 5465_NIJ_13-17_52_18S_noprimer_table.qza 6398_NIJ_18-21_24_25_18S_noprimer_table.qza \
  --o-merged-table 18S_noprimer_complete_table.qza

qiime feature-table merge-seqs \
  --i-data 4733_NIJ_1-6_18S_noprimer_rep-seqs.qza 4733_NIJ_7-12_18S_noprimer_rep-seqs.qza 5465_NIJ_13-17_52_18S_noprimer_rep-seqs.qza 6398_NIJ_18-21_24_25_18S_noprimer_rep-seqs.qza 6353_NIJ_1_18S_noprimer_rep-seqs.qza 6351_NIJ_2_18S_noprimer_rep-seqs.qza 5403_NIJ_32_37-41_18S_noprimer_rep-seqs.qza 5402_Metcalf_NIJ_42-47_18S_noprimer_rep-seqs.qza \
  --o-merged-data 18S_noprimer_complete_rep-seqs.qza

qiime feature-table summarize \
  --i-table 18S_noprimer_complete_table.qza \
  --o-visualization 18S_noprimer_complete_table.qzv \
  --m-sample-metadata-file ../combined-metadata-simple-nov2020.txt

qiime feature-table tabulate-seqs \
  --i-data 18S_noprimer_complete_rep-seqs.qza \
  --o-visualization 18S_noprimer_complete_rep-seqs.qzv
