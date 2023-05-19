#!/bin/sh
#SBATCH --job-name=cutadapt
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --partition=shas
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL


#Activate qiime
source activate qiime2-2022.8

#Command
#1-6
qiime cutadapt trim-single \
  --i-demultiplexed-sequences 1-6_18S_demux.qza \
  --p-adapter GTAGGTGAACCTGCAGAAGGATCA \
  --o-trimmed-sequences 1-6_18S_noprimer_sequences.qza \
  --p-cores 24
  
#7-12
qiime cutadapt trim-single \
  --i-demultiplexed-sequences 7-12_18S_demux.qza \
  --p-adapter GTAGGTGAACCTGCAGAAGGATCA \
  --o-trimmed-sequences 7-12_18S_noprimer_sequences.qza \
  --p-cores 24

#13-17_52
qiime cutadapt trim-single \
  --i-demultiplexed-sequences 13-17_52_18S_demux.qza \
  --p-adapter GTAGGTGAACCTGCAGAAGGATCA \
  --o-trimmed-sequences 13-17_52_18S_noprimer_sequences.qza \
  --p-cores 24

#18-21_24-25
qiime cutadapt trim-single \
  --i-demultiplexed-sequences 18-21_24-25_18S_demux.qza \
  --p-adapter GTAGGTGAACCTGCAGAAGGATCA \
  --o-trimmed-sequences 18-21_24-25_18S_noprimer_sequences.qza \
  --p-cores 24

#26-28_30-31
qiime cutadapt trim-single \
  --i-demultiplexed-sequences 26-28_30-31_18S_demux.qza \
  --p-adapter GTAGGTGAACCTGCAGAAGGATCA \
  --o-trimmed-sequences 26-28_30-31_18S_noprimer_sequences.qza \
  --p-cores 24

#29_48_51-53
qiime cutadapt trim-single \
  --i-demultiplexed-sequences 29_48_51-53_18S_demux.qza \
  --p-adapter GTAGGTGAACCTGCAGAAGGATCA \
  --o-trimmed-sequences 29_48_51-53_18S_noprimer_sequences.qza \
  --p-cores 24
  
#32_37-41
qiime cutadapt trim-single \
  --i-demultiplexed-sequences 32_37-41_18S_demux.qza \
  --p-adapter GTAGGTGAACCTGCAGAAGGATCA \
  --o-trimmed-sequences 32_37-41_18S_noprimer_sequences.qza \
  --p-cores 24

#42-47
qiime cutadapt trim-single \
  --i-demultiplexed-sequences 42-47_18S_demux.qza \
  --p-adapter GTAGGTGAACCTGCAGAAGGATCA \
  --o-trimmed-sequences 42-47_18S_noprimer_sequences.qza \
  --p-cores 24


