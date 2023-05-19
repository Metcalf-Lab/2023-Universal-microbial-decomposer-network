#!/bin/sh
#SBATCH --job-name=import_18S
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --partition=shas
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL


#Activate qiime
source activate qiime2-2022.8


#1-6
qiime tools import \
  --type EMPSingleEndSequences \
  --input-path 18S_raw_data/emp_se_seqs/4733_NIJ_1-6_18S \
  --output-path 03_artifacts-and-visualizations/completed_18S_sequences/1-6_sequences.qza

#7-12
qiime tools import \
  --type EMPSingleEndSequences \
  --input-path 18S_raw_data/emp_se_seqs/4733_NIJ_7-12_18S \
  --output-path 03_artifacts-and-visualizations/completed_18S_sequences/7-12_sequences.qza

#13-17, 52
qiime tools import \
  --type EMPSingleEndSequences \
  --input-path 18S_raw_data/emp_se_seqs/5465_Metcalf_NIJ_13-17_52 \
  --output-path 03_artifacts-and-visualizations/completed_18S_sequences/13-17_52_18S_sequences.qza
  
#18-21_24-25
qiime tools import \
  --type EMPSingleEndSequences \
  --input-path 18S_raw_data/emp_se_seqs/6398_Metcalf_NIJ_18-21_24_25 \
  --output-path 03_artifacts-and-visualizations/completed_18S_sequences/18-21_24-25_sequences.qza
  
#26-28_30-31_49
qiime tools import \
  --type EMPSingleEndSequences \
  --input-path 18S_raw_data/emp_se_seqs/6353_M_NIJ_1 \
  --output-path 03_artifacts-and-visualizations/completed_18S_sequences/26-28_30-31_49_sequences.qza
  
#29_48_51-53
qiime tools import \
  --type EMPSingleEndSequences \
  --input-path 18S_raw_data/emp_se_seqs/6351_M_NIJ_2 \
  --output-path 03_artifacts-and-visualizations/completed_18S_sequences/29_48_51-53_sequences.qza
  
#32_37-41
qiime tools import \
  --type EMPSingleEndSequences \
  --input-path 18S_raw_data/emp_se_seqs/5403_Metcalf_NIJ_32_37-41 \
  --output-path 03_artifacts-and-visualizations/completed_18S_sequences/32_37-41_sequences.qza

#42-47
qiime tools import \
  --type EMPSingleEndSequences \
  --input-path 18S_raw_data/emp_se_seqs/5402_Metcalf_NIJ_42-47 \
  --output-path 03_artifacts-and-visualizations/completed_18S_sequences/42-47_sequences.qza
