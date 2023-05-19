#!/bin/sh
#SBATCH --job-name=pmi3_import
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=ssky
#SBATCH --qos=long
#SBATCH --time=148:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aeriel.belk@gmail.com


#Activate qiime
source activate qiime2-2020.2

#Command
qiime tools import \
  --type EMPSingleEndSequences \
  --input-path 16S-1 \
  --output-path 16S-1_sequences.qza

qiime tools import \
  --type EMPSingleEndSequences \
  --input-path 16S-2 \
  --output-path 16S-2_sequences.qza
  
qiime tools import \
  --type EMPSingleEndSequences \
  --input-path 16S-3 \
  --output-path 16S-3_sequences.qza

qiime tools import \
  --type EMPSingleEndSequences \
  --input-path 16S-4 \
  --output-path 16S-4_sequences.qza

qiime tools import \
  --type EMPSingleEndSequences \
  --input-path 16S-5 \
  --output-path 16S-5_sequences.qza

qiime tools import \
  --type EMPSingleEndSequences \
  --input-path 16S-6 \
  --output-path 16S-6_sequences.qza

qiime tools import \
  --type EMPSingleEndSequences \
  --input-path 16S-7 \
  --output-path 16S-7_sequences.qza

qiime tools import \
  --type EMPSingleEndSequences \
  --input-path 16S-8 \
  --output-path 16S-8_sequences.qza

qiime tools import \
  --type EMPSingleEndSequences \
  --input-path 16S-9 \
  --output-path 16S-9_sequences.qza
