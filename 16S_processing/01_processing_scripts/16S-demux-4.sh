#!/bin/sh
#SBATCH --job-name=16S-demux-4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=smem
#SBATCH --qos=long
#SBATCH --time=148:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aeriel.belk@gmail.com

#Activate qiime
source activate qiime2-2020.2

#Command
qiime demux emp-single \
  --i-seqs 16S-sequences/16S-7_sequences.qza \
  --m-barcodes-file 03_metadata/16S-7-full-metadata.txt \
  --m-barcodes-column 16S_barcode \
  --o-per-sample-sequences 04_artifacts-and-visualizations/demux-files/16S-7_demux.qza \
  --o-error-correction-details 04_artifacts-and-visualizations/demux-files/16S-7_demux_details.qza \
  --p-rev-comp-mapping-barcodes \
  --p-rev-comp-barcodes

qiime demux emp-single \
  --i-seqs 16S-sequences/16S-8_sequences.qza \
  --m-barcodes-file 03_metadata/16S-8-full-metadata.txt \
  --m-barcodes-column 16S_barcode \
  --o-per-sample-sequences 04_artifacts-and-visualizations/demux-files/16S-8_demux.qza \
  --o-error-correction-details 04_artifacts-and-visualizations/demux-files/16S-8_demux_details.qza \
  --p-rev-comp-mapping-barcodes \
  --p-rev-comp-barcodes
