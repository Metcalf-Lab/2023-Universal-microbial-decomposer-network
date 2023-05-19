#!/bin/sh
#SBATCH --job-name=16S-demux-2
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
  --i-seqs 16S-sequences/16S-3_sequences.qza \
  --m-barcodes-file 03_metadata/16S-3-full-metadata.txt \
  --m-barcodes-column 16S_barcode \
  --o-per-sample-sequences 04_artifacts-and-visualizations/demux-files/16S-3_demux.qza \
  --o-error-correction-details 04_artifacts-and-visualizations/demux-files/16S-3_demux_details.qza \
  --p-rev-comp-mapping-barcodes \
  --p-rev-comp-barcodes

qiime demux emp-single \
  --i-seqs 16S-sequences/16S-4_sequences.qza \
  --m-barcodes-file 03_metadata/16S-4-full-metadata.txt \
  --m-barcodes-column 16S_barcode \
  --o-per-sample-sequences 04_artifacts-and-visualizations/demux-files/16S-4_demux.qza \
  --o-error-correction-details 04_artifacts-and-visualizations/demux-files/16S-4_demux_details.qza \
  --p-rev-comp-mapping-barcodes \
  --p-rev-comp-barcodes
