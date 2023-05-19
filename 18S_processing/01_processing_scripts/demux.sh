#!/bin/sh
#SBATCH --job-name=demux1_18S
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --partition=shas
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL

#qiime
source activate qiime2-2022.8

#1-6
qiime demux emp-single \
  --i-seqs 03_artifacts-and-visualizations/completed_18S_sequences/1-6_sequences.qza \
  --m-barcodes-file 02_metadata/plates_1-6_metadata.txt \
  --m-barcodes-column barcode \
  --o-per-sample-sequences 03_artifacts-and-visualizations/completed_18S_demux/1-6_18S_demux.qza \
  --o-error-correction-details 03_artifacts-and-visualizations/completed_18S_demux/1-6_18S_demux_details.qza \
  --p-rev-comp-mapping-barcodes
  
#7-12
qiime demux emp-single \
  --i-seqs 03_artifacts-and-visualizations/completed_18S_sequences/7-12_sequences.qza \
  --m-barcodes-file 02_metadata/plates_7-12_metadata.txt \
  --m-barcodes-column barcode \
  --o-per-sample-sequences 03_artifacts-and-visualizations/completed_18S_demux/7-12_18S_demux.qza \
  --o-error-correction-details 03_artifacts-and-visualizations/completed_18S_demux/7-12_18S_demux_details.qza \
  --p-rev-comp-mapping-barcodes

#13-17, 52
qiime demux emp-single \
  --i-seqs 03_artifacts-and-visualizations/completed_18S_sequences/13-17_52_18S_sequences.qza \
  --m-barcodes-file 02_metadata/plates_13-17_52_metadata.txt \
  --m-barcodes-column barcode \
  --o-per-sample-sequences 03_artifacts-and-visualizations/completed_18S_demux/13-17_52_18S_demux.qza \
  --o-error-correction-details 03_artifacts-and-visualizations/completed_18S_demux/13-17_52_18S_demux_details.qza \
  --p-rev-comp-mapping-barcodes

#18-21_24-25
qiime demux emp-single \
  --i-seqs 03_artifacts-and-visualizations/completed_18S_sequences/18-21_24-25_sequences.qza \
  --m-barcodes-file 02_metadata/plates_18-21_24-25_metadata.txt \
  --m-barcodes-column barcode \
  --o-per-sample-sequences 03_artifacts-and-visualizations/completed_18S_demux/18-21_24-25_18S_demux.qza \
  --o-error-correction-details 03_artifacts-and-visualizations/completed_18S_demux/18-21_24-25_18S_demux_details.qza \
  --p-rev-comp-mapping-barcodes

#26-28_30-31_49
qiime demux emp-single \
  --i-seqs 03_artifacts-and-visualizations/completed_18S_sequences/26-28_30-31_49_sequences.qza \
  --m-barcodes-file 02_metadata/plates_26-28_30-31_49_metadata.txt \
  --m-barcodes-column barcode \
  --o-per-sample-sequences 03_artifacts-and-visualizations/completed_18S_demux/26-28_30-31_49_18S_demux.qza \
  --o-error-correction-details 03_artifacts-and-visualizations/completed_18S_demux/26-28_30-31_49_18S_demux_details.qza \
  --p-rev-comp-mapping-barcodes

#29_48_51-53
qiime demux emp-single \
  --i-seqs 03_artifacts-and-visualizations/completed_18S_sequences/29_48_51-53_sequences.qza \
  --m-barcodes-file 02_metadata/plates_29_48_51-53_metadata.txt \
  --m-barcodes-column barcode \
  --o-per-sample-sequences 03_artifacts-and-visualizations/completed_18S_demux/29_48_51-53_18S_demux.qza \
  --o-error-correction-details 03_artifacts-and-visualizations/completed_18S_demux/29_48_51-53_18S_demux_details.qza \
  --p-rev-comp-mapping-barcodes

#32_37-41
qiime demux emp-single \
  --i-seqs 03_artifacts-and-visualizations/completed_18S_sequences/32_37-41_sequences.qza \
  --m-barcodes-file 02_metadata/plates_32_37-41_metadata.txt \
  --m-barcodes-column barcode \
  --o-per-sample-sequences 03_artifacts-and-visualizations/completed_18S_demux/32_37-41_18S_demux.qza \
  --o-error-correction-details 03_artifacts-and-visualizations/completed_18S_demux/32_37-41_18S_demux_details.qza \
  --p-rev-comp-mapping-barcodes

#42-47
qiime demux emp-single \
  --i-seqs 03_artifacts-and-visualizations/completed_18S_sequences/42-47_sequences.qza \
  --m-barcodes-file 02_metadata/plates_42-47_metadata.txt \
  --m-barcodes-column barcode \
  --o-per-sample-sequences 03_artifacts-and-visualizations/completed_18S_demux/42-47_18S_demux.qza \
  --o-error-correction-details 03_artifacts-and-visualizations/completed_18S_demux/42-47_18S_demux_details.qza \
  --p-rev-comp-mapping-barcodes

