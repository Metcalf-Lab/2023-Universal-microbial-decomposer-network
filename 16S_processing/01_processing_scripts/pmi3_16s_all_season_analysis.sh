## Qiime2 Import Spring Data
# Import qiita biom 58345
qiime tools import \
  --input-path qiita_files/11271_Spring/58345_reference-hit.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path qiita_files/11271_Spring/58345_reference-hit.qza

# Import the 58345 reference seqs
qiime tools import \
  --input-path qiita_files/11271_Spring/58345_reference-hit.seqs.fa \
  --output-path qiita_files/11271_Spring/58345_reference-hit.seqs.qza \
  --type 'FeatureData[Sequence]'
  
## Qiime2 Import Summer,Winter,Fall Data
# Import qiita biom 66414
qiime tools import \
  --input-path qiita_files/11489_Summer-Fall-Winter/66414_reference-hit.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path qiita_files/11489_Summer-Fall-Winter/66414_reference-hit.qza

# Import the 66414 reference seqs
qiime tools import \
  --input-path qiita_files/11489_Summer-Fall-Winter/66414_reference-hit.seqs.fa \
  --output-path qiita_files/11489_Summer-Fall-Winter/66414_reference-hit.seqs.qza \
  --type 'FeatureData[Sequence]'

# Import qiita biom 66370
qiime tools import \
  --input-path qiita_files/11489_Summer-Fall-Winter/66370_reference-hit.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path qiita_files/11489_Summer-Fall-Winter/66370_reference-hit.qza

# Import the 66370 reference seqs
qiime tools import \
  --input-path qiita_files/11489_Summer-Fall-Winter/66370_reference-hit.seqs.fa \
  --output-path qiita_files/11489_Summer-Fall-Winter/66370_reference-hit.seqs.qza \
  --type 'FeatureData[Sequence]'
  
# Import qiita biom 66419
qiime tools import \
  --input-path qiita_files/11489_Summer-Fall-Winter/66419_reference-hit.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path qiita_files/11489_Summer-Fall-Winter/66419_reference-hit.qza

# Import the 66419 reference seqs
qiime tools import \
  --input-path qiita_files/11489_Summer-Fall-Winter/66419_reference-hit.seqs.fa \
  --output-path qiita_files/11489_Summer-Fall-Winter/66419_reference-hit.seqs.qza \
  --type 'FeatureData[Sequence]'
  
# Import qiita biom 66359
qiime tools import \
  --input-path qiita_files/11489_Summer-Fall-Winter/66359_reference-hit.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path qiita_files/11489_Summer-Fall-Winter/66359_reference-hit.qza

# Import the 66359 reference seqs
qiime tools import \
  --input-path qiita_files/11489_Summer-Fall-Winter/66359_reference-hit.seqs.fa \
  --output-path qiita_files/11489_Summer-Fall-Winter/66359_reference-hit.seqs.qza \
  --type 'FeatureData[Sequence]'
  
# Import qiita biom 66404
qiime tools import \
  --input-path qiita_files/11489_Summer-Fall-Winter/66404_reference-hit.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path qiita_files/11489_Summer-Fall-Winter/66404_reference-hit.qza

# Import the 66404 reference seqs
qiime tools import \
  --input-path qiita_files/11489_Summer-Fall-Winter/66404_reference-hit.seqs.fa \
  --output-path qiita_files/11489_Summer-Fall-Winter/66404_reference-hit.seqs.qza \
  --type 'FeatureData[Sequence]'
  
# Import qiita biom 66406
qiime tools import \
  --input-path qiita_files/11489_Summer-Fall-Winter/66406_reference-hit.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path qiita_files/11489_Summer-Fall-Winter/66406_reference-hit.qza

# Import the 66406 reference seqs
qiime tools import \
  --input-path qiita_files/11489_Summer-Fall-Winter/66406_reference-hit.seqs.fa \
  --output-path qiita_files/11489_Summer-Fall-Winter/66406_reference-hit.seqs.qza \
  --type 'FeatureData[Sequence]'
  
# Import qiita biom 66410
qiime tools import \
  --input-path qiita_files/11489_Summer-Fall-Winter/66410_reference-hit.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path qiita_files/11489_Summer-Fall-Winter/66410_reference-hit.qza

# Import the 66410 reference seqs
qiime tools import \
  --input-path qiita_files/11489_Summer-Fall-Winter/66410_reference-hit.seqs.fa \
  --output-path qiita_files/11489_Summer-Fall-Winter/66410_reference-hit.seqs.qza \
  --type 'FeatureData[Sequence]'
  
# Import qiita biom 66412
qiime tools import \
  --input-path qiita_files/11489_Summer-Fall-Winter/66412_reference-hit.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path qiita_files/11489_Summer-Fall-Winter/66412_reference-hit.qza

# Import the 66412 reference seqs
qiime tools import \
  --input-path qiita_files/11489_Summer-Fall-Winter/66412_reference-hit.seqs.fa \
  --output-path qiita_files/11489_Summer-Fall-Winter/66412_reference-hit.seqs.qza \
  --type 'FeatureData[Sequence]'

## Merge tables into single table
qiime feature-table merge \
	--i-tables qiita_files/11271_Spring/58345_reference-hit.qza \
	--i-tables qiita_files/11489_Summer-Fall-Winter/66359_reference-hit.qza \
	--i-tables qiita_files/11489_Summer-Fall-Winter/66370_reference-hit.qza \
	--i-tables qiita_files/11489_Summer-Fall-Winter/66404_reference-hit.qza \
	--i-tables qiita_files/11489_Summer-Fall-Winter/66406_reference-hit.qza \
	--i-tables qiita_files/11489_Summer-Fall-Winter/66410_reference-hit.qza \
	--i-tables qiita_files/11489_Summer-Fall-Winter/66412_reference-hit.qza \
	--i-tables qiita_files/11489_Summer-Fall-Winter/66414_reference-hit.qza \
	--i-tables qiita_files/11489_Summer-Fall-Winter/66419_reference-hit.qza \
	--o-merged-table qiita_files/merged_table.qza

## Merge seqs into single file
qiime feature-table merge-seqs \
	--i-data qiita_files/11271_Spring/58345_reference-hit.seqs.qza \
	--i-data qiita_files/11489_Summer-Fall-Winter/66359_reference-hit.seqs.qza \
	--i-data qiita_files/11489_Summer-Fall-Winter/66370_reference-hit.seqs.qza \
	--i-data qiita_files/11489_Summer-Fall-Winter/66404_reference-hit.seqs.qza \
	--i-data qiita_files/11489_Summer-Fall-Winter/66406_reference-hit.seqs.qza \
	--i-data qiita_files/11489_Summer-Fall-Winter/66410_reference-hit.seqs.qza \
	--i-data qiita_files/11489_Summer-Fall-Winter/66412_reference-hit.seqs.qza \
	--i-data qiita_files/11489_Summer-Fall-Winter/66414_reference-hit.seqs.qza \
	--i-data qiita_files/11489_Summer-Fall-Winter/66419_reference-hit.seqs.qza \
	--o-merged-data qiita_files/merged_seqs.qza

## Feature table and data summaries
qiime feature-table summarize \
  --i-table qiita_files/merged_table.qza \
  --o-visualization qiita_files/merged_table.qzv \
  --m-sample-metadata-file 16s_all_seasons_map.txt

qiime feature-table tabulate-seqs \
  --i-data qiita_files/merged_seqs.qza \
  --o-visualization qiita_files/merged_seqs.qzv














