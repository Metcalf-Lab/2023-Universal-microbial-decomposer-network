########################################
# Qiime2 analysis of metabolomics data
########################################
conda activate qiime2-2022.2
cd /Users/zacharyburcham/Dropbox/PMI_3_analyses/multi-omics_data/metabolomic/pieters_fbmn_run/qiime2-2022.2
mkdir   metadata feature_tables


# input normalized, scaled table
biom convert -i ../quantification_table/norm-scaled-feature-table.txt -o ../quantification_table/norm-scaled-feature-table.biom --table-type="OTU table" --to-hdf5

qiime tools import \
  --input-path ../quantification_table/norm-scaled-feature-table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path feature_tables/qiime2_table_normalized.qza


########################################
### Create PCOAs to check where controls fall
mkdir control_check
## Jaccard
qiime diversity beta \
	--i-table feature_tables/qiime2_table_normalized.qza \
	--p-metric 'jaccard' \
	--o-distance-matrix control_check/jd_matrix.qza
 
# pcoa
qiime diversity pcoa \
	--i-distance-matrix control_check/jd_matrix.qza \
	--o-pcoa control_check/jd_pcoa.qza

# emperor	
qiime emperor plot \
	--i-pcoa control_check/jd_pcoa.qza \
	--m-metadata-file metadata/pmi3_metab_meta_q2.txt \
	--o-visualization control_check/jd_pcoa_emp.qzv

## Bray-Curtis
qiime diversity beta \
	--i-table feature_tables/qiime2_table_normalized.qza \
	--p-metric 'braycurtis' \
	--o-distance-matrix control_check/bc_matrix.qza
 
# pcoa
qiime diversity pcoa \
	--i-distance-matrix control_check/bc_matrix.qza \
	--o-pcoa control_check/bc_pcoa.qza

# emperor	
qiime emperor plot \
	--i-pcoa control_check/bc_pcoa.qza \
	--m-metadata-file metadata/pmi3_metab_meta_q2.txt \
	--o-visualization control_check/bc_pcoa_emp.qzv

## Skin and hip are clearly seperated in PCOAs
######################################

# Filter all controls and samples falling near controls in PCOA
qiime feature-table filter-samples \
  --i-table feature_tables/qiime2_table_normalized.qza \
  --m-metadata-file metadata/pmi3_metab_meta_q2.txt \
  --p-where "[proccessing_control]='y' OR [pcoa_removals]='y'" \
  --p-exclude-ids \
  --o-filtered-table feature_tables/filtered_table.qza

# Summarize table
qiime feature-table summarize \
  --i-table feature_tables/filtered_table.qza \
  --o-visualization feature_tables/filtered_table.qzv \
  --m-sample-metadata-file metadata/pmi3_metab_meta_q2.txt

# Filter features in less than 30 samples
qiime feature-table filter-features \
  --i-table feature_tables/filtered_table.qza \
  --p-min-samples 30 \
  --o-filtered-table feature_tables/final_filtered_table.qza

# Summarize table
qiime feature-table summarize \
  --i-table feature_tables/final_filtered_table.qza \
  --o-visualization feature_tables/final_filtered_table.qzv \
  --m-sample-metadata-file metadata/pmi3_metab_meta_q2.txt

# Export biom table
qiime tools export \
  --input-path feature_tables/final_filtered_table.qza \
  --output-path feature_tables/final_filtered_table
  
mv feature_tables/final_filtered_table/feature-table.biom feature_tables/final_filtered_table.biom
rmdir feature_tables/final_filtered_table

# Richness
qiime diversity alpha \
	--i-table feature_tables/final_filtered_table.qza \
	--p-metric 'observed_features' \
	--o-alpha-diversity control_check/observed_features-vector.qza
	
### Final full PCoAs
## Jaccard
qiime diversity beta \
	--i-table feature_tables/final_filtered_table.qza \
	--p-metric 'jaccard' \
	--o-distance-matrix control_check/jd_matrix_filtered.qza
 
# pcoa
qiime diversity pcoa \
	--i-distance-matrix control_check/jd_matrix_filtered.qza \
	--o-pcoa control_check/jd_pcoa_filtered.qza

# emperor	
qiime emperor plot \
	--i-pcoa control_check/jd_pcoa_filtered.qza \
	--m-metadata-file metadata/pmi3_metab_meta_q2.txt \
	--o-visualization control_check/jd_pcoa_emp_filtered.qzv

# ADD on axis
qiime emperor plot \
  --i-pcoa control_check/jd_pcoa_filtered.qza \
  --m-metadata-file metadata/pmi3_metab_meta_q2.txt \
  --p-custom-axes add_0c \
  --o-visualization control_check/jd-emp-add0c.qzv

# beta group significance between sample types
qiime diversity beta-group-significance \
	--i-distance-matrix control_check/jd_matrix_filtered.qza \
	--m-metadata-file metadata/pmi3_metab_meta_q2.txt \
	--m-metadata-column 'sample_type' \
	--p-pairwise \
	--o-visualization control_check/jd_sample_type_permanova.qzv

## Bray Curtis
qiime diversity beta \
	--i-table feature_tables/final_filtered_table.qza \
	--p-metric 'braycurtis' \
	--o-distance-matrix control_check/bc_matrix_filtered.qza
 
# pcoa
qiime diversity pcoa \
	--i-distance-matrix control_check/bc_matrix_filtered.qza \
	--o-pcoa control_check/bc_pcoa_filtered.qza

# emperor	
qiime emperor plot \
	--i-pcoa control_check/bc_pcoa_filtered.qza \
	--m-metadata-file metadata/pmi3_metab_meta_q2.txt \
	--o-visualization control_check/bc_pcoa_emp_filtered.qzv

# ADD on axis
qiime emperor plot \
  --i-pcoa control_check/bc_pcoa_filtered.qza \
  --m-metadata-file metadata/pmi3_metab_meta_q2.txt \
  --p-custom-axes add_0c \
  --o-visualization control_check/bc-emp-add0c.qzv

# beta group significance between sample types
qiime diversity beta-group-significance \
	--i-distance-matrix control_check/bc_matrix_filtered.qza \
	--m-metadata-file metadata/pmi3_metab_meta_q2.txt \
	--m-metadata-column 'sample_type' \
	--p-pairwise \
	--o-visualization control_check/bc_sample_type_permanova.qzv


### FILTER TABLE TO METABS WITH SIRIUS FORMULA
qiime tools import \
  --input-path filtering/for-filtering-sirius.txt \
  --output-path filtering/for-filtering-sirius.qza \
  --type 'FeatureData[filtering]'
  
qiime metadata tabulate \
  --m-input-file filtering/for-filtering-sirius.qza \
  --o-visualization filtering/for-filtering-sirius.qzv
  
# less strict
qiime taxa filter-table \
  --i-table feature_tables/final_filtered_table.qza \
  --i-filtering filtering/for-filtering-sirius.qza \
  --p-exclude 'none' \
  --p-mode 'exact' \
  --o-filtered-table feature_tables/library_less_strict_filtered_table.qza

qiime feature-table summarize \
  --i-table feature_tables/library_less_strict_filtered_table.qza \
  --o-visualization feature_tables/library_less_strict_filtered_table.qzv \
  --m-sample-metadata-file metadata/pmi3_metab_meta_q2.txt
  
# remove soil con
qiime feature-table filter-samples \
  --i-table feature_tables/library_less_strict_filtered_table.qza \
  --m-metadata-file metadata/pmi3_metab_meta_q2.txt \
  --p-where "[soil_control]='n'" \
  --o-filtered-table feature_tables/library_less_strict_filtered_nocon_table.qza
  
qiime feature-table summarize \
  --i-table feature_tables/library_less_strict_filtered_nocon_table.qza \
  --o-visualization feature_tables/library_less_strict_filtered_nocon_table.qzv \
  --m-sample-metadata-file metadata/pmi3_metab_meta_q2.txt

## Jaccard
qiime diversity beta \
	--i-table feature_tables/library_less_strict_filtered_nocon_table.qza \
	--p-metric 'jaccard' \
	--o-distance-matrix control_check/library_less_strict_filtered_nocon_jd_matrix.qza
 
# pcoa
qiime diversity pcoa \
	--i-distance-matrix control_check/library_less_strict_filtered_nocon_jd_matrix.qza \
	--o-pcoa control_check/library_less_strict_filtered_nocon_jd_pcoa.qza

# emperor	
qiime emperor plot \
	--i-pcoa control_check/library_less_strict_filtered_nocon_jd_pcoa.qza \
	--m-metadata-file metadata/pmi3_metab_meta_q2.txt \
	--o-visualization control_check/library_less_strict_filtered_nocon_jd_emp.qzv

## Bray Curtis
qiime diversity beta \
	--i-table feature_tables/library_less_strict_filtered_nocon_table.qza \
	--p-metric 'braycurtis' \
	--o-distance-matrix control_check/library_less_strict_filtered_nocon_bc_matrix.qza
 
# pcoa
qiime diversity pcoa \
	--i-distance-matrix control_check/library_less_strict_filtered_nocon_bc_matrix.qza \
	--o-pcoa control_check/library_less_strict_filtered_nocon_bc_pcoa.qza

# emperor	
qiime emperor plot \
	--i-pcoa control_check/library_less_strict_filtered_nocon_bc_pcoa.qza \
	--m-metadata-file metadata/pmi3_metab_meta_q2.txt \
	--o-visualization control_check/library_less_strict_filtered_nocon_bc_emp.qzv

#### Soil metabolites
# keep soil
qiime feature-table filter-samples \
  --i-table feature_tables/library_less_strict_filtered_nocon_table.qza \
  --m-metadata-file metadata/pmi3_metab_meta_q2.txt \
  --p-where "[sample_type]='soil'" \
  --o-filtered-table feature_tables/library_less_strict_filtered_nocon_table_soil.qza
  
qiime feature-table summarize \
  --i-table feature_tables/library_less_strict_filtered_nocon_table_soil.qza \
  --o-visualization feature_tables/library_less_strict_filtered_nocon_table_soil.qzv \
  --m-sample-metadata-file metadata/pmi3_metab_meta_q2.txt

## Jaccard
qiime diversity beta \
	--i-table feature_tables/library_less_strict_filtered_nocon_table_soil.qza \
	--p-metric 'jaccard' \
	--o-distance-matrix diversity/library_less_strict_filtered_nocon_soil_jd_mat.qza
 
# pcoa
qiime diversity pcoa \
	--i-distance-matrix diversity/library_less_strict_filtered_nocon_soil_jd_mat.qza \
	--o-pcoa diversity/library_less_strict_filtered_nocon_soil_jd_pcoa.qza

# emperor	
qiime emperor plot \
	--i-pcoa diversity/library_less_strict_filtered_nocon_soil_jd_pcoa.qza \
	--m-metadata-file metadata/pmi3_metab_meta_q2.txt \
	--o-visualization diversity/library_less_strict_filtered_nocon_soil_jd_emp.qzv

## Bray Curtis
qiime diversity beta \
	--i-table feature_tables/library_less_strict_filtered_nocon_table_soil.qza \
	--p-metric 'braycurtis' \
	--o-distance-matrix diversity/library_less_strict_filtered_nocon_soil_bc_mat.qza
 
# pcoa
qiime diversity pcoa \
	--i-distance-matrix diversity/library_less_strict_filtered_nocon_soil_bc_mat.qza \
	--o-pcoa diversity/library_less_strict_filtered_nocon_soil_bc_pcoa.qza

# emperor	
qiime emperor plot \
	--i-pcoa diversity/library_less_strict_filtered_nocon_soil_bc_pcoa.qza \
	--m-metadata-file metadata/pmi3_metab_meta_q2.txt \
	--o-visualization diversity/library_less_strict_filtered_nocon_soil_bc_emp.qzv


#### Skin metabolites
# keep skin
qiime feature-table filter-samples \
  --i-table feature_tables/library_less_strict_filtered_nocon_table.qza \
  --m-metadata-file metadata/pmi3_metab_meta_q2.txt \
  --p-where "[sample_type]='skin'" \
  --o-filtered-table feature_tables/library_less_strict_filtered_nocon_table_skin.qza
  
qiime feature-table summarize \
  --i-table feature_tables/library_less_strict_filtered_nocon_table_skin.qza \
  --o-visualization feature_tables/library_less_strict_filtered_nocon_table_skin.qzv \
  --m-sample-metadata-file metadata/pmi3_metab_meta_q2.txt

## Jaccard
qiime diversity beta \
	--i-table feature_tables/library_less_strict_filtered_nocon_table_skin.qza \
	--p-metric 'jaccard' \
	--o-distance-matrix diversity/library_less_strict_filtered_nocon_skin_jd_mat.qza
 
# pcoa
qiime diversity pcoa \
	--i-distance-matrix diversity/library_less_strict_filtered_nocon_skin_jd_mat.qza \
	--o-pcoa diversity/library_less_strict_filtered_nocon_skin_jd_pcoa.qza

# emperor	
qiime emperor plot \
	--i-pcoa diversity/library_less_strict_filtered_nocon_skin_jd_pcoa.qza \
	--m-metadata-file metadata/pmi3_metab_meta_q2.txt \
	--o-visualization diversity/library_less_strict_filtered_nocon_skin_jd_emp.qzv

## Bray Curtis
qiime diversity beta \
	--i-table feature_tables/library_less_strict_filtered_nocon_table_skin.qza \
	--p-metric 'braycurtis' \
	--o-distance-matrix diversity/library_less_strict_filtered_nocon_skin_bc_mat.qza
 
# pcoa
qiime diversity pcoa \
	--i-distance-matrix diversity/library_less_strict_filtered_nocon_skin_bc_mat.qza \
	--o-pcoa diversity/library_less_strict_filtered_nocon_skin_bc_pcoa.qza

# emperor	
qiime emperor plot \
	--i-pcoa diversity/library_less_strict_filtered_nocon_skin_bc_pcoa.qza \
	--m-metadata-file metadata/pmi3_metab_meta_q2.txt \
	--o-visualization diversity/library_less_strict_filtered_nocon_skin_bc_emp.qzv



