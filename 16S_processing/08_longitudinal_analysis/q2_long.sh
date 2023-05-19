# Create rarefied species level table
qiime taxa collapse \
  --i-table all/rarefied_table.qza \
  --i-taxonomy ../07_taxonomy/pmi3-16S-silva-taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table all/species_rarefied_table.qza

qiime taxa collapse \
  --i-table all/PMI-16S-nochlomito-filtered-table.qza \
  --i-taxonomy ../07_taxonomy/pmi3-16S-silva-taxonomy.qza \
  --p-level 7 \
  --o-collapsed-table all/species_PMI-16S-nochlomito-filtered-table.qza


# Split out soil and skin
qiime feature-table filter-samples \
  --i-table all/species_rarefied_table.qza \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --p-where "[sample_site]='skin.face' OR [sample_site]='skin.hip'" \
  --o-filtered-table all/skin-rarefied-table.qza

qiime feature-table filter-samples \
  --i-table all/species_rarefied_table.qza \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --p-where "[sample_site]='soil.face' OR [sample_site]='soil.hip'" \
  --o-filtered-table all/soil-rarefied-table.qza

qiime feature-table filter-samples \
  --i-table all/species_PMI-16S-nochlomito-filtered-table.qza \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --p-where "[sample_site]='skin.face' OR [sample_site]='skin.hip'" \
  --o-filtered-table all/skin-table.qza

qiime feature-table filter-samples \
  --i-table all/species_PMI-16S-nochlomito-filtered-table.qza \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --p-where "[sample_site]='soil.face' OR [sample_site]='soil.hip'" \
  --o-filtered-table all/soil-table.qza

# Add 0c
qiime longitudinal feature-volatility \
  --i-table all/soil-rarefied-table.qza \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --p-state-column add_0c \
  --p-individual-id-column host_subject_id \
  --p-estimator 'RandomForestRegressor' \
  --p-parameter-tuning \
  --p-importance-threshold 'q1' \
  --p-feature-count 100 \
  --p-n-estimators 1000 \
  --p-random-state 999 \
  --p-n-jobs -1 \
  --output-dir all/add_0c-soil-feature-volatility-species 


# Add 0c
qiime longitudinal feature-volatility \
  --i-table all/skin-rarefied-table.qza \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --p-state-column add_0c \
  --p-individual-id-column host_subject_id \
  --p-estimator 'RandomForestRegressor' \
  --p-parameter-tuning \
  --p-importance-threshold 'q1' \
  --p-feature-count 100 \
  --p-n-estimators 1000 \
  --p-random-state 999 \
  --p-n-jobs -1 \
  --output-dir all/add_0c-skin-feature-volatility-species 

### Create LME with rel freq table of time series and metric of important features 
## https://github.com/caporaso-lab/longitudinal-notebooks/blob/master/notebooks/analysis.ipynb
# gOTU chosen based on their importance in our ML models

# Filter time series to 0-50 ADD at gOTU level
qiime feature-table filter-samples \
  --i-table all/soil-rarefied-table.qza \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --p-where "[decomp_stage]=='early'" \
  --o-filtered-table all/soil-rarefied-early-table.qza

# Filter time series to 0-50 ADD at gOTU level
qiime feature-table filter-samples \
  --i-table all/soil-rarefied-table.qza \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --p-where "[decomp_stage]=='active'" \
  --o-filtered-table all/soil-rarefied-active-table.qza

# Filter time series to 0-50 ADD at gOTU level
qiime feature-table filter-samples \
  --i-table all/soil-rarefied-table.qza \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --p-where "[decomp_stage]=='advanced'" \
  --o-filtered-table all/soil-rarefied-advanced-table.qza


# Convert to relative abundances
qiime feature-table relative-frequency \
	--i-table all/soil-rarefied-table.qza \
	--o-relative-frequency-table all/soil-rarefied-RF-table.qza

qiime feature-table relative-frequency \
	--i-table all/soil-rarefied-early-table.qza \
	--o-relative-frequency-table all/soil-rarefied-early-RF-table.qza

qiime feature-table relative-frequency \
	--i-table all/soil-rarefied-active-table.qza \
	--o-relative-frequency-table all/soil-rarefied-active-RF-table.qza

qiime feature-table relative-frequency \
	--i-table all/soil-rarefied-advanced-table.qza \
	--o-relative-frequency-table all/soil-rarefied-advanced-RF-table.qza


qiime feature-table relative-frequency \
	--i-table all/skin-rarefied-table.qza \
	--o-relative-frequency-table all/skin-rarefied-RF-table.qza

qiime feature-table relative-frequency \
	--i-table all/soil-table.qza \
	--o-relative-frequency-table all/soil-RF-table.qza

qiime feature-table relative-frequency \
	--i-table all/skin-table.qza \
	--o-relative-frequency-table all/skin-RF-table.qza

### LME models
mkdir all/lme

## LME Helcococcus seattlensis
# Full ADD
qiime longitudinal linear-mixed-effects \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --i-table all/soil-rarefied-RF-table.qza \
  --p-metric 'D_0__Bacteria;D_1__Firmicutes;D_2__Clostridia;D_3__Clostridiales;D_4__Family XI;D_5__Helcococcus;D_6__Helcococcus seattlensis' \
  --p-group-columns facility \
  --p-state-column add_0c \
  --p-individual-id-column host_subject_id \
  --o-visualization all/lme/Helcococcus_seattlensis-lme-soil-full.qzv

# early ADD
qiime longitudinal linear-mixed-effects \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --i-table all/soil-rarefied-early-RF-table.qza \
  --p-metric 'D_0__Bacteria;D_1__Firmicutes;D_2__Clostridia;D_3__Clostridiales;D_4__Family XI;D_5__Helcococcus;D_6__Helcococcus seattlensis' \
  --p-group-columns facility \
  --p-state-column add_0c \
  --p-individual-id-column host_subject_id \
  --o-visualization all/lme/Helcococcus_seattlensis-lme-soil-early.qzv

# active ADD
qiime longitudinal linear-mixed-effects \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --i-table all/soil-rarefied-active-RF-table.qza \
  --p-metric 'D_0__Bacteria;D_1__Firmicutes;D_2__Clostridia;D_3__Clostridiales;D_4__Family XI;D_5__Helcococcus;D_6__Helcococcus seattlensis' \
  --p-group-columns facility \
  --p-state-column add_0c \
  --p-individual-id-column host_subject_id \
  --o-visualization all/lme/Helcococcus_seattlensis-lme-soil-active.qzv

# advanced ADD
qiime longitudinal linear-mixed-effects \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --i-table all/soil-rarefied-advanced-RF-table.qza \
  --p-metric 'D_0__Bacteria;D_1__Firmicutes;D_2__Clostridia;D_3__Clostridiales;D_4__Family XI;D_5__Helcococcus;D_6__Helcococcus seattlensis' \
  --p-group-columns facility \
  --p-state-column add_0c \
  --p-individual-id-column host_subject_id \
  --o-visualization all/lme/Helcococcus_seattlensis-lme-soil-advanced.qzv





# Soil hip + control
qiime diversity filter-distance-matrix \
  --i-distance-matrix /Users/Zach/Dropbox/PMI_3_analyses/multi-omics_data/amplicon/16S/06_core-metrics/filtered-core-metrics-5000/generalized-unifrac05-distance-matrix.qza \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt  \
  --p-where "([soil_control]='y'AND [facility]='SHSU') OR ([sample_site]='soil.hip' AND [facility]='SHSU')" \
  --o-filtered-distance-matrix soilhip/soilhip-con-unifrac05-SHSU-mat.qza

qiime diversity filter-distance-matrix \
  --i-distance-matrix /Users/Zach/Dropbox/PMI_3_analyses/multi-omics_data/amplicon/16S/06_core-metrics/filtered-core-metrics-5000/generalized-unifrac05-distance-matrix.qza \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt  \
  --p-where "([soil_control]='y'AND [facility]='UTK') OR ([sample_site]='soil.hip' AND [facility]='UTK')" \
  --o-filtered-distance-matrix soilhip/soilhip-con-unifrac05-UTK-mat.qza

qiime diversity filter-distance-matrix \
  --i-distance-matrix /Users/Zach/Dropbox/PMI_3_analyses/multi-omics_data/amplicon/16S/06_core-metrics/filtered-core-metrics-5000/generalized-unifrac05-distance-matrix.qza \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt  \
  --p-where "([soil_control]='y'AND [facility]='CMU') OR ([sample_site]='soil.hip' AND [facility]='CMU')" \
  --o-filtered-distance-matrix soilhip/soilhip-con-unifrac05-CMU-mat.qza
 
qiime diversity pcoa \
	--i-distance-matrix soilhip/soilhip-con-unifrac05-SHSU-mat.qza \
	--o-pcoa soilhip/soilhip-con-unifrac05-SHSU-pcoa.qza

qiime diversity pcoa \
	--i-distance-matrix soilhip/soilhip-con-unifrac05-UTK-mat.qza \
	--o-pcoa soilhip/soilhip-con-unifrac05-UTK-pcoa.qza

qiime diversity pcoa \
	--i-distance-matrix soilhip/soilhip-con-unifrac05-CMU-mat.qza \
	--o-pcoa soilhip/soilhip-con-unifrac05-CMU-pcoa.qza

qiime emperor plot \
	--i-pcoa soilhip/soilhip-con-unifrac05-CMU-pcoa.qza \
	--m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
	--o-visualization soilhip/soilhip-con-unifrac05-CMU-emp.qzv

qiime longitudinal volatility \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --m-metadata-file soilhip/soilhip-con-unifrac05-SHSU-pcoa.qza \
  --p-state-column add_0c \
  --p-individual-id-column host_subject_id \
  --p-default-group-column 'soil_control' \
  --p-default-metric 'Axis 1' \
  --o-visualization soilhip/pc_SHSU_vol.qzv

qiime longitudinal volatility \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --m-metadata-file soilhip/soilhip-con-unifrac05-UTK-pcoa.qza \
  --p-state-column add_0c \
  --p-individual-id-column host_subject_id \
  --p-default-group-column 'soil_control' \
  --p-default-metric 'Axis 1' \
  --o-visualization soilhip/pc_UTK_vol.qzv

qiime longitudinal volatility \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --m-metadata-file soilhip/soilhip-con-unifrac05-CMU-pcoa.qza \
  --p-state-column add_0c \
  --p-individual-id-column host_subject_id \
  --p-default-group-column 'soil_control' \
  --p-default-metric 'Axis 1' \
  --o-visualization soilhip/pc_CMU_vol.qzv


### LME models
mkdir soilhip/lme

qiime longitudinal linear-mixed-effects \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --m-metadata-file soilhip/soilhip-con-unifrac05-SHSU-pcoa.qza \
  --p-metric 'Axis 1' \
  --p-group-columns sample_site \
  --p-state-column add_0c \
  --p-individual-id-column host_subject_id \
  --o-visualization soilhip/lme/SHSU-unifrac05-axis1.qzv

qiime longitudinal linear-mixed-effects \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --m-metadata-file soilhip/soilhip-con-unifrac05-UTK-pcoa.qza \
  --p-metric 'Axis 1' \
  --p-group-columns sample_site \
  --p-state-column add_0c \
  --p-individual-id-column host_subject_id \
  --o-visualization soilhip/lme/UTK-unifrac05-axis1.qzv

qiime longitudinal linear-mixed-effects \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --m-metadata-file soilhip/soilhip-con-unifrac05-CMU-pcoa.qza \
  --p-metric 'Axis 1' \
  --p-group-columns sample_site \
  --p-state-column add_0c \
  --p-individual-id-column host_subject_id \
  --o-visualization soilhip/lme/CMU-unifrac05-axis1.qzv

qiime longitudinal linear-mixed-effects \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --m-metadata-file soilhip/soilhip-con-unifrac05-SHSU-pcoa.qza \
  --p-metric 'Axis 2' \
  --p-group-columns sample_site \
  --p-state-column add_0c \
  --p-individual-id-column host_subject_id \
  --o-visualization soilhip/lme/SHSU-unifrac05-axis2.qzv
			
# richness lme of decomp soil per facility
# SHSU
qiime feature-table filter-samples \
  --i-table soilhip/soilhip-con-rarefied-table.qza \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --p-where "([soil_control]='y'AND [facility]='SHSU') OR ([sample_site]='soil.hip' AND [facility]='SHSU')" \
  --o-filtered-table soilhip/soilhip-SHSU-rarefied-table.qza
  
qiime diversity alpha \
	--i-table	soilhip/soilhip-SHSU-rarefied-table.qza \
	--p-metric 'observed_otus' \
	--o-alpha-diversity soilhip/soilhip-SHSU-richness.qza
	
qiime longitudinal linear-mixed-effects \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --m-metadata-file soilhip/soilhip-SHSU-richness.qza \
  --p-metric 'observed_otus' \
  --p-group-columns sample_site \
  --p-state-column add_0c \
  --p-individual-id-column host_subject_id \
  --o-visualization soilhip/lme/richness-SHSU-lme.qzv
  
# UTK
qiime feature-table filter-samples \
  --i-table soilhip/soilhip-con-rarefied-table.qza \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --p-where "([soil_control]='y'AND [facility]='UTK') OR ([sample_site]='soil.hip' AND [facility]='UTK')" \
  --o-filtered-table soilhip/soilhip-UTK-rarefied-table.qza
  
qiime diversity alpha \
	--i-table	soilhip/soilhip-UTK-rarefied-table.qza \
	--p-metric 'observed_otus' \
	--o-alpha-diversity soilhip/soilhip-UTK-richness.qza
	
qiime longitudinal linear-mixed-effects \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --m-metadata-file soilhip/soilhip-UTK-richness.qza \
  --p-metric 'observed_otus' \
  --p-group-columns sample_site \
  --p-state-column add_0c \
  --p-individual-id-column host_subject_id \
  --o-visualization soilhip/lme/richness-UTK-lme.qzv
  
# CMU
qiime feature-table filter-samples \
  --i-table soilhip/soilhip-con-rarefied-table.qza \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --p-where "([soil_control]='y'AND [facility]='CMU') OR ([sample_site]='soil.hip' AND [facility]='CMU')" \
  --o-filtered-table soilhip/soilhip-CMU-rarefied-table.qza
  
qiime diversity alpha \
	--i-table	soilhip/soilhip-CMU-rarefied-table.qza \
	--p-metric 'observed_otus' \
	--o-alpha-diversity soilhip/soilhip-CMU-richness.qza
	
qiime longitudinal linear-mixed-effects \
  --m-metadata-file ../03_metadata/combined-metadata-simple-nov2020.txt \
  --m-metadata-file soilhip/soilhip-CMU-richness.qza \
  --p-metric 'observed_otus' \
  --p-group-columns sample_site \
  --p-state-column add_0c \
  --p-individual-id-column host_subject_id \
  --o-visualization soilhip/lme/richness-CMU-lme.qzv