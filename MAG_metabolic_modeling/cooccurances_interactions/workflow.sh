
### Cooccurances and SMETANA
qiime feature-table filter-features \
	--i-table feature_tables/MAG_table_final.qza \
	--m-metadata-file models.txt \
	--p-exclude-ids \
	--o-filtered-table feature_tables/MAG_table_final_models.qza
	
qiime feature-table relative-frequency \
	--i-table feature_tables/MAG_table_final_models.qza \
	--o-relative-frequency-table feature_tables/MAG_table_final_models-relfreq.qza

## Separate tables by decomp stage and calculate co-occurances
# FIRS early
qiime feature-table filter-samples \
  --i-table feature_tables/MAG_table_final_models-relfreq.qza \
  --m-metadata-file metadata/shotgun_metadata_q2.txt \
  --p-where "[add_facility]='earlyFIRS'" \
  --o-filtered-table feature_tables/FIRS-early-relfreq-table.qza

qiime tools export \
  --input-path feature_tables/FIRS-early-relfreq-table.qza \
  --output-path feature_tables/FIRS-early-relfreq-table
  
mv feature_tables/FIRS-early-relfreq-table/feature-table.biom feature_tables/FIRS-early-relfreq-table.biom
rmdir feature_tables/FIRS-early-relfreq-table

biom convert -i feature_tables/FIRS-early-relfreq-table.biom -o feature_tables/FIRS-early-relfreq-table.tsv --to-tsv

sed -i .bak '1d' feature_tables/FIRS-early-relfreq-table.tsv
sed -i .bak 's/#//' feature_tables/FIRS-early-relfreq-table.tsv
rm feature_tables/*.bak
mkdir cooccurances FIRS_cooccurances/ FIRS_cooccurances//early

hiorco -k 20 --cutoff 0.001 -o FIRS_cooccurances//early feature_tables/FIRS-early-relfreq-table.tsv

# FIRS active
qiime feature-table filter-samples \
  --i-table feature_tables/MAG_table_final_models-relfreq.qza \
  --m-metadata-file metadata/shotgun_metadata_q2.txt \
  --p-where "[add_facility]='activeFIRS'" \
  --o-filtered-table feature_tables/FIRS-active-relfreq-table.qza

qiime tools export \
  --input-path feature_tables/FIRS-active-relfreq-table.qza \
  --output-path feature_tables/FIRS-active-relfreq-table
  
mv feature_tables/FIRS-active-relfreq-table/feature-table.biom feature_tables/FIRS-active-relfreq-table.biom
rmdir feature_tables/FIRS-active-relfreq-table

biom convert -i feature_tables/FIRS-active-relfreq-table.biom -o feature_tables/FIRS-active-relfreq-table.tsv --to-tsv

sed -i .bak '1d' feature_tables/FIRS-active-relfreq-table.tsv
sed -i .bak 's/#//' feature_tables/FIRS-active-relfreq-table.tsv
rm feature_tables/*.bak
mkdir FIRS_cooccurances//active

hiorco -k 20 --cutoff 0.001 -o FIRS_cooccurances//active feature_tables/FIRS-active-relfreq-table.tsv

# FIRS advanced
qiime feature-table filter-samples \
  --i-table feature_tables/MAG_table_final_models-relfreq.qza \
  --m-metadata-file metadata/shotgun_metadata_q2.txt \
  --p-where "[add_facility]='advancedFIRS'" \
  --o-filtered-table feature_tables/FIRS-advanced-relfreq-table.qza

qiime tools export \
  --input-path feature_tables/FIRS-advanced-relfreq-table.qza \
  --output-path feature_tables/FIRS-advanced-relfreq-table
  
mv feature_tables/FIRS-advanced-relfreq-table/feature-table.biom feature_tables/FIRS-advanced-relfreq-table.biom
rmdir feature_tables/FIRS-advanced-relfreq-table

biom convert -i feature_tables/FIRS-advanced-relfreq-table.biom -o feature_tables/FIRS-advanced-relfreq-table.tsv --to-tsv

sed -i .bak '1d' feature_tables/FIRS-advanced-relfreq-table.tsv
sed -i .bak 's/#//' feature_tables/FIRS-advanced-relfreq-table.tsv
rm feature_tables/*.bak
mkdir FIRS_cooccurances//advanced

hiorco -k 20 --cutoff 0.001 -o FIRS_cooccurances//advanced feature_tables/FIRS-advanced-relfreq-table.tsv

## Using co-occurrences for SMETANA - global
# FIRS community size 20 - aerobic
smetana /Users/Zach/Dropbox/PMI_3_analyses/multi-omics_data/shotgun/wrighton_results/Models/*.xml -c FIRS_cooccurances//early/size_20_part_1.tsv -g -v --flavor ucsd --molweight --aerobic --exclude cooccurances/inorganic.txt -o FIRS_cooccurances/20_early_aerobic
smetana /Users/Zach/Dropbox/PMI_3_analyses/multi-omics_data/shotgun/wrighton_results/Models/*.xml -c FIRS_cooccurances//active/size_20_part_1.tsv -g -v --flavor ucsd --molweight --aerobic --exclude cooccurances/inorganic.txt -o FIRS_cooccurances/20_active_aerobic

# FIRS community size 20
smetana Models/*.xml -c FIRS_cooccurances//early/size_20_part_1.tsv -g -v --flavor ucsd --molweight --exclude cooccurances/inorganic.txt -o FIRS_cooccurances/20_early
smetana Models/*.xml -c FIRS_cooccurances//active/size_20_part_1.tsv -g -v --flavor ucsd --molweight --exclude cooccurances/inorganic.txt -o FIRS_cooccurances/20_active

# FIRS active detailed - aerobic
smetana Models/*.xml -c FIRS_cooccurances//active/size_20_part_1.tsv -d -v --flavor ucsd --molweight --aerobic --exclude cooccurances/inorganic.txt -o FIRS_cooccurances/20_active_aerobic

# FIRS active detailed 
smetana Models/*.xml -c FIRS_cooccurances//active/size_20_part_1.tsv -d -v --flavor ucsd --molweight --exclude cooccurances/inorganic.txt -o FIRS_cooccurances/20_active


## Separate tables by decomp stage and calculate co-occurances
# STAFS early
qiime feature-table filter-samples \
  --i-table feature_tables/MAG_table_final_models-relfreq.qza \
  --m-metadata-file metadata/shotgun_metadata_q2.txt \
  --p-where "[add_facility]='earlySTAFS'" \
  --o-filtered-table feature_tables/STAFS-early-relfreq-table.qza

qiime tools export \
  --input-path feature_tables/STAFS-early-relfreq-table.qza \
  --output-path feature_tables/STAFS-early-relfreq-table
  
mv feature_tables/STAFS-early-relfreq-table/feature-table.biom feature_tables/STAFS-early-relfreq-table.biom
rmdir feature_tables/STAFS-early-relfreq-table

biom convert -i feature_tables/STAFS-early-relfreq-table.biom -o feature_tables/STAFS-early-relfreq-table.tsv --to-tsv

sed -i .bak '1d' feature_tables/STAFS-early-relfreq-table.tsv
sed -i .bak 's/#//' feature_tables/STAFS-early-relfreq-table.tsv
rm feature_tables/*.bak
mkdir cooccurances STAFS_cooccurances STAFS_cooccurances/early

hiorco -k 20 --cutoff 0.001 -o STAFS_cooccurances/early feature_tables/STAFS-early-relfreq-table.tsv

# STAFS active
qiime feature-table filter-samples \
  --i-table feature_tables/MAG_table_final_models-relfreq.qza \
  --m-metadata-file metadata/shotgun_metadata_q2.txt \
  --p-where "[add_facility]='activeSTAFS'" \
  --o-filtered-table feature_tables/STAFS-active-relfreq-table.qza

qiime tools export \
  --input-path feature_tables/STAFS-active-relfreq-table.qza \
  --output-path feature_tables/STAFS-active-relfreq-table
  
mv feature_tables/STAFS-active-relfreq-table/feature-table.biom feature_tables/STAFS-active-relfreq-table.biom
rmdir feature_tables/STAFS-active-relfreq-table

biom convert -i feature_tables/STAFS-active-relfreq-table.biom -o feature_tables/STAFS-active-relfreq-table.tsv --to-tsv

sed -i .bak '1d' feature_tables/STAFS-active-relfreq-table.tsv
sed -i .bak 's/#//' feature_tables/STAFS-active-relfreq-table.tsv
rm feature_tables/*.bak
mkdir STAFS_cooccurances/active

hiorco -k 20 --cutoff 0.001 -o STAFS_cooccurances/active feature_tables/STAFS-active-relfreq-table.tsv

# STAFS advanced
qiime feature-table filter-samples \
  --i-table feature_tables/MAG_table_final_models-relfreq.qza \
  --m-metadata-file metadata/shotgun_metadata_q2.txt \
  --p-where "[add_facility]='advancedSTAFS'" \
  --o-filtered-table feature_tables/STAFS-advanced-relfreq-table.qza

qiime tools export \
  --input-path feature_tables/STAFS-advanced-relfreq-table.qza \
  --output-path feature_tables/STAFS-advanced-relfreq-table
  
mv feature_tables/STAFS-advanced-relfreq-table/feature-table.biom feature_tables/STAFS-advanced-relfreq-table.biom
rmdir feature_tables/STAFS-advanced-relfreq-table

biom convert -i feature_tables/STAFS-advanced-relfreq-table.biom -o feature_tables/STAFS-advanced-relfreq-table.tsv --to-tsv

sed -i .bak '1d' feature_tables/STAFS-advanced-relfreq-table.tsv
sed -i .bak 's/#//' feature_tables/STAFS-advanced-relfreq-table.tsv
rm feature_tables/*.bak
mkdir STAFS_cooccurances/advanced

hiorco -k 20 --cutoff 0.001 -o STAFS_cooccurances/advanced feature_tables/STAFS-advanced-relfreq-table.tsv

## Using co-occurrences for SMETANA - global
# STAFS community size 20 - aerobic
smetana /Users/Zach/Dropbox/PMI_3_analyses/multi-omics_data/shotgun/wrighton_results/Models/*.xml -c STAFS_cooccurances/early/size_20_part_1.tsv -g -v --flavor ucsd --molweight --aerobic --exclude cooccurances/inorganic.txt -o STAFS_cooccurances20_early_aerobic
smetana /Users/Zach/Dropbox/PMI_3_analyses/multi-omics_data/shotgun/wrighton_results/Models/*.xml -c STAFS_cooccurances/active/size_20_part_1.tsv -g -v --flavor ucsd --molweight --aerobic --exclude cooccurances/inorganic.txt -o STAFS_cooccurances20_active_aerobic
smetana /Users/Zach/Dropbox/PMI_3_analyses/multi-omics_data/shotgun/wrighton_results/Models/*.xml -c STAFS_cooccurances/advanced/size_20_part_1.tsv -g -v --flavor ucsd --molweight --aerobic --exclude cooccurances/inorganic.txt -o STAFS_cooccurances20_advanced_aerobic

# STAFS community size 20
smetana Models/*.xml -c STAFS_cooccurances/early/size_20_part_1.tsv -g -v --flavor ucsd --molweight --exclude cooccurances/inorganic.txt -o STAFS_cooccurances20_early
smetana Models/*.xml -c STAFS_cooccurances/active/size_20_part_1.tsv -g -v --flavor ucsd --molweight --exclude cooccurances/inorganic.txt -o STAFS_cooccurances20_active
smetana Models/*.xml -c STAFS_cooccurances/advanced/size_20_part_1.tsv -g -v --flavor ucsd --molweight --exclude cooccurances/inorganic.txt -o STAFS_cooccurances20_advanced

# STAFS active detailed
smetana /Users/Zach/Dropbox/PMI_3_analyses/multi-omics_data/shotgun/wrighton_results/Models/*.xml -c STAFS_cooccurances/active/size_20_part_1.tsv -d -v --flavor ucsd --molweight --aerobic --exclude cooccurances/inorganic.txt -o STAFS_cooccurances20_active_aerobic 

# STAFS advanced detailed
smetana /Users/Zach/Dropbox/PMI_3_analyses/multi-omics_data/shotgun/wrighton_results/Models/*.xml -c STAFS_cooccurances/advanced/size_20_part_1.tsv -d -v --flavor ucsd --molweight --aerobic --exclude cooccurances/inorganic.txt -o STAFS_cooccurances20_advanced_aerobic1 
smetana /Users/Zach/Dropbox/PMI_3_analyses/multi-omics_data/shotgun/wrighton_results/Models/*.xml -c STAFS_cooccurances/advanced/size_20_part_2.tsv -d -v --flavor ucsd --molweight --aerobic --exclude cooccurances/inorganic.txt -o STAFS_cooccurances20_advanced_aerobic2 
smetana /Users/Zach/Dropbox/PMI_3_analyses/multi-omics_data/shotgun/wrighton_results/Models/*.xml -c STAFS_cooccurances/advanced/size_20_part_3.tsv -d -v --flavor ucsd --molweight --aerobic --exclude cooccurances/inorganic.txt -o STAFS_cooccurances20_advanced_aerobic3 
smetana /Users/Zach/Dropbox/PMI_3_analyses/multi-omics_data/shotgun/wrighton_results/Models/*.xml -c STAFS_cooccurances/advanced/size_20_part_4.tsv -d -v --flavor ucsd --molweight --aerobic --exclude cooccurances/inorganic.txt -o STAFS_cooccurances20_advanced_aerobic4 


## Separate tables by decomp stage and calculate co-occurances
# ARF early
qiime feature-table filter-samples \
  --i-table feature_tables/MAG_table_final_models-relfreq.qza \
  --m-metadata-file metadata/shotgun_metadata_q2.txt \
  --p-where "[add_facility]='earlyARF'" \
  --o-filtered-table feature_tables/ARF-early-relfreq-table.qza

qiime tools export \
  --input-path feature_tables/ARF-early-relfreq-table.qza \
  --output-path feature_tables/ARF-early-relfreq-table
  
mv feature_tables/ARF-early-relfreq-table/feature-table.biom feature_tables/ARF-early-relfreq-table.biom
rmdir feature_tables/ARF-early-relfreq-table

biom convert -i feature_tables/ARF-early-relfreq-table.biom -o feature_tables/ARF-early-relfreq-table.tsv --to-tsv

sed -i .bak '1d' feature_tables/ARF-early-relfreq-table.tsv
sed -i .bak 's/#//' feature_tables/ARF-early-relfreq-table.tsv
rm feature_tables/*.bak
mkdir cooccurances ARF_cooccurances ARF_cooccurances/early

hiorco -k 20 --cutoff 0.001 -o ARF_cooccurances/early feature_tables/ARF-early-relfreq-table.tsv

# ARF active
qiime feature-table filter-samples \
  --i-table feature_tables/MAG_table_final_models-relfreq.qza \
  --m-metadata-file metadata/shotgun_metadata_q2.txt \
  --p-where "[add_facility]='activeARF'" \
  --o-filtered-table feature_tables/ARF-active-relfreq-table.qza

qiime tools export \
  --input-path feature_tables/ARF-active-relfreq-table.qza \
  --output-path feature_tables/ARF-active-relfreq-table
  
mv feature_tables/ARF-active-relfreq-table/feature-table.biom feature_tables/ARF-active-relfreq-table.biom
rmdir feature_tables/ARF-active-relfreq-table

biom convert -i feature_tables/ARF-active-relfreq-table.biom -o feature_tables/ARF-active-relfreq-table.tsv --to-tsv

sed -i .bak '1d' feature_tables/ARF-active-relfreq-table.tsv
sed -i .bak 's/#//' feature_tables/ARF-active-relfreq-table.tsv
rm feature_tables/*.bak
mkdir ARF_cooccurances/active

hiorco -k 20 --cutoff 0.001 -o ARF_cooccurances/active feature_tables/ARF-active-relfreq-table.tsv

# ARF advanced
qiime feature-table filter-samples \
  --i-table feature_tables/MAG_table_final_models-relfreq.qza \
  --m-metadata-file metadata/shotgun_metadata_q2.txt \
  --p-where "[add_facility]='advancedARF'" \
  --o-filtered-table feature_tables/ARF-advanced-relfreq-table.qza

qiime tools export \
  --input-path feature_tables/ARF-advanced-relfreq-table.qza \
  --output-path feature_tables/ARF-advanced-relfreq-table
  
mv feature_tables/ARF-advanced-relfreq-table/feature-table.biom feature_tables/ARF-advanced-relfreq-table.biom
rmdir feature_tables/ARF-advanced-relfreq-table

biom convert -i feature_tables/ARF-advanced-relfreq-table.biom -o feature_tables/ARF-advanced-relfreq-table.tsv --to-tsv

sed -i .bak '1d' feature_tables/ARF-advanced-relfreq-table.tsv
sed -i .bak 's/#//' feature_tables/ARF-advanced-relfreq-table.tsv
rm feature_tables/*.bak
mkdir ARF_cooccurances/advanced

hiorco -k 20 --cutoff 0.001 -o ARF_cooccurances/advanced feature_tables/ARF-advanced-relfreq-table.tsv

## Using co-occurrences for SMETANA - global

# ARF community size 20 
smetana Models/*.xml -c ARF_cooccurances/early/size_20_part_1.tsv -g -v --flavor ucsd --molweight --exclude cooccurances/inorganic.txt -o ARF_cooccurances20_early
smetana Models/*.xml -c ARF_cooccurances/active/size_20_part_1.tsv -g -v --flavor ucsd --molweight --exclude cooccurances/inorganic.txt -o ARF_cooccurances20_active
smetana Models/*.xml -c ARF_cooccurances/advanced/size_20_part_1.tsv -g -v --flavor ucsd --molweight --exclude cooccurances/inorganic.txt -o ARF_cooccurances20_advanced

# ARF advanced detailed
smetana /Users/Zach/Dropbox/PMI_3_analyses/multi-omics_data/shotgun/wrighton_results/Models/*.xml -c ARF_cooccurances/advanced/size_20_part_1.tsv -d -v --flavor ucsd --molweight --exclude cooccurances/inorganic.txt -o ARF_cooccurances20_advanced

### NULL MODEL ANALYSIS
# Each facility and stage within had 20 MAGs randomly subsampled from them to test if a random subsample produces signfiicant difference similar to use co-occuring MAGs
python random_selection.py null/CMU-early-relfreq-table.tsv null/CMU-early-rand.tsv
python random_selection.py null/CMU-acive-relfreq-table.tsv null/CMU-active-rand.tsv
python random_selection.py null/CMU-advanced-relfreq-table.tsv null/CMU-advanced-rand.tsv

python random_selection.py null/UTK-early-relfreq-table.tsv null/UTK-early-rand.tsv
python random_selection.py null/UTK-acive-relfreq-table.tsv null/UTK-active-rand.tsv
python random_selection.py null/UTK-advanced-relfreq-table.tsv null/UTK-advanced-rand.tsv

python random_selection.py null/SHSU-early-relfreq-table.tsv null/SHSU-early-rand.tsv
python random_selection.py null/SHSU-acive-relfreq-table.tsv null/SHSU-active-rand.tsv
python random_selection.py null/SHSU-advanced-relfreq-table.tsv null/SHSU-advanced-rand.tsv

# run smetana on random null communitites
cd null/
for i in *-rand.tsv;                                                                                                                                
do
        filename=$(echo $i | cut -d. -f1)
        smetana Models/*.xml -c $i -g -v --flavor ucsd --molweight --exclude cooccurances/inorganic.txt -o $filename
done
