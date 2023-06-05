#!/bin/bash
set -x 
#set -e
curl -O http://ftp.microbio.me/emp/release1/mapping_files/emp_qiime_mapping_release1.tsv
ctx=Deblur_2021.09-Illumina-16S-V4-150nt-ac8c0b
for asv in $(cut -f 3 ASVs_repseq.txt | grep -v 150)
do 
    echo $asv
    mkdir -p results/${asv}
    echo ${asv} | \
        redbiom search features --context $ctx --min-count 100 > results/${asv}/samples.ids
    if [[ $(wc -l < "results/${asv}/samples.ids") -eq 0 ]];
    then
        continue
    fi

    cat results/${asv}/samples.ids | \
        grep -v "_714." | \
        redbiom fetch sample-metadata \
            --force-category env_package \
            --force-category empo_1 \
            --force-category empo_2 \
            --force-category empo_3 \
            --output results/${asv}/samples.metadata.tsv
done

python summarize.py
