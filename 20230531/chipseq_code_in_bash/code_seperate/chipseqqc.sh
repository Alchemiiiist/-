#!/bin/bash

mkdir 05_chipseq_qc
cd 05_chipseq_qc

echo "chip-seq quality control started at $(date)"

cat ../rawdata/chipseq/config2 | while read id

do
plotFingerprint -b ../03_mapping/${id}_filtered_dedup.bam ../03_mapping/ctcf_chip_input_filtered_dedup.bam \
 --labels CTCF input --plotFileFormat svg --plotTitle CTCF_Enrichment -o CTCF_Enrichment.svg -p 10
done

cd /home/wuhangrui/practice/chipseq

echo "chip-seq quality control finished at $(date)"


### 
# 方法1, 看call到的motif和数据库里记录的motif是不是一致的