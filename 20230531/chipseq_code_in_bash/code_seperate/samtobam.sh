#!/bin/bash

## step4 比对结果的排序、过滤、去重、索引

echo "start to convert sam to bam, sort, and filtered at $(date)"
### convert SAM to BAM, and sorted
cat rawdata/chipseq/config | while read id
do
samtools sort -O BAM -o ./03_mapping/${id}.bam -@ 10 -m 10G 03_mapping/${id}.sam
done

# 按照染色体坐标排序
# 如果需要构建BAM文件的index文件，必须按照染色体坐标排序

### BAM 文件过滤
cat rawdata/chipseq/config | while read id
do
samtools view -bS -q 30 -@ 10 -m 10G -f 1 -f 2 -o ./03_mapping/${id}_filtered.bam 03_mapping/${id}.bam
done

echo "sam to bam finished at $(date)"