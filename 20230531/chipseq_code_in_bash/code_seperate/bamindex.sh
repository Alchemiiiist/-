#!/bin/bash

cd /home/wuhangrui/practice/chipseq
### 索引
# 比对软件产生的序列通常是随机的。然而，比对后的分析步骤通常要求sam/bam文件被进一步处理，例如在IGV查看比对结果时，常需要输入的bam文件已经被index
echo "index bam files started at $(date)"
cat rawdata/chipseq/config | while read id
do
samtools index -@ 10 ./03_mapping/${id}_filtered_dedup.bam
done
echo "index bam files finished at $(date)"