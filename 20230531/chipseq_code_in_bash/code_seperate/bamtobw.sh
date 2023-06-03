#!/bin/bash
# 使用bamCoverage进行 bam -> bigwig 文件转换（一种压缩方式：将基因组压缩为很多个bin区域）
# bigwig 文件记录的是“每个碱基”的coverage情况，可以理解为ChIP-Seq的信号
cd /home/wuhangrui/practice/chipseq/06_downstream
mkdir co_location

# 将input和rep1都进行bam -> bw 文件的转换

cat /home/wuhangrui/practice/chipseq/rawdata/chipseq/config | while read id
do
bamCoverage --bam /home/wuhangrui/practice/chipseq/03_mapping/${id}_filtered_dedup.bam \
 -o ./co_location/${id}_filtered_dedup.bw \
 --binSize 10 \
 --normalizeUsing RPGC \
 --effectiveGenomeSize 2913022398 \
 --extendReads \
 -p 10
done