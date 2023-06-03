#!/bin/bash


cd /home/wuhangrui/practice/chipseq

mkdir 06_downstream


##### 1. call motif, peak注释和notif search ########
prefix=CTCF_ChIP # 该prefix与peak calling的prefix相同
### 首先提取用于寻找motif的fasta文件
## 取峰值前后100个碱基，第二列-50得到第六列，第三列+49得到第七列，print第1（染色体列），6，7，4列
awk -v OFS="\t" '{$6=$2-50;$7=$3+49;print $1','$6','$7','$4}' ./04_peak_calling/${prefix}_summits.bed > ./06_downstream/${prefix}_motif.bed
bedtools getfasta -fi /home/wuhangrui/database/ref/hg38/GRCh38.p13.genome.fa -bed ./06_downstream/${prefix}_motif.bed > ./06_downstream/${prefix}_motif.fasta

cd /home/wuhangrui/practice/chipseq/06_downstream

echo "calling motif started at $(date)"

mkdir memechip_out

meme-chip -db /home/wuhangrui/database/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme -meme-p 10 -dna ${prefix}_motif.fasta -oc /home/wuhangrui/practice/chipseq/06_downstream/memechip_out

echo "calling motif finished at $(date)"
cd /home/wuhangrui/practice/chipseq

##### 2. peakannotation #####
cd /home/wuhangrui/practice/chipseq/06_downstream

echo "peak annotation started at $(date)"
mkdir peak_annotation
annotatePeaks.pl /home/wuhangrui/practice/chipseq/04_peak_calling/${prefix}_peaks.narrowPeak hg38 > peak_annotation/peak_annotation.tsv
echo "peak_annotation finished at $(date)"

## 之后可以使用peak_annotation.tsv文件丢在R里做GO分析等。


##### 3. 不适合看motif的共定位分析 #######
##### 3.1 与TSS, TES的共定位情况 #####
# 使用bamCoverage进行 bam -> bigwig 文件转换（一种压缩方式：将基因组压缩为很多个bin区域）
# bigwig 文件记录的是“每个碱基”的coverage情况，可以理解为ChIP-Seq的信号
cd /home/wuhangrui/practice/chipseq/06_downstream
mkdir co_location

# 将input和rep1都进行bam -> bw 文件的转换

cat /home/wuhangrui/practice/chipseq/rawdata/config | while read id
do
bamCoverage --bam /home/wuhangrui/practice/chipseq/03_mapping/${id}_filtered_dedup.bam \
 -o ./co_location/${id}_filtered_dedup.bw \
 --binSize 10 \
 --normalizeUsing RPGC \
 --effectiveGenomeSize 2913022398 \
 --extendReads \
 -p 10
done

# --bam输入的bam文件
# -o 输出的bigwig文件
# --binSize 计算coverage时的bin的大小
# --normalizeUsing 选择normalization的方法 RPGC按照bin和测序量大小归一化
# --effectiveGenomesize 如果需要normalization则需有设置有效基因组大小
# --extendReads 设置该值以拓展reads的长度，以计算出用于测序的DNA片段的长度；pair-end的测序可以自动估计DNA片段的长度
# -p 使用的线程数



### TSS bed 文件准备
#grep "+" uscs_refseq.bed > uscs_refseq.plus.bed
#awk '{$7=$2-2500;$8=$2+2500;print $1,$7,$8,$4,$5,$6}' uscs_refseq.plus.bed > uscs_refseq.plus.tss.bed
#grep "-" uscs_refseq.bed > uscs_refseq.minor.bed
#awk '{$7=$3-2500;$8=$3+2500;print $1,$7,$8,$4,$5,$6}' uscs_refseq.minor.bed > uscs_refseq.minor.tss.bed
#cat uscs_refseq.plus.tss.bed uscs_refseq.minor.tss.bed > uscs_refseq.tss.bed
#sort -n -k2 uscs_refseq.tss.bed >uscs_refseq.tss.sorted.bed


##### 与其他区域的共定位情况（高甲基化区域，低甲基化区域，5'UTR，3'UTR等）
# 使用“computeMatrix”计算指定位置的信号矩阵
cd /home/wuhangrui/practice/chipseq/06_downstream

computeMatrix reference-point \
 -S ./co_location/ctcf_chip_input_filtered_dedup.bw ./co_location/ctcf_chip_rep1_filtered_dedup.bw ./co_location/ctcf_chip_rep2_filtered_dedup.bw \
 -R /home/wuhangrui/database/ref/hg38/uscs_refseq.bed \
 -a 2500 -b 2500 \
 --samplesLabel input ctcf_rep1 ctcf_rep2 \
 --sortRegions descend \
 -o ./co_location/Rep1_Input_TSS_Matrix.gz \
 -p 10 > computeMatrix_input.log 2>&1

# 使用reference-point模式，reference位点为TSS位置
# -S 输入的bigwig文件
# -R 记录reference位点信息的bed文件
# -a 起始位点上游的距离
# -b 起始位点下游的距离
# --sortRegions 在矩阵中，信号按什么顺序排列？
# --samplesLabel 在矩阵中，样品如何命名
# -o 输出的矩阵名
# -p 线程数

# 使用plotHeatmap作图
cd /home/wuhangrui/practice/chipseq/06_downstream
plotHeatmap -m ./co_location/Rep1_Input_TSS_Matrix.gz \
 -out ./co_location/Rep1_Input_TSS_Matrix.svg \
 --plotFileFormat svg --yAxisLabel RPGC --regionsLabel all_tss \
 --legendLocation none