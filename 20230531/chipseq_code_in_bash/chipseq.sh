#!/bin/bash

########
# code for chipseq
########

### 当前工作路径为
# pwd
# /home/wuhangrui/practice/chipseq

### 创建rawdata的软连接
# 由于下载的rawdata通常存放于组内的nas中，因此可以将目标原文件建立软连接到自己的工作目录下
ln -s /nas_data/whr/practice/chipseq/ /home/wuhangrui/practice/chipseq/rawdata/

### 创建文件夹
mkdir 01_fastqc 02_cutadapt 03_mapping 04_peak_calling

## 读取样本名字
cd ./rawdata/chipseq

ls *_r1.fastq.gz >config # 把样本名字信息存放在config文件
sed -i "s/_r1.fastq.gz//g" config # 去掉后缀，仅保留名字信息
cd /home/wuhangrui/practice/chipseq


## step1 fastqc
## 从nas里直接读取数据进行fastq，速度有点慢，建议在同一局域网下进行或储存条件允许时将数据存放在服务器本地
echo "fastq started at $(date)"

cat rawdata/chipseq/config | while read id
do
fastqc -t 10 -o 01_fastqc/ rawdata/${id}_r1.fastq.gz rawdata/${id}_r2.fastq.gz > 01_fastqc/fastqc.log 2>&1
done

echo "fastq finished at $(date)"

# fastqc [options] fastaq1 [fastq2 ...] 
# -t 指定使用的线程数
# -o 指定输出文件夹


## step2 cutadapt
echo "cutadapt started at $(date)"

cat rawdata/chipseq/config | while read id
do
cutadapt -j 10 \
 --time 1 -e 0.1 -O 3 --quality-cutoff 25 -m 55 \
 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
 -o 02_cutadapt/fixed_${id}_read1.fastq.gz -p 02_cutadapt/fixed_${id}_read2.fastq.gz rawdata/chipseq/${id}_r1.fastq.gz rawdata/chipseq/${id}_r2.fastq.gz > 02_cutadapt/cutadapt.log 2>&1
done

echo "cutadapt finished at $(date)"

# -j 2 设置线程数
# --times 1 只去处一次接头，因为接头只出现一次
# -e 0.1  去除的序列与adaptor相比，missmatch率低于该值，则认为是adaptor，一般设置为0.1 
# -O 3  当与adaptor overlap的碱基数大于等于该值时，才进行去除
# --quality-cutoff 25  小于等于该qulity的碱基去除
# -m 55 trim之后低于该值的一对reads丢弃，一般要大于40
# -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA read1 3'的adaptor
# -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT read2 3'的adaptor
# -o fix.fastq/test_R1_cutadapt.temp.fq.gz read1的输出文件
# -p fix.fastq/test_R2_cutadapt.temp.fq.gz read2的输出文件
# Read1输入文件
# Read2输入文件



## step3 Mapping

### build index

# 建立目录专门存放各种比对软件的index
# cd /home/wuhangrui/database/ref/hg38
mkdir bowtie_index_gencode
cd bowtie_index_gencode

bowtie2-build -f ../GRCh38.p13.genome.fa hg38_human_index > build_index.log 2>&1 
# bowtie2-build [options]* <reference_in> <bt2_base>
# 第二个参数[hg38_human_index]为index前缀
# -f 基因组为fasta格式

cd /home/wuhangrui/practice/chipseq

### mapping
echo "bowtie2 mapping started at $(date)"

cat rawdata/chipseq/config | while read id
do
bowtie2 -x /home/wuhangrui/database/ref/hg38/bowtie_index_gencode/hg38_human_index \
 -1 ./02_cutadapt/fixed_${id}_read1.fastq.gz \
 -2 ./02_cutadapt/fixed_${id}_read2.fastq.gz \
 -p 10 \
 -S ./03_mapping/${id}.sam > ./03_mapping/bowtie2_mapping.log 2>&1
done

echo "bowtie2 mapping finished at $(date)"

# -x index file
# -1 read1 fastq file
# -2 read2 fastq file
# -p threads
# -S out SAM file

## step4 比对结果的排序、过滤、去重、索引
### convert SAM to BAM, and sorted
cat rawdata/chipseq/config | while read id
do
samtools sort -O BAM -o 03_mapping/${id}.bam -@ 10 -m 10G 03_mapping/${id}.sam
done

# 按照染色体坐标排序
# 如果需要构建BAM文件的index文件，必须按照染色体坐标排序

### BAM 文件过滤
cat rawdata/chipseq/config | while read id
do
samtools view -bS -q 30 -@ 10 -m 10G -f 1 -f 2 -o 03_mapping/${id}_filtered.bam 03_mapping/${id}.bam
done

### 去重

echo "deduplicate started at $(date)"
cat rawdata/chipseq/config | while read id
do
picard MarkDuplicates -REMOVE_DUPLICATES true -I ./03_mapping/${id}_filtered.bam \
 -O ./03_mapping/${id}_filtered_dedup.bam -M deduplication.log > picard.log 2>&1
done
echo "deduplicate finished at $(date)"

### 索引
# 比对软件产生的序列通常是随机的。然而，比对后的分析步骤通常要求sam/bam文件被进一步处理，例如在IGV查看比对结果时，常需要输入的bam文件已经被index
echo "index bam files started at $(date)"
cat rawdata/chipseq/config | while read id
do
samtools index -@ 10 ./03_mapping/${id}_filtered_dedup.bam
done
echo "index bam files finished at $(date)"


## step4 call peaks
echo "call peaks started at $(date)"

cp rawdata/chipseq/config rawdata/chipseq/config2
cat rawdata/chipseq/config2 | sed "/input/d" > rawdata/chipseq/config2

cat rawdata/chipseq/config2 | while read id
do
macs2 callpeak -t ./03_mapping/${id}_filtered_dedup.bam \
 -c ./03_mapping/ctcf_chip_input_filtered_dedup.bam \
 -f BAMPE -g hs -n CTCF_ChIP -B -q 0.01 \
 --outdir ./04_peak_calling > $./04_peak_calling/peak_calling.log 2>&1
done
echo "call peaks finished at $(date)"

# -t 处理组的比对文件
# -c 对照组的比对文件
# -f 输入文件的格式, 设为auto则表示自动识别输入文件的格式；双端测序用BAMPE
# -g 有效基因组大小，hs表示人的有效基因组大小，约为2.7E9
# -n name string,会成为输出文件的prefix
# --outdir 输出文件路径
# -B/--BDG 保存堆积对齐片段到bedGraph文件里
# -q FDR Cutoff


# broad peak常用于组蛋白修饰的peak calling，因为组蛋白修饰往往是一大片；
# narrow peak用于转录因子的peak calling，因此本课题应该使用narrow peak calling。