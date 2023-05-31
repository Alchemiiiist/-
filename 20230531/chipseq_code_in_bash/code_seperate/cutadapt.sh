#!/bin/bash
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