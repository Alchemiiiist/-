#!/bin/bash
## step1 fastqc
cat rawdata/chipseq/config | while read id
do
fastqc -t 10 -o 01_fastqc/ rawdata/chipseq/${id}_r1.fastq.gz rawdata/chipseq/${id}_r2.fastq.gz > 01_fastqc/fastqc.log 2>&1
done

# fastqc [options] fastaq1 [fastq2 ...] 
# -t 指定使用的线程数
# -o 指定输出文件夹