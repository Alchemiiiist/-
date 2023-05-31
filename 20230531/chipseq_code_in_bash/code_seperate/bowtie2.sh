#!/bin/bash
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