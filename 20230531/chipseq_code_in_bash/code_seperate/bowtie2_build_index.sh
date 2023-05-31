## step3 Mapping
### build index

# 建立目录专门存放各种比对软件的index
cd /home/wuhangrui/database/ref/hg38
mkdir bowtie_index_gencode
cd bowtie_index_gencode

bowtie2-build -f ../GRCh38.p13.genome.fa hg38_human_index > build_index.log 2>&1 
# bowtie2-build [options]* <reference_in> <bt2_base>
# 第二个参数[hg38_human_index]为index前缀
# -f 基因组为fasta格式

cd /home/wuhangrui/practice/chipseq