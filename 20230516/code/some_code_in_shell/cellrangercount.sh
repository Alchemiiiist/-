cellranger count --id=sample1 \
                   --transcriptome=/home/wuhangrui/database/ref/hg38/refdata-gex-GRCh38-2020-A \
                   --fastqs=/home/wuhangrui/R/data/20230516/rawdata/SRR7722937 \
                   --sample=SRR7722937 \
                   --nosecondary
# id指定输出文件存放目录名
# transcriptome指定与CellRanger兼容的参考基因组
# fastqs指定mkfastq或者自定义的测序文件
# sample要和fastq文件的前缀中的sample保持一致，作为软件识别的标志,S1
# expect-cells指定复现的细胞数量，这个要和实验设计结合起来
# nosecondary 只获得表达矩阵，不进行后续的降维、聚类和可视化分析(因为后期会自行用R包去做)