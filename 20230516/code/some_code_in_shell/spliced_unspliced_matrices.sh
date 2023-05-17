#### step0 Constructing spliced and unspliced counts matrices
# we need to have a matrix for spliced and unspliced transcripts
# provide a .gtf to mask repeat regions (这个文件从https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=&db=mm39&hgta_group=allTracks&hgta_track=rmsk&hgta_table=rmsk&hgta_regionType=genome&position=&hgta_outputType=gff&hgta_outFileName=mm10_rmsk.gtf网站生成)
repeats="/home/wuhangrui/database/ref/hg38/hg38_rmsk.gtf"
transcriptome="/home/wuhangrui/database/ref/hg38/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
cellranger_output="/home/wuhangrui/R/data/20230516/rawdata/SRR7722937/sample/sample1"


velocyto run10x -m $repeats \
                $cellranger_output \
                $transcriptome