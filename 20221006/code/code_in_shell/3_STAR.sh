cd clean2

mkdir aligndata_star

echo "STAR started at $(date)"

#STAR \
#--runThreadN 40 \
#--runMode genomeGenerate \
#--genomeDir /home/wuhangrui/database/ref/hg38/star \
#--genomeFastaFiles /home/wuhangrui/database/ref/hg38/GRCh38.p13.genome.fa \
#--sjdbGTFfile /home/wuhangrui/database/ref/hg38/gencode.v42.chr_patch_hapl_scaff.annotation.gtf \
#--sjdbOverhang 149 \

for i in *_1P.fq

do
filename=`echo $(basename $i _1P.fq)`

fq11=${filename}_1P.fq
fq22=${filename}_2P.fq

STAR \
--runThreadN 40 \
--genomeDir /home/wuhangrui/database/ref/hg38/star \
--readFilesIn $fq11 $fq22 \
--sjdbOverhang 149 \
--outFileNamePrefix ./aligndata_star/${filename} 

cd aligndata_star

samtools view -S ${filename}Aligned.out.sam -b > ${filename}.bam

samtools sort ${filename}.bam -o ${filename}_sorted.bam

samtools index ${filename}_sorted.bam

cd ../

done
echo "STAR finished at $(date)"


