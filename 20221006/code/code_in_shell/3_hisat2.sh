cd clean

mkdir aligndata_hisat2

echo "hisat2 started at $(date)"

cat config1  |while read id
do
arr=(${id})
fq1=${arr[0]}
fq2=${arr[1]}
filename=`echo $(basename $fq1 _1.fq.gz)`

fq11=${filename}_1P.fq.gz
fq22=${filename}_2P.fq.gz

hisat2 -p 20 -x /home/wuhangrui/database/ref/hg38/hisat2/grch38/genome -1 $fq11 -2 $fq22 -S ${filename}.sam  
 
mv *.sam aligndata_hisat2
cd aligndata_hisat2

samtools view -S ${filename}.sam -b > ${filename}.bam

samtools sort ${filename}.bam -o ${filename}_sorted.bam

samtools index ${filename}_sorted.bam

cd ../

done
echo "hisat2 finished at $(date)"


