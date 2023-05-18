cd clean2

mkdir htseq
mkdir salmon_quant

echo "htseq+salmon started at $(date)"

find /home/wuhangrui/rawdata/20221024/01.RawData -name "*fq.gz" |grep 1.fq.gz | sort > 1.txt
find /home/wuhangrui/rawdata/20221024/01.RawData -name "*fq.gz" |grep 2.fq.gz | sort > 2.txt
paste 1.txt 2.txt > config1

cat config1 |while read id
do
arr=(${id})
fq1=${arr[0]}
fq2=${arr[1]}
filename=`echo $(basename $fq1 _1.fq.gz)`

fq11=${filename}_1P.fq
fq22=${filename}_2P.fq

nohup htseq-count -f bam -r name -s no -a 10 -t exon -i gene_id -m union ./aligndata_star/${filename}_sorted.bam /home/wuhangrui/database/ref/hg38/gencode.v42.chr_patch_hapl_scaff.annotation.gtf > ./htseq/${filename}_star_htseq_count.txt 2>log &


cd salmon_quant
salmon quant --gcBias -i /home/wuhangrui/database/ref/hg38/salmon_index_gencode -l A -1 ../$fq11 -2 ../$fq22 -g /home/wuhangrui/database/ref/hg38/gencode.v42.chr_patch_hapl_scaff.annotation.gtf -o ./${filename}_salmon_quant
#salmon quant -i /home/wuhangrui/database/ref/hg38/salmon_index_gencode -l A -1 ../$fq11 -2 ../$fq22 -g /home/wuhangrui/database/ref/hg38/Homo_sapines.GRCh38.108.chr_patch_hapl_scaff.annotation.gtf -o ./${filename}_salmon_quant

awk '$4>1{print$0}' ${filename}_salmon_quant/quant.genes.sf > ${filename}_tpm.txt
wc -l ${filename}_tpm.txt >> tpm.txt
cd ../

done

echo "htseq+salmon finished at $(date)"
