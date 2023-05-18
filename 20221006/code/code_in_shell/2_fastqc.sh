cd clean

mkdir {fastqc,multiqc_result,fastqc_filtered,multiqc_result_filtered}
echo "fastqc started at $(date)"

cat config1 |while read id
do
arr=(${id})
fq1=${arr[0]}
fq2=${arr[1]}
filename=`echo $(basename $fq1 _1.fq.gz)`

fq11=${filename}_1P.fq.gz
fq22=${filename}_2P.fq.gz
fastqc $fq1 $fq2 -q --noextract -t 20 -o fastqc
fastqc $fq11 $fq22 -q --noextract -t 20 -o fastqc_filtered

done

multiqc fastqc -o multiqc_result
multiqc fastqc_filtered -o multiqc_result_filtered

echo "fastqc finished at $(date)"
