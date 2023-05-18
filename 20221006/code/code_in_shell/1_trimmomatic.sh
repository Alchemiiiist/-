mkdir clean
cd clean

find /home/wuhangrui/rawdata/20221024/01.RawData -name "*fq.gz" |grep 1.fq.gz | sort > 1.txt
find /home/wuhangrui/rawdata/20221024/01.RawData -name "*fq.gz" |grep 2.fq.gz | sort > 2.txt
paste 1.txt 2.txt > config1

#mkdir {fastqc,multiqc_result,fastqc_filtered,multiqc_result_filtered}

echo "trimmomatic cut adapters starts at $(date)"
cat config1 |while read id
do
arr=(${id})
fq1=${arr[0]}
fq2=${arr[1]}
filename=`echo $(basename $fq1 _1.fq.gz)`
trimmomatic PE -threads 20 -phred33 $fq1 $fq2 \
-baseout ${filename}.fq.gz \
ILLUMINACLIP:/home/wuhangrui/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:40 HEADCROP:13
echo "trimomatic cut adapters finished at $(date)"
done
