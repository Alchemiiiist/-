cat /home/wuhangrui/R/data/20230516/rawdata/SRR_Acc_List.txt |while read i
do
/home/wuhangrui/installApp/sratools/sratoolkit.3.0.5-ubuntu64/bin/fastq-dump --gzip --split-files -A $i ${i}.sra && echo "** ${i}.sra to fastq done **"
done