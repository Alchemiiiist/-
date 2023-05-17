cat /home/wuhangrui/R/data/20230516/rawdata/SRR_Acc_List.txt | while read i
do
mv ${i}_1*.gz ${i}_S1_L001_I1_001.fastq.gz
mv ${i}_2*.gz ${i}_S1_L001_R1_001.fastq.gz
mv ${i}_3*.gz ${i}_S1_L001_R2_001.fastq.gz
done