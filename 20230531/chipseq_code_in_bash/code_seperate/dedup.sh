### 去重
cat rawdata/chipseq/config | while read id
do
java -jar /home/wuhangrui/installApp/picard/build/libs/picard.jar MarkDuplicates -REMOVE_DUPLICATES true -I ./03_mapping/${id}_filtered.bam \
 -O ./03_mapping/${id}_filtered_dedup.bam -M deduplication.log > ./03_mapping/picard.log 2>&1
done