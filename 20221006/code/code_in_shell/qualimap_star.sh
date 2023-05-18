cd clean2

mkdir markup
cp -r ./aligndata_star/*_sorted.bam markup

cd markup

echo 'markup+qualimap started at $(date)'
#analysis duplication rate 
#analysis exon rate
work_path=`pwd`
ls $work_path|while read dir
do
( sample=`echo $dir`
picard MarkDuplicates I=${sample} O=${sample}_markdup.bam M=${sample}_markdup_metrics.txt

qualimap rnaseq -bam ${sample} -gtf /home/wuhangrui/database/ref/hg38/gencode.v42.chr_patch_hapl_scaff.annotation.gtf -outdir qualimap/${sample}_rnaseq_qc_results --java-mem-size=20G
)
done
