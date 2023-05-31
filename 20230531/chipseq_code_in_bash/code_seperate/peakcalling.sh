## step4 call peaks
echo "call peaks started at $(date)"

# cp rawdata/chipseq/config rawdata/chipseq/config2
# cat rawdata/chipseq/config2 | sed "/input/d" > rawdata/chipseq/config2

cat rawdata/chipseq/config2 | while read id
do
macs2 callpeak -t ./03_mapping/${id}_filtered_dedup.bam \
 -c ./03_mapping/ctcf_chip_input_filtered_dedup.bam \
 -f BAMPE -g hs -n CTCF_ChIP -B -q 0.05 \
 --outdir ./04_peak_calling > ./04_peak_calling/peak_calling.log 2>&1
done
echo "call peaks finished at $(date)"

# -t 处理组的比对文件
# -c 对照组的比对文件
# -f 输入文件的格式, 设为auto则表示自动识别输入文件的格式；双端测序用BAMPE
# -g 有效基因组大小，hs表示人的有效基因组大小，约为2.7E9
# -n name string,会成为输出文件的prefix
# --outdir 输出文件路径
# -B/--BDG 保存堆积对齐片段到bedGraph文件里
# -q FDR Cutoff


# broad peak常用于组蛋白修饰的peak calling，因为组蛋白修饰往往是一大片；
# narrow peak用于转录因子的peak calling，因此本课题应该使用narrow peak calling。