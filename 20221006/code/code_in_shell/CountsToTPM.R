### 提取count表达矩阵
#load library
library("DESeq2")
library('tximeta')
library("GEOquery")
library("GenomicFeatures")
library("BiocParallel")
library("GenomicAlignments")
library("dplyr")
#########################################
##1. 提取表达矩阵
#generate a bam list
sample_name = readLines("./sorted_data/namelist.txt")
dir = "./sorted_data"
fls = file.path(dir,paste0(sample_name, "_sorted.bam"))
file.exists(fls)
bamlist = BamFileList(fls, yieldSize = 2000000)

#generate features
txdb = makeTxDbFromGFF("/home/wuhangrui/database/ref/hg38/Homo_sapiens.GRCh38.108.chr_patch_hapl_scaff.gtf")
exonsByGene = exonsBy(txdb, by = "gene")


#generate SE project
register(MulticoreParam(workers = 16))
SE = summarizeOverlaps(features = exonsByGene, reads = bamlist, mode = "Union",
                       singleEnd = F, ignore.strand = T)

#add coldata
colnames(SE) = sub(colnames(SE), pattern = ".bam", replacement = '')

#savedata
saveRDS(SE, file = "./SE.rds")

#提取表达矩阵
express_matrix = assay(SE)

##将counts表达矩阵导出
write.csv(express_matrix,"./count_matrix.csv")

#######################################
##2. 将表达矩阵的行名由(ensembl)转化为(symbol)
#导入注释文件
human_anno = rtracklayer::import("/home/wuhangrui/database/ref/hg38/Homo_sapiens.GRCh38.108.chr_patch_hapl_scaff.gtf")
human_anno = rtracklayer::as.data.frame(human_anno)

id_name = human_anno[,c("gene_id","gene_name")]
id_name = unique(id_name)

#change gene_id to gene_name
rownames(id_name) = id_name$gene_id

#按照express_matrix的顺序排序
id_name = id_name[rownames(express_matrix),]
id_name = id_name[,2]

#将表达矩阵的行名由gene_id更改为gene_name
express_matrix_genename = express_matrix
rownames(express_matrix_genename) = id_name

######################################
##3. 将counts表达矩阵归一化得到Tpm
#获取基因外显子长度
# get total exon length
exons_gene_lens <- lapply(exonsByGene,function(x){sum(width(reduce(x)))})
# 将外显子长度转化为data.frame
gene_length = sapply(exons_gene_lens,function(x){x})
id_length = as.data.frame(gene_length)

#进行归一化——得到TPM表达矩阵
##TPM calculation
kb = id_length$gene_length / 1000
rpk = express_matrix_genename / kb

calcu_tpm1 = function(rpkdata){
  rpk = as.data.frame(rpkdata)
  tpm_tem = rpk %>%
    mutate_if(is.numeric, funs(./sum(.)))
  tpm = tpm_tem*1000000
  # tpm = log(tpm)
  return(tpm)
}

tpm1 = calcu_tpm1(rpk)

#another function to calculate tpm
calcu_tpm2 = function(data){
  data = as.data.frame(data)
  total_data = data %>%
    summarise_all(sum)
  for(i in 1:ncol(data)){
    data[,i] = (data[,i])/total_data[,i]*1e6
  }
  # data <- log(data)
  return(data)
}

tpm2 = calcu_tpm2(rpk)

##将TPM表达矩阵导出
write.csv(tpm1, "./tpm_matrix.csv")
