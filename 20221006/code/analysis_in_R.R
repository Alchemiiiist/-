##################
#第二次作业
#组学分析
##################
#差异分析
library(DESeq2)
library(dplyr)
# heatmap
#install.packages("pheatmap")
#?pheatmap::pheatmap
library("pheatmap")


#load data
control_treat.data = read.table("../../homework2_1006/control_treat.count.txt", header = T, row.names = 1)

#构建分组信息
condition = factor(c(rep("CONTROL",3),rep("TREAT",3)))
coldata = data.frame(row.names = colnames(control_treat.data),condition)
coldata#检查分组信息是否和真实数据一致

#构建dds矩阵
#利用control_treat.count文件（即对象control_treat.data）和分组信息（即对象coldata）构建
dds = DESeqDataSetFromMatrix(countData = control_treat.data,
                             colData = coldata,
                             design = ~condition)

#数据清洗，过滤在所有重复样本中小于1的基因(低质量数据过滤)
dds = dds[rowSums(counts(dds))>1,]


#PCA查看样本间的差异
dds_rld = rlog(dds)
plotPCA(dds_rld, intgroup = "condition")


#差异表达分析
#提取差异分析结果
dds = DESeq(dds) #将数据标准化
resultsNames(dds) #查看结果名称
dds$condition #默认后者的处理组比前面的对照组
res = results(dds)
summary(res) #查看实验对照找到的差异基因分布的情况

table(res$padj < 0.05) # 矫正后的P值,看看有多少差异基因满足所设的P值要求
res = res[order(res$padj),] # 按照pagj的升值排序
resdata = merge(as.data.frame(res), as.data.frame(counts(dds,normalized=TRUE)), by="row.names",sort=FALSE)
write.csv(resdata,file = "../control_treat_dseq.csv")

# 筛选出自己要求的（比较有意义的）差异基因
diff_gene_deseq2 <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(diff_gene_deseq2,file = "../diff_gene.csv") #生成的差异性结果分析文件
# 基因上调与下调的数目
# 并提取显著性差异的基因
upregulate = subset(res, padj < 0.05 & log2FoldChange > 1)
dim(upregulate)
downregulate = subset(res, padj < 0.05 & log2FoldChange < -1)
dim(downregulate)
write.csv(upregulate, "../up_gene.csv")
write.csv(downregulate, "../down_gene.csv")


# 差异表达可视化


#筛选出来的所有差异基因分析可视化
# 提取经过deseq2归一化后的表达矩阵
counts.norm.data = resdata[,-2:-7]
# 将行名改为基因名
rownames(counts.norm.data) = counts.norm.data[,1]
counts.norm.data = counts.norm.data[,-1]
# 在原矩阵中寻找差异表达基因，原表达矩阵的基因和筛选出的差异基因交集
df = counts.norm.data[intersect(rownames(counts.norm.data),rownames(diff_gene_deseq2)),]
df2 = as.matrix(df)
pheatmap(df2,
         show_rownames = F,
         show_colnames = T,
         cluster_cols = F,
         cluster_rows = T,
         height = 10,
         scale = "row",
         fontsize = 10,
         angle_col = 45,
         clustering_method = "single",
)

#最显著高表达和最显著地表达的20个基因差异分析可视化
# 提取treat组最显著高表达和最低表达的基因,并排序
up_gene.data = as.data.frame(upregulate)
up_gene.data = up_gene.data[order(up_gene.data$padj),]
down_gene.data = as.data.frame(downregulate)
down_gene.data = down_gene.data[order(down_gene.data[,6]),]
# 提取排序后最显著表达的20个基因
up_gene.data20 = up_gene.data[1:20,]
down_gene.data20 = down_gene.data[1:20,]
regulate_gene40 = rbind(up_gene.data20,down_gene.data20)
df40 = counts.norm.data[intersect(rownames(counts.norm.data),rownames(regulate_gene40)),]
pheatmap(df40,
         show_rownames = T,
         show_colnames = T,
         cluster_cols = F,
         cluster_rows = T,
         height = 10,
         scale = "row",
         fontsize = 10,
         angle_col = 45,
         clustering_method = "single",
)

##对差异显著的基因进行GO富集分析
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db") #人类注释数据库
BiocManager::install("enrichplot")
BiocManager::install("GOplot")
library("clusterProfiler")
library("org.Hs.eg.db")
library("GOplot")

diff_gene = as.data.frame(diff_gene_deseq2)
diff_gene["ENSEMBL"] = rownames(diff_gene)
#转化基因id
gene.df = bitr(diff_gene$ENSEMBL,
               fromType = "ENSEMBL",
               toType = c("ENTREZID"),
               OrgDb = org.Hs.eg.db)
enrich = merge(diff_gene,
               gene.df,
               by.x = "ENSEMBL",
               by.y = "ENSEMBL",
               all = TRUE)
#过滤ENTRIZID中有缺失值的行
enrich = enrich[complete.cases(enrich[,8]),]

enrichment = enrich %>%
  as.data.frame()%>%
  dplyr::select(ENSEMBL,ENTREZID,log2FoldChange)

enrichment = enrichment[!duplicated(enrichment[,1]),]

gene = enrichment$ENTREZID

#GO富集分析
ego = enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db,
               pvalueCutoff = 0.01,
               pAdjustMethod = "BH",
               qvalueCutoff = 0.05,
               ont = "all",
               readable = T)

barplot(ego, showCategory = 5, color = "pvalue")
dotplot(ego, showCategory = 5)
cnetplot(ego, showCategory = 3, circular = TRUE)


# 引入差异倍数
go = data.frame(Category = ego$ONTOLOGY,
                ID = ego$ID,
                Term = ego$Description,
                Genes = gsub("/", ", ", ego$geneID),
                adj_pval = ego$p.adjust)

# 将ENSEMBL转化为symbol，与ego的geneID名对应
names = bitr(enrichment$ENSEMBL,
               fromType = "ENSEMBL",
               toType = c("SYMBOL"),
               OrgDb = org.Hs.eg.db)
enrichment2 = merge(enrichment, names,
                    by.x = "ENSEMBL",
                    by.y = "ENSEMBL",
                    )
enrichment2 = enrichment2[complete.cases(enrichment2[,3]),]


genelist = data.frame(ID = enrichment2$SYMBOL,
                      logFC = enrichment2$log2FoldChange)

#获取
circ = circle_dat(go, genelist)
termNum = 5
geneNum = nrow(genelist)
chord = chord_dat(circ, genelist[1:geneNum,],go$Term[1:termNum])

GOChord(chord,
        space = 0.02,        
        gene.order = 'logFC',    
        gene.space = 0.25,      
        gene.size = 2)

GOBubble(circ, title = 'Bubble_plot',
         colour = c('skyblue', 'pink', 'red'),
         display = 'multiple', labels = 3)

GOHeat(chord, nlfc = 1, fill.col = c('gold', 'white', 'purple'))
#只用显著高表达或低表达的40个基因
#regulate_gene40_2 = regulate_gene40
#regulate_gene40_2["ENSEMBL"] = rownames(regulate_gene40)
#symbol = bitr(regulate_gene40_2$ENSEMBL,
             #fromType = "ENSEMBL",
             #toType = c("SYMBOL"),
             #OrgDb = org.Hs.eg.db)
#enrichment_40 = merge(regulate_gene40_2, symbol,
                    #by.x = "ENSEMBL",
                    #by.y = "ENSEMBL",
#)
#genelist_40 = data.frame(ID = enrichment_40$SYMBOL,
                      #logFC = enrichment_40$log2FoldChange)
#circ_40 = circle_dat(go, genelist_40)
#geneNum_40 = nrow(genelist_40)
#chord = chord_dat(circ_40, genelist[1:geneNum_40,],go$Term[1:termNum])
#GOChord(chord,
       # space = 0.02,        
        #gene.order = 'logFC',    
        #gene.space = 0.25,      
        #gene.size = 2)
