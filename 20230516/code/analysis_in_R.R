###


library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyverse)
###### 0. 创建Seurat对象##########
# load the dataset 读入cellranger产生的标准输出文件
# data.dir Directory - containing the matrix.mtx, gene.tsv, barcodes.tsv files provided by 10x
# gene.column - gene列在Dimnames中的位置
# unique.features - Make feature names unique
pro_data <- Read10X(
  data.dir = "/home/wuhangrui/R/data/20230516/rawdata/SRR7722937/sample/sample1/outs/filtered_feature_bc_matrix",
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE # 重复的基因名自动区分
)

# 使用原始数据创建Seurat对象
# min.cell 一个基因最少在多少个细胞中表达
# min.features 一个细胞中最少表达多少个基因
# 筛除在少于20个细胞中表达的基因
# 筛除表达量小于300个基因的细胞
pro <- CreateSeuratObject(
  counts = pro_data,
  project = "SRR7722937",
  min.cells = 20,
  min.features = 300)


# 查看barcode的属性信息
# nCount表示UMI数量，nfeature表示捕获基因的数量
head(pro@meta.data)
# 在metadata中增加一列barcode属性,属性明CB
#pro[["CB"]] <- rownames(pro@meta.data)
#pro@meta.data$sample <- str_replace(pro@meta.data$CB,"_[ATGC]{16}.*$","")
#pro@meta.data$animal <- str_replace(pro@meta.data$sample,"_.*$","")
#pro@meta.data$condition <- str_replace(pro@meta.data$sample,"^.*_","")

##### 1. 质控###############
# 使用PercentageFeatureSet 对某些基因表达量的和取百分比，得到线粒体基因所占的比例
pro[["percent.mt"]] <- PercentageFeatureSet(
  pro,
  pattern = "^MT-"
)

#head(pro@meta.data)

# 可视化QC需要的可视化指标
VlnPlot(
  pro,
  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
  ncol = 3,
  pt.size = 0.002
)
# 可视化barcode属性中两个特征的相关关系
#plot1 <- FeatureScatter(pro, feature1 = "nCount_RNA",feature2 = "percent.mt")
#plot2 <- FeatureScatter(pro, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1+plot2

# 质控
pro <- subset(pro, subset = nFeature_RNA > 1000 & percent.mt < 20)
VlnPlot(
  pro,
  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
  ncol = 3,
  pt.size = 0
)
plot1 <- FeatureScatter(pro, feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2 <- FeatureScatter(pro, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2


#####2. 归一化(对数)#################################
# 使用NormalizeData()进行归一化
pro <- NormalizeData(pro, normalization.method = "LogNormalize", scale.factor = 10000)
# 将归一化结果储存在pro@assays$RNA@data
head(pro@assays$RNA@data)
dim(pro@assays$RNA@data)


#####3. 寻找高变基因################################
#高变基因：在不同细胞中表达差异大的基因
pro <- FindVariableFeatures(pro, selection.method = "vst", nfeatures = 4000)
#结果存放在pro@assays$RNA@var.features
head(pro@assays$RNA@var.features)
# 对高变基因可视化
plot3 <- VariableFeaturePlot(pro)
plot3

##### 4.标准化##################################
# 对上述4000个高变基因进行标准化
pro <- ScaleData(pro)
# 结果存放在pro@assays$RNA@scale.data
#dim(pro@assays$RNA@scale.data)

##### 5.PCA 降维################################
pro <-RunPCA(pro,npcs=50)
# 降维结果存储在pro@meta.data@pca降维结果可视化
DimPlot(pro,reduction = "pca")
#DimHeatmap(pro,dims = 1:2, cells = 500, balanced =TRUE)



##### 6.选取前20个PCs对样品进行非线性降维
# 选择非线性降维所使用的PCA维数
ElbowPlot(pro)
# 前20维度能解释几乎所有的变异，选择20维
pro <- RunUMAP(pro, dims = 1:20)
Embeddings(pro,"umap") %>% head()
Embeddings(pro,"umap") %>% as.data.frame() %>%
  ggplot(aes(x=UMAP_1,y=UMAP_2)) + geom_point()



##### 7.聚类######################################
pro <- FindNeighbors(pro,dims = 1:20)
# resolution 越高聚类越多
pro <- FindClusters(pro,resolution = 0.5)
# 聚类信息储存在pro@mata.data$seurat_clusters中
table(pro@meta.data$seurat_clusters)
# cluster标签信息已经赋予每个细胞
head(pro@meta.data)
# 聚类可视化
DimPlot(pro, reduction = "umap", pt.size=0.2) 


##### 8.将seruat对象转化为anndata格式（用于scVelo）#################
###导出相应的矩阵
# save metadata table
pro$UMAP_1 <- pro@reductions$umap@cell.embeddings[,1]
pro$UMAP_2 <- pro@reductions$umap@cell.embeddings[,2]
pro$barcode <- colnames(pro)
head(pro@meta.data)
write.csv(pro@meta.data, file='metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
count_matrix <- GetAssayData(pro, assay="RNA", slot="counts")
writeMM(count_matrix, file="/home/wuhangrui/R/data/20230516/analysis_in_R/counts.mtx")

# write dimensionality reduction matrix
write.csv(pro@reductions$pca@cell.embeddings, file="pca.csv", quote=F, row.names=F)

# write gene names
write.table(
  data.frame("gene"=rownames(count_matrix)),file="gene_names.csv",
  quote=F,row.names=F,col.names=F 
)

### 接下来的步骤在python环境中进行
