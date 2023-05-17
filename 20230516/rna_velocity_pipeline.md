***
<font face = "times new roman" size = 3><center>
scRNA-sequencing  
RNA velocity analysis  

</center>
</font>

***
#### 0. 准备工作，配置scRNA分析的环境
- 数据下载工具-sra-tools
- 10X 单细胞上游分析工具 cellranger
- RNA velocity分析工具 scvelo
- R 环境中需要Seurat及其依赖包
  - 最好单独创建一个conda环境
  - 使用源码下载上述软件时，可以添加至环境变量中
  - cellranger除了官方requirements，还需要samtools软件进行内置的比对。
  - 其他依赖请见官方文档
    - https://github.com/ncbi/sra-tools/wiki
    - https://github.com/theislab/scvelo
    - https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation
    - https://satijalab.org/seurat/
#### 1. 下载数据集
- 该tutorial以SRR7722937为例
```
### 下载内容较大，可用tmux放在后台进行
conda install sra-tools
### SRR_Acc_List 可以在数据库中找到SRA生成并下载
prefetch --option-file SRR_Acc_List.txt

### 将下载的sra文件转为fastq文件
fastq-dump --gzip --split-files SRR7722937
```


#### 2. cellranger分析

- sra文件是从数据库中下载的原始数据，首先从该文件中提取出相应的序列信息
```
## sra2fastq.sh
## 将下载的sra文件进行拆分，得到fastq文件，包含reads和index(barcode)信息
cat /home/wuhangrui/R/data/20230516/rawdata/SRR_Acc_List.txt |while read i
do
/home/wuhangrui/installApp/sratools/sratoolkit.3.0.5-ubuntu64/bin/fastq-dump --gzip --split-files -A $i ${i}.sra && echo "** ${i}.sra to fastq done **"
done

## 改名，改成cellranger可识别的标准格式 
cat /home/wuhangrui/R/data/20230516/rawdata/SRR_Acc_List.txt | while read i
do
mv ${i}_1*.gz ${i}_S1_L001_I1_001.fastq.gz
mv ${i}_2*.gz ${i}_S1_L001_R1_001.fastq.gz
mv ${i}_3*.gz ${i}_S1_L001_R2_001.fastq.gz
done
```
- 使用cellranger进行定量，得到相应的bam文件和标准输出文件（表达矩阵，barcodes，features）等一系列outputs
```
### 比对时间可能较长，可使用tmux放在后台进行
cellranger count --id=sample1 \
                   --transcriptome=/home/wuhangrui/database/ref/hg38/refdata-gex-GRCh38-2020-A \
                   --fastqs=/home/wuhangrui/R/data/20230516/rawdata/SRR7722937 \
                   --sample=SRR7722937 \
                   --nosecondary
# id指定输出文件存放目录名
# transcriptome指定与CellRanger兼容的参考基因组
# fastqs指定mkfastq或者自定义的测序文件
# sample要和fastq文件的前缀中的sample保持一致，作为软件识别的标志,S1
# expect-cells指定复现的细胞数量，这个要和实验设计结合起来
# nosecondary 只获得表达矩阵，不进行后续的降维、聚类和可视化分析(因为后期会自行用R包去做)
```
- 该步骤完成后会得到一个比较复杂的output结构：
  - 包含比对信息--bam文件--实际上调用了samtools工具进行比对
  - 包含表达矩阵信息，barcode信息，和基因信息，这三者是下游分析所使用的标准输出文件
  - ![](2023-05-17-12-58-08.png)
#### 3. 使用Seurat进行下游分析
- 使用outs内filtered_feeature_bc_matrix中的标准输出文件进行下游的seurat分析
```
##### -1. 导入相关的R包#######
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


# 基本的信息储存在meta.data中，同时如果有需要可以先向meta中添加一些不要的信息
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

# 质控
pro <- subset(pro, subset = nFeature_RNA > 1000 & percent.mt < 20)
VlnPlot(
  pro,
  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
  ncol = 3,
  pt.size = 0
)
# 可视化barcode属性中两个特征的相关关系
plot1 <- FeatureScatter(pro, feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2 <- FeatureScatter(pro, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2
```
  - 质控后可通过小提琴图可视化质控结果，也可以使用相关分析查看质控结果是否理想
  - ![](2023-05-17-15-18-31.png)
  - ![](2023-05-17-15-19-26.png)

```
#####2. 归一化(对数)############################
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
```
  - 此处挑选了4000个高变基因用于后续的降维
![](2023-05-17-15-25-37.png)
```
#####3.5 去除批次效应###########
#由于该数据集中只有单个样品，因此，无需进行去批次处理
#如果有多个生物学重复合并成为一个样品，或分析不同样本间的差异，则需要考虑去除批次

### 1. 首先对于不同的样本分别运行Seurat标准流程到找高变基因这一步
### 2. 然后进行整合 intergration, 该例子中有testA.seu和testB.seu两个样品数据
##### Integration ----
testAB.anchors <- FindIntegrationAnchors(object.list = list(testA.seu,testB.seu), dims = 1:20)
testAB.integrated <- IntegrateData(anchorset = testAB.anchors, dims = 1:20)
# 该步骤以后多了一个整合后的assay，还有一个RNA的assay存放整合前的数据
### 3. 之后基于后续的标准流程进行分析
DefaultAssay(testAB.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
testAB.integrated <- ScaleData(testAB.integrated, features = rownames(testAB.integrated))
testAB.integrated <- RunPCA(testAB.integrated, npcs = 50, verbose = FALSE)
testAB.integrated <- FindNeighbors(testAB.integrated, dims = 1:30)
testAB.integrated <- FindClusters(testAB.integrated, resolution = 0.5)
testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:30)
testAB.integrated <- RunTSNE(testAB.integrated, dims = 1:30)
### 4. 可以查看去批次后的结果
library(cowplot)
testAB.integrated$patient=str_replace(testAB.integrated$orig.ident,"_.*$","")
p1 <- DimPlot(testAB.integrated, reduction = "tsne", group.by = "patient", pt.size=0.5)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
p2 <- DimPlot(testAB.integrated, reduction = "tsne", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
fig_tsne <- plot_grid(p1, p2, labels = c('patient','ident'),align = "v",ncol = 2)
ggsave(filename = "tsne2.pdf", plot = fig_tsne, device = 'pdf', width = 27, height = 12, units = 'cm')
### 上述方法使用与较小的数据集，细胞数尽量不要超过几万，且整合中的任意一个seurat对象细胞数不能少于200
```
- （SRR7722937数据集）寻找高变基因后，继续进行seurat标准流程
```
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
```
- PCA降维结果，在后续流程中，使用PCA降维后的高变异方向再进行后续的非线性降维，以提升聚类效果
- ![](2023-05-17-15-39-52.png)

```
##### 5.5 使用harmony去除批次，如果数据量较大，不建议使用第一中去批次效应方法，
##### 可以在PCA这一步使用harmony整合，这可以理解成一种新的降维，通过harmony进行去批次
### 1. 在harmony流程中首先将两个数据集直接合并
library(harmony)
testdf=cbind(testA,testB)
### 2. 进行seurat标准流程
test.seu <- CreateSeuratObject(counts = testdf) %>%
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData()
test.seu <- RunPCA(test.seu, npcs = 50, verbose = FALSE)
test.seu@meta.data$patient=str_replace(test.seu$orig.ident,"_.*$","")
### 3. 完成PCA以后进行Harmony整合
test.seu=test.seu %>% RunHarmony("patient", plot_convergence = TRUE)
# 可以查看整合后的结果
> test.seu
An object of class Seurat 
33538 features across 6746 samples within 1 assay 
Active assay: RNA (33538 features)
 2 dimensional reductions calculated: pca, harmony
### 4. 基于harmony的Embedding矩阵进行常规的聚类和降维
test.seu <- test.seu %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
```
- 在使用PCA（和Harmony）后，基于该降维的embeddings矩阵进行进一步的非线性降维
```
##### 6.选取前20个PCs对样品进行非线性降维
# 选择非线性降维所使用的PCA维数
ElbowPlot(pro)
# 前20维度能解释几乎所有的变异，选择20维（基于PCA的Embedding矩阵进行非线性降维）
pro <- RunUMAP(pro, dims = 1:20)
Embeddings(pro,"umap") %>% head()
Embeddings(pro,"umap") %>% as.data.frame() %>%
  ggplot(aes(x=UMAP_1,y=UMAP_2)) + geom_point()
```
![](2023-05-17-15-49-32.png)
![](2023-05-17-15-49-49.png)

```
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
```
![](2023-05-17-15-51-04.png)
- 该数据集中的细胞被分为11类

```
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
```
#### 4. 使用scVelo进行RNA velocity分析
- 首先需要根据cellranger的output文件得到一个spliced and unspliced 矩阵
```
#spliced_unspliced_matrices.sh
##### 0. Constructing spliced and unspliced counts matrices #####
# we need to have a matrix for spliced and unspliced transcripts
# provide a .gtf to mask repeat regions (这个文件从https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=&db=mm39&hgta_group=allTracks&hgta_track=rmsk&hgta_table=rmsk&hgta_regionType=genome&position=&hgta_outputType=gff&hgta_outFileName=hg38_rmsk.gtf网站生成)
repeats="/home/wuhangrui/database/ref/hg38/hg38_rmsk.gtf"
transcriptome="/home/wuhangrui/database/ref/hg38/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
cellranger_output="/home/wuhangrui/R/data/20230516/rawdata/SRR7722937/sample/sample1"

velocyto run10x -m $repeats \
                $cellranger_output \
                $transcriptome
```
- 该代码可以得到一个loom文件
- ![](2023-05-17-16-00-39.png)

```
##### 1.加载数据，并赋值给相应变量########## 
# load sparse matrix
X = io.mmread("/home/wuhangrui/R/data/20230516/analysis_in_R/counts.mtx")
# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv("/home/wuhangrui/R/data/20230516/analysis_in_R/metadata.csv")

# load gene names:
with open("/home/wuhangrui/R/data/20230516/analysis_in_R/gene_names.csv","r") as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("/home/wuhangrui/R/data/20230516/analysis_in_R/pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm["X_pca"] = pca.to_numpy()
adata.obsm["X_umap"] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata,color=["seurat_clusters"], frameon=False, save=True)

# save dataset as anndata format
adata.write("pro.h5ad")

# reload dataset
adata = sc.read_h5ad("pro.h5ad")
```
![](2023-05-17-16-02-05.png)

```
##### 2. 将表达矩阵和可变剪切矩阵整合
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2

adata = sc.read_h5ad('pro.h5ad')

# load loom files for spliced/unspliced matrices for each sample:
# 主要是利用barcode进行整合，因此两个矩阵的相对应的barcodes需要完全相同
ldata1 = scv.read('/home/wuhangrui/R/data/20230516/rawdata/SRR7722937/sample/sample1/velocyto/sample1.loom', cache=True)
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1' for bc in barcodes]
ldata1.obs.index = barcodes

# make variable names unique
ldata1.var_names_make_unique()

# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata1)

# plot umap to check
sc.pl.umap(adata, color='seurat_clusters', frameon=False, legend_loc='on data', title='', save='_seurat_clusters.pdf')

```
![](2023-05-17-16-07-28.png)

```
##### 3.compute RNA velocity #################

scv.pl.proportions(adata)

# pre-process
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

# 计算 velocity
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

# 可视化 velocity fields
scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='embedding.pdf')
scv.pl.velocity_embedding_grid(adata, basis='umap', color='seurat_clusters', save='embedding_grid.pdf', title='', scale=0.25)
```
![](2023-05-17-16-11-42.png)
![](2023-05-17-16-11-49.png)
- 可以看见大致的细胞伪时间轨迹的相图

```
### 选择一个基因看一看
# plot velocity of a selected gene
scv.pl.velocity(adata, var_names=['APOE'], color='seurat_clusters')
```
![](2023-05-17-16-12-15.png)
- 这个基因似乎在这群细胞中没有发生什么动态变化
