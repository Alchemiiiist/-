library(dplyr)
library(Seurat)
library(ggplot2)
##### 0.创建Seurat对象 ###########
# load the bal dataset
# data.dir Directory containing the matrix.mtx, genes.tsv (or features.tsv) and barcodes.tsv files provided by 10X
# gene.column gene列在Dimnames中的位置 default 2
# unique.features Make feature names unique default TRUE
bal.data = Read10X(data.dir = "./data_dir",
                   gene.column = 2,
                   cell.column = 1,
                   unique.features = TRUE)
# 使用原始数据创建Seurat对象
# min.cell 一个基因最少在多少个细胞中表达，default 3
# min.features 一个细胞中最少表达多少个基因,default 200
bal <- CreateSeuratObject(counts = bal.data, project = "bal-w1",
                          min.cells = 3, min.features = 200)

# 查看barcode 的属性信息
# nCount表示UMI数量，nfeature表示捕获的基因数量
head(bal@meta.data)
# 之后会在metadata内继续添加barcode的其他属性信息
# 在metadata中增加一列barcode属性
bal[["CB"]] <- rownames(bal@meta.data)

##### 1.质控#####################################
#三个质控指标
#线粒体中UMI占总UMI比例--不能过高
# nCount_RNA
# nFeature_RNA
# PercentageFeatureSet 对某些基因表达量的和取百分比
# The [[ operator can add columns to object metadata.
bal[["percent.mt"]] <- PercentageFeatureSet(bal,
                                            #pattern = "^MT-"
                                            features = c("ND1","ND2","COX1",
                                                         "COX2","ATP8","ATP6",
                                                         "COX3","ND3","ND4L",
                                                         "ND4","ND5","ND6","CYTB"))
head(bal@meta.data)

# 使用violinplot{Seurat}可视化QC需要QC的指标
VlnPlot(bal,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3, pt.size = 0.005)

# 使用FeatureScatter{Seurat}可视化barcode属性中两个特征的相关关系
plot1 <- FeatureScatter(bal,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2 <- FeatureScatter(bal,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
plot1+plot2

# 进行质量控制
bal <- subset(bal, subset =  nFeature_RNA > 363 & percent.mt < 10)
VlnPlot(bal,features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size = 0)

#####2.归一化(对数)######################################
# 使用NormalizedData(){Seurat}进行归一化
bal <- NormalizeData(bal,normalization.method = "LogNormalize",scale.factor = 10000)
# 将归一化结果储存在了bal@assays$RNA@data中
head(bal@assays$RNA@data)
dim(bal@assays$RNA@data)
# 目前有19829个基因，6233个细胞

##### 3.找高变基因##################################
#高变基因：在不同（例如2000个）细胞中表达差异大的基因（特征选择）
#使用FindFeatures筛选高变基因
# FindVariableFeatures nfeatures default 2000 找2000个高变基因
# FindVariableFeatures vst 筛选方法：2000个方差最大的数据点（那些6233个细胞种表达方差最大的2000个基因）
# 使用lognormalized的矩阵，不使用原始表达矩阵和scaledata矩阵进行筛选
bal = FindVariableFeatures(bal, selection.method = "vst", nfeatures = 2000)
# 没有更改原始的对数阵的data，筛选结果存放在bal@assays$RNA@var.features
head(bal@assays$RNA@var.features)

# 对高变基因可视化
# 使用VariableFeaturePlot(){Seurat}
plot3<-VariableFeaturePlot(bal)
plot3

##### 4.标准化 ###########################################
# 对上述2000个高变基因进行标准化
# 使所有细胞的平均基因表达为0，细胞间的表达方差为1（进行了一次线性变换）
# 进行PCA前的预处理步骤
# 使用ScaleData(){Seurat}进行标准化
bal <- ScaleData(bal)
# 没有更改原始的对数归一阵的data,标准化结果储存在bal@assays$RNA@scale.data
dim(bal@assays$RNA@scale.data)

##### 5.PCA降维 ##########################################
# 使用RunPCA(){Seurat}线性降维
bal <- RunPCA(bal,features = VariableFeatures(object=bal),npcs=50)
# 降维结果储存在bal@meta.data@reductions$pca中
# PCA降维可视化
DimPlot(bal,reduction = "pca")
DimHeatmap(bal,dims = 1:2,cells = 500,balanced = TRUE)

# 确定选择的降维维数
#jackstraw 方法：随机置换一部分数据（默认为1％），然后重新 PCA，重复此过程。
#将包含较多低 P 值特征的主成分为「重要的」主成分。
#JackStraw法相当于计算每个主成分的p值，根据p值选择显著性的PCs纳入下游分析，
#比较科学，但是当数据量比较大时，计算非常慢
#结果储存在reduction pca 下面的jackstraw中
bal <- JackStraw(bal,num.replicate = 100)
# ScoreJackStraw用于量化主成分的显著性强度，富含低P值基因较多的主成分更有统计学意义
# 分数结果储存在jackstraw中的p.values
bal <- ScoreJackStraw(bal,dims = 1:20)
#使用JackStrawPlot()函数可视化比较每个主成分的 p 值分布和均匀分布（虚线）
JackStrawPlot(bal,dims = 1:20)

##### 6. 非线性降维 #####################################################
# 选择前15 个PC 对样本进行t-SNE 降维
bal <- RunTSNE(bal, dims = 1:15)
Embeddings(bal,"tsne") %>% head()
Embeddings(bal,"tsne") %>% as.data.frame() %>% 
  ggplot(aes(x=tSNE_1,y=tSNE_2)) + geom_point()

# bal <- RunUMAP(bal, dims = 1:15)

##### 7.聚类 ############################################################
bal <- FindNeighbors(bal, dims = 1:15)
# res越高，聚类的族数越多
bal <- FindClusters(bal, resolution = 1.0)
# 聚类信息储存在bal@meta.data$seurat_clusters
table(bal@meta.data$seurat_clusters)
# 每个细胞赋了一个新的身份--cluster信息
head(Idents(bal))
# 聚类可视化
DimPlot(bal, reduction = "tsne",pt.size = 1)+
  DimPlot(bal, reduction = "tsne",pt.size = 1,label = T,repel = T,label.size = 6)
#DimPlot(bal_cluster, reduction = "umap",pt.size = 1)+
  #DimPlot(bal_cluster, reduction = "umap",pt.size = 1,label = T,repel = T,label.size = 6)

##### 8.细胞注释 ###################################################
# 找marker基因
# 寻找差异表达的特征（簇生物标志物）
# findmarkers为所有集群自动执行此过程，也可以测试集群组之间的对比
# min.pct: 在两组（检测组和其他组）细胞中的任何一组中以最小百分比检测到一个特征
# thresh: 一个特征在两组之间差异表达（平均）阈值筛选
# FinAllMarkers(){Seurat}找出每个细胞簇的标记物，与所有剩余的细胞进行比较，只报告阳性细胞
# min.pct 表达比例（在15%以上的细胞有表达的基因），logfc差异倍数（）。检测的cluster和其他cluster
bal.makers <- FindAllMarkers(bal, only.pos = T,min.pct = 0.1, logfc.threshold = 0.1, test.use = "wilcox")
bal.makers <- bal.makers%>%filter(p_val_adj < 0.01)
write.csv(bal.makers,file = "14cluster_log2fc1_padj0.01.csv",quote = F,row.names = F)

# 使用celltypist进行注释
# 生成
gene.num <- as.data.frame(rowSums(bal@assays$RNA@counts > 0 ))
colnames(gene.num) <- "num"
gene.num$gene <- rownames(gene.num)
gene.num <- gene.num %>%filter(num > 5)

rawcountmat <- bal@assays$RNA@counts[gene.num$gene,]
rawcountmat <- t(as.matrix(rawcountmat))
write.csv(rawcountmat,file = "rawcountmat.csv",quote = F)

# 读取注释结果
celltypist.anno<-read.csv("predicted_labels_modify.csv",header = T, row.names = 1)
celltypist.anno$CB <- rownames(celltypist.anno)
celltypist.anno <- celltypist.anno[,c("CB","majority_voting_modify")]
# 使用inner_join{dplyr}将注释结果合并进入metadata
bal@meta.data <- bal@meta.data%>%inner_join(celltypist.anno, by="CB")
head(bal@meta.data)
rownames(bal@meta.data) <- bal@meta.data$CB

DimPlot(bal, reduction = "tsne",group.by = "majority_voting_modify.y",pt.size = 1,label = T,repel = T,label.size = 4)+
  DimPlot(bal, reduction = "tsne",group.by = "seurat_clusters",pt.size = 1,label = T,repel = T,label.size = 6)

##### 9.检查结果 ##############################################################
# 使用marker基因检查
FeaturePlot(bal,reduction = "tsne", features = c("THBS1","GZMB","LTB","C1QA"),pt.size = 0.05,ncol = 2)
FeaturePlot(bal,reduction = "tsne", features = c("CPA3","CPVL","CAMP","SCD"),pt.size = 0.05,ncol = 2)
FeaturePlot(bal,reduction = "tsne", features = c("ISG15","CCL22","NFIB","MS4A1"),pt.size = 0.05,ncol = 2)
FeaturePlot(bal,reduction = "tsne", features = c("LOC720839","GATA2"),pt.size = 0.05,ncol = 2)
# 根据文献进行marker归属
FeaturePlot(bal,reduction = "tsne", features = c("CD3D","PTPRC","NKG7","NKG2D"),pt.size = 0.05,ncol = 2)#T_NK
FeaturePlot(bal,reduction = "tsne", features = c("CLEC9A","CD1A","CD1C","MAMU-DRA"),pt.size = 0.05,ncol = 2) #DC
FeaturePlot(bal,reduction = "tsne", features = c("CD19","MS4A1","CD79A","CD1C"),pt.size = 0.05,ncol = 2) #B

VlnPlot(bal,features = c("NKG7","NKG2D"),pt.size = 0,ncol = 1)
VlnPlot(bal,features = c("CLEC9A","CD1A","CD1C","MAMU-DRA"),pt.size = 0,ncol = 1)
VlnPlot(bal,features = c("ISG15","CCL22","NFIB","MS4A1"),pt.size = 0,ncol = 1)


##### 10.整合前面的结果确定细胞类型
bal@meta.data$celltype <- as.numeric(as.character(bal@meta.data$seurat_clusters))
bal@meta.data$celltype[bal@meta.data$celltype %in% c(0,1,2,3,7,9,10,13)] =  "Alveolar_Macrophages"
bal@meta.data$celltype[bal@meta.data$celltype == 11] =  "DC"
bal@meta.data$celltype[bal@meta.data$celltype == 12] =  "memory_B_cells"
bal@meta.data$celltype[bal@meta.data$celltype == 5] =  "mast_cells"
bal@meta.data$celltype[bal@meta.data$celltype %in% c(6,8)] = "T_NK_cells"
bal@meta.data$celltype[bal@meta.data$celltype == 4] =  "ILC_T_cells"
head(bal@meta.data)

DimPlot(bal, reduction = "tsne",group.by = "celltype",pt.size = 1,label = T,repel = T,label.size = 4)

### 11.保存rds #################################################################
saveRDS(bal,file = "bal.1122.rds")

##### 12.可视化结果#####################################################
library(RColorBrewer)
library(scales)
#### 降维图
#color_celltype=c("#bdb9da","#b3de69","#1b9e77","#fa8071","#377db8",)
color_celltype <- c("#1e5670","#4198b9","#6bb3c0","#91cfc9","#cde8f3","#377db8")
# show_col(color_celltype)
names(color_celltype) <- sort(unique(bal@meta.data$celltype))


DimPlot(bal, reduction = "tsne",group.by = "celltype",pt.size = 0.5)+
  scale_color_manual(values = color_celltype)+
  theme_bw()+ #这个主题默认会在图中加几条横竖的灰线
  theme(
    panel.grid = element_blank(), #element_blank()这个函数表示把灰线都去掉
    plot.title = element_text(hjust = 0.5,size = 20), #element_text()设置文本；hjust=0.5表示水平方向处在0.5的位置，也就是中间位置
    axis.text = element_text(size = 14), #axis.text表示坐标轴的文本
    axis.title = element_text(size = 18), #axis.title表示坐标轴的标题
    legend.text = element_text(size = 13) #legend.text表示图例的文本
  )+ guides(colour = guide_legend(override.aes = list(size=5))) #用于修改图例中点的大小
ggsave("tsne_celltype.pdf",width = 19,height = 13.5,units = "cm")

#### 差异基因
bal@meta.data$celltype <- factor(
  bal@meta.data$celltype,
  levels = sort(unique(bal@meta.data$celltype))
) #因子型变量可以控制元素的顺序，这个顺序会体现在画图中
Idents(bal)="celltype" #FindAllMarkers这一步会对identity来操作，所以要提前设置想去比较的那个变量
##### 1. 热图
#每一类选前5个基因，avg_log2FC是权重依据
top5marker=bal.makers %>% group_by(cluster) %>% top_n(wt = avg_log2FC,n = 5) 
DoHeatmap(bal,features = top5marker$gene,size = 5,group.bar.height = 0.05,angle = 30) 

##### 2. 气泡图
DotPlot(bal, features = top5marker$gene)+RotatedAxis()+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  theme_bw()+
  theme(
    axis.title = element_blank(), #element_blank()表示清空坐标轴title
    axis.text.x.bottom = element_text(hjust = 1, angle = 45,color = "black"), #angle表示x轴的文本旋转一定角度
    axis.text.y.left = element_text(color = "black"),
    
    legend.position = "top" #图例放在上方
  )
ggsave("marker_bubble.pdf",width = 20,height = 9,units = "cm")


##### 3. 基因表达梯度图
FeaturePlot(bal,features = "CD3D",reduction = "tsne",pt.size = 1)+
  scale_color_gradient(low = "lightblue",high = "red")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    
    plot.title = element_text(hjust = 0.5,size=14)
  )
FeaturePlot(bal,features = "MAMU-DRA",reduction = "tsne",pt.size = 1)+
  scale_color_gradient(low = "lightblue",high = "red")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    
    plot.title = element_text(hjust = 0.5,size=14)
  )


##### 4. 小提琴图
VlnPlot(bal,features = c("CD3D","PTPRC"),
        group.by = "celltype",pt.size = 0,ncol = 1)&
  scale_y_continuous(expand = c(0.02,0))& #expand表示延伸，在y轴方向，图的主体部分除了画到y的最小值和最大值，还向上下两侧延伸一点点
  theme(
    axis.title = element_blank(),
    axis.text.x.bottom = element_text(size = 14),
    axis.text.y.left = element_text(size = 16),
    axis.ticks.length = unit(0.2,"cm"), #坐标轴刻度线的长度
    
    plot.title = element_text(size = 20,hjust = 0.5)
  )

##### 5. 各个样本中不同细胞类型的比例 
bar.df=bal@meta.data
text.df=as.data.frame(table(bar.df$orig.ident))

bar.df%>%ggplot(aes(x=orig.ident))+
  geom_bar(aes(fill=celltype))+
  scale_x_discrete("")+
  scale_y_continuous("cell number",expand = c(0.01,0))+
  scale_fill_manual(values = color_celltype)+
  # 在每个bar上面添加样本细胞数
  geom_text(data = text.df,aes(x=Var1,y=Freq+100,label=Freq),size=4)+
  theme_classic()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,color = "black",size = 12),
    axis.text.y.left = element_text(color = "black",size = 12)
  )
