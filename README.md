# Need Not to Read Me
ThiS is just a record
一名做湿实验的研究僧的生信学习之路
- 20230516 完成了一次单细胞转录组学的生信分析基础流程，该流程以SRR7722937为例，主要包括
-   cellranger上游处理
-   seurat分析cellranger的标准输出
-   尝试进行RNA velocity 分析
-   流程代码与结果输出
</br>

- 20230518 补充了一次20221006初次进行bulk RNA seq的代码与分析流程，主要包括
-   从fastq文件质控开始至获得表达矩阵的上游分析
-   DEseq2差异分析
-   GO富集分析
-   代码流程与结果输出（一次初步的分析）
-   补充 -- GO富集分析可以通过tb-tools可视化工具实现
-   tb-tools整合了相当多的分析功能，是湿实验人的福音
</br>

- 20230531 完成了一次ChIP-Seq生信分析流程，主要包括
- ChIP-Seq的上游分析，fastq文件的质控，过滤，比对
- 比对文件的排序，过滤，去重，索引
- 查看了chipseq分析的质量
- 进行了少许下游分析，call motif 和 TSS 共定位
