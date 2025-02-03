###The original data is too large to be uploaded. Please contact the corresponding author to obtain it after the article is published.



rm(list = ls())
options(stringsAsFactors = F)
library(doParallel) 
registerDoParallel(cores=detectCores())

#install.packages("Seurat")  #用来做图
#install.packages("outliers")
#install.packages("pbmcapply")
#install.packages("doFuture")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("singscore")
#BiocManager::install("GSVA")
#BiocManager::install("GSEABase")
#BiocManager::install("limma")

#install.packages("devtools")
#library(devtools)
#devtools::install_github('dviraran/SingleR')#此包是用来做细胞注释



library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(CellChat)
library(patchwork)
library(celldex)
library(SingleR)


#BiocManager::install("limma")
#BiocManager::install("celldex")
BiocManager::install("SingleR")


rm(list = ls())
options(stringsAsFactors = F)

###################################04.数据前期处理和矫正###################################
#读取数据
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(CellChat)
library(patchwork)

data_treat1=read.csv('data_treat1.csv',row.names=1,header=T)  #第一列设置为行名，第一行设置为列名
data_treat2=read.csv('data_treat2.csv',row.names=1,header=T)
data_treat3=read.csv('data_treat3.csv',row.names=1,header=T)
data_control1=read.csv('data_control1.csv',row.names=1,header=T)
data_control2=read.csv('data_control2.csv',row.names=1,header=T)
data_control3=read.csv('data_control3.csv',row.names=1,header=T)



#查看前6行、6列，即查看数据的格式
data_treat1[1:6,1:6 ]
data_treat2[1:6,1:6 ]
data_treat3[1:6,1:6 ]
data_control1[1:6,1:6 ]
data_control2[1:6,1:6 ]
data_control3[1:6,1:6 ]



#将矩阵转换为Seurat对象，并对数据进行过滤。某个基因少于3个细胞表达、某个细胞测得的基因数少于200的删掉。names.delim为样品名的分隔符
#min.cells = 3 意味着只包括至少有3个细胞表达的基因，这有助于过滤掉稀有的或低表达的基因，以简化分析。
#min.features = 200 意味着只包括至少在200个细胞中具有表达的基因，这有助于过滤掉稀有的或低表达的基因，以简化分析。
#names.delim = "_"    。names.delim为样品名的分隔符
t1 = CreateSeuratObject(counts = data_treat1, min.cells = 500, min.features = 500) 
head(t1@meta.data)
t2 = CreateSeuratObject(counts = data_treat2, min.cells = 500, min.features = 500) 
head(t2@meta.data)
t3 = CreateSeuratObject(counts = data_treat3, min.cells = 500, min.features = 500) 
head(t3@meta.data)

c1 = CreateSeuratObject(counts = data_control1, min.cells = 500, min.features = 500) 
head(c1@meta.data)
c2 = CreateSeuratObject(counts = data_control2, min.cells = 500, min.features = 500) 
head(c2@meta.data)
c3 = CreateSeuratObject(counts = data_control3, min.cells = 500, min.features = 500) 
head(c3@meta.data)



#查看行、列数
dim(data_treat1)
dim(data_treat2)
dim(data_treat3)
dim(data_control1)
dim(data_control2)
dim(data_control3)
#发现c1基因数最多，而c3细胞数最多

dim(t1@meta.data)
dim(t2@meta.data)
dim(t3@meta.data)
dim(c1@meta.data)
dim(c2@meta.data)
dim(c3@meta.data)
#发现c3行数最多，即细胞数目最多


             #######################以下几行代码不用运行#######################
#将以上几个Seurat合并成一个
m = merge(x=c1,y=c(c2,c3,t1,t2,t3), merge.data = TRUE)
#n = merge(x=c3,y=c(c1,c2,t1,t2,t3),merge.data = TRUE)

dim(m@meta.data)
#dim(n@meta.data)
#发现两个的行数，即细胞数是一样的

head(m@meta.data)
#head(n@meta.data)
#查看列名即细胞，与原数据一致，即c1为_1结尾，c3为_6结尾。

dim(m@assays$RNA@counts)
#dim(n@assays$RNA@counts)
#查看合并后表达矩阵，发现基因数、细胞数是一致的，所以在上面合并的代码中，x无论是哪个样本的数据都没有影响

list(m)
pbmc=m
head(m@assays$RNA@counts[1:10,1:10])
                    #######################以上几行代码不用运行#######################


pbmc = merge(x=c1,y=c(c2,c3,t1,t2,t3), merge.data = TRUE)

#将矩阵转换为Seurat对象，并对数据进行过滤。某个基因少于3个细胞表达、某个细胞测得的基因数少于50的删掉。names.delim为样品名的分隔符
#pbmc <- CreateSeuratObject(counts = data,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_",)

#使用PercentageFeatureSet函数计算线粒体基因的百分比。PercentageFeatureSet: 这是Seurat包中的函数，用于计算单细胞数据中特定基因集的百分比。
#  在这里，它将计算线粒体基因的百分比。object = pbmc: 这是PercentageFeatureSet函数的参数，指定要计算百分比的Seurat对象，即pbmc。
#  pattern = "^MT-": 这是PercentageFeatureSet函数的参数，指定要计算百分比的基因模式。在这里，它使用正则表达式模式 "^MT-" 匹配所有以 "MT-" 开头的基因，通常表示线粒体基因。
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")


#save(pbmc,file="pbmc  cells500   features500.Rdata") 
#load("pbmc  cells500   features500.Rdata")

#画小提琴图
pdf(file="04.featureViolin.pdf",width=10,height=6)           #保存基因特征小提琴图
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf(file="04.featureViolin  无细胞点.pdf",width=10,height=6)           #保存基因特征小提琴图
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()


pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 &  percent.mt < 10)  
#save(pbmc,file="pbmc  cells500   features大于500小于10000  percent.mt小于10.Rdata") 
load("pbmc  cells500   features大于500小于10000  percent.mt小于10.Rdata")
#对pbmc进一步质控过滤
#  nFeature_RNA在每个细胞中检测到的独特基因的数量。低质量的细胞或空液滴通常只有很少的基因，细胞双联体或多重细胞可能表现出异常高的基因计数。总结：过高和过低的nFeature_RNA要滤过
#  nCount_RNA在细胞内检测到的分子总数（与独特基因密切相关）。一般文献里不对nCount_RNA进行过滤。总结：过低的nCount_RNA要滤过。
#  percent.mt映射到线粒体基因组的读段百分比，低质量/垂死的细胞通常表现出广泛的线粒体污染。


#以下三行代码为后续细胞互作、细胞通讯分析做准备
#data.input总的文件=as.matrix(pbmc@assays$RNA@data)
#data.input总的文件[1:6,1:6 ]
#save(data.input总的文件,file="data.input总的文件.Rdata")
##load("data.input总的文件.Rdata")


pdf(file="04.featureViolin  无细胞点  线粒体小于10  基因小于1万大于五百.pdf",width=10,height=6)           #保存基因特征小提琴图
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()


#pbmc1 <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 &  percent.mt > 10)  
#pdf(file="04.featureViolin  无细胞点  线粒体大于10  基因小于1万大于五百.pdf",width=10,height=6)           #保存基因特征小提琴图
#VlnPlot(object = pbmc1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
#dev.off()


#测序深度的相关性绘图
pdf(file="04.featureCor.pdf",width=10,height=6)              #保存基因特征相关性图
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",,pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

#对数据进行标准化
#pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#提取那些在细胞间变异系数较大的基因
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)   
#输出特征方差图
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="04.featureVar.pdf",width=10,height=6)              #保存基因特征方差图
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)  #将前10个基因进行标注
CombinePlots(plots = list(plot1, plot2))
dev.off()


#save(pbmc,file="pbmc  cells500   features500  percent.mt小于10  标准化之后.Rdata") 
#load("pbmc  cells500   features500  percent.mt小于10  标准化之后.Rdata")


###################################05.PCA主成分分析###################################
##PCA分析
pbmc=ScaleData(pbmc)                     #PCA降维之前的标准预处理步骤
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = bpmc))     #PCA分析。#npcs = 20为对前20维基因进行分析。#VariableFeatures(object = bpmc)为前面分析出的，波动最大的基因进行组成分分析

#绘制每个PCA成分的相关基因
pdf(file="05.pcaGene.pdf",width=10,height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)  # dims = 1:4把前4个PC的图进行绘制。  #reduction = "pca" 降维方法为pca。  #nfeatures = 20  把前面20个基因呈现出来。
dev.off()

#主成分分析图形
pdf(file="05.PCA.pdf",width=6.5,height=6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()

#主成分分析热图
pdf(file="05.pcaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)#dims = 1:4为画前4个PC，nfeatures = 30为画前30个基因，ncol=2为每一行放2个热图
dev.off()

#每个PC的p值分布和均匀分布
pbmc <- JackStraw(object = pbmc, num.replicate = 100) #这步耗时较长，约23分钟
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
pdf(file="05.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object = pbmc, dims = 1:20)
dev.off()



###################################06.TSNE聚类分析和marker基因###################################
###TSNE聚类和UMAP聚类只选其中一种就可以了


##TSNE聚类分析
pcSelect=20    #要看有多少个PC的p值是小于0.05的就写多少，这里是20个都有意义
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                #计算邻接距离    #此步耗时较长
pbmc <- FindClusters(object = pbmc, resolution = 0.5)                  #对细胞分组,优化标准模块化
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)                      #TSNE聚类    #此步耗时较长
DimPlot(pbmc, reduction = "tsne")
pdf(file="06.TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = pbmc, label = TRUE, pt.size = 0.5)    #TSNE可视化#pt.size = 0.5为图形中点的大小  
dev.off()
write.table(pbmc$seurat_clusters,file="06.tsneCluster.txt",quote=F,sep="\t",col.names=F)



##UMAP聚类分析
#pcSelect=20    #要看有多少个PC的p值是小于0.05的就写多少，这里是20个都有意义
#pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)                #计算邻接距离    
#pbmc <- FindClusters(object = pbmc, resolution = 0.5)                  #对细胞分组,优化标准模块化
pbmc <- RunUMAP(object = pbmc, dims = 1:pcSelect)                      #UMAP聚类
DimPlot(pbmc, reduction = "umap")
pdf(file="06.UMAP.pdf",width=6.5,height=6)
UMAPPlot(object = pbmc, label = TRUE, pt.size = 0.5)    #TSNE可视化#pt.size = 2为图形中点的大小  #不懂这步为啥画不了图。直接从上一步绘制的图导出PDF文件
dev.off()
write.table(pbmc$seurat_clusters,file="06.umapCluster.txt",quote=F,sep="\t",col.names=F)


##寻找差异表达的特征
logFCfilter=2    #绝对值大于1
adjPvalFilter=0.05  #校正的p值小于0.05
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,       #定义最小的PCT值
                               logfc.threshold = logFCfilter)          #此步耗时较长
#save(pbmc,file="pbmc  cells500   features500  筛选差异基因后.Rdata") 
#load("pbmc  cells500   features500  筛选差异基因后.Rdata")

#   FindAllMarkers 是Seurat包中的一个函数，用于查找在不同细胞群中具有差异表达的基因。
#   object 是包含单细胞RNA测序数据的Seurat对象，其中包含了细胞表达矩阵、细胞群信息和其他相关信息。
#   only.pos 是一个逻辑参数，控制是否只查找在某一细胞群中具有正向差异表达的基因。
#            如果设置为 FALSE，则会查找在所有细胞群中具有差异表达的基因。
#   min.pct 是最小百分比阈值，用于指定一个基因在某个细胞群中的表达比例必须大于等于该阈值才会被考虑。
#   logfc.threshold 是对数折叠变化（log fold change）的阈值，用于筛选差异表达的基因。只有具有大于等于该阈值的对数折叠变化的基因才会被保留。

write.table(pbmc.markers,file="06.pbmc.markers cells500   features500  未进行第二次标准化  logFC绝对值大于等于2  adjp  0.05.xls",sep="\t",row.names=F,quote=F)
pbmc.markers=read.table("06.pbmc.markers cells500   features500  未进行第二次标准化  logFC绝对值大于等于2  adjp  0.05.xls",header=T,sep="\t")
dim(pbmc.markers)  #查看有多少
length(unique(pbmc.markers$gene))   #查看有多少不同名的gene symbol

sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]   #进一步过滤
dim(sig.markers)      #查看有多少不同名的gene symbol
length(unique(sig.markers$gene)) #查看有多少不同名的gene symbol，再根据基因数目进行调整上面logFCfilter的值
write.table(sig.markers,file="06.sig.markers cells500   features500  未进行第二次标准化   logFC绝对值大于等于2  adjP小于0.05.xls",sep="\t",row.names=F,quote=F)
#sig.markers=read.table ("06.sig.markers cells500   features500  未进行第二次标准化   logFC绝对值大于等于2  adjP小于0.05.xls",header=T,sep="\t")

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#绘制marker在各个cluster的热图
pdf(file="06.tsneHeatmap  各个cluster前10个基因热图.pdf",width=40,height=20)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()   #使用了前十个基因进行绘制热图  ##有很多基因不存在，不懂为啥
dev.off()

top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#绘制marker在各个cluster的热图
pdf(file="06.tsneHeatmap  各个cluster前5个基因热图.pdf",width=40,height=20)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()   #使用了前五个基因进行绘制热图  ##有很多基因不存在，不懂为啥
dev.off()


#绘制marker的小提琴图
pdf(file="06.markerViolin.pdf", width=10,height=6)
VlnPlot(object = pbmc, features = c("AREG", "AZU1"))  #"HLA−DRA", "HLA−DQA1","HLA−DQB1"，定义基因，关注什么基因就把什么基因绘制出来。绘制的基因在列表中不能有相同的，否则不能绘制
dev.off()

#绘制marker在各个cluster的散点图
pdf(file="06.markerScatter.pdf",width=10,height=6)
FeaturePlot(object = pbmc, features = c("AREG", "AZU1"),cols = c("green", "red"))
dev.off()

#绘制marker在各个cluster的气泡图
pdf(file="06.markerBubble.pdf",width=15,height=6)
cluster10Marker=c("MZB1","HMOX1","TIMP1","AZU1","FABP5","HLA-C","SUB1","H3F3A","PPIB","NAMPT",
                  "ACTG1","RPS4X","MT-ND1","RPL3","S100A6","RPS3A","HIST1H4C","TPT1","RPL28","NEAT1")  #定义基因。这里是找出logFC最大的5个基因和最小的5个基因进行绘制。
DotPlot(object = pbmc, features = cluster10Marker)
dev.off()

#拿每个cluster前5个基因记性绘制气泡图
pdf(file="06.markerBubble  拿每个cluster前5个基因记性绘制气泡图 top5.pdf",width=15,height=6)
cluster10Marker=top5  #定义基因。这里是找出logFC最大的5个基因和最小的5个基因进行绘制。
unique_levels <- unique(cluster10Marker)

features_to_check <- cluster10Marker
present_features <- features_to_check[features_to_check %in% colnames(pbmc@assays$RNA@data)]

if (length(present_features) > 0) {
  cat("在数据集中存在的变量：", paste(present_features, collapse = ", "), "\n")
} else {
  cat("在数据集中找不到指定的任何变量。\n")
}




DotPlot(object = pbmc, features = cluster10Marker)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



#拿每个cluster前5个基因记性绘制气泡图
pdf(file="06.markerBubble  拿每个cluster前5个基因记性绘制气泡图.pdf",width=15,height=6)
cluster10Marker=c("LINC01138","PITPNA","SELENOW","RPF1","RWDD3",
                  "CCL4","XPO7","PAM16","MYBL1","YES1",
                  "IL32","CNOT6L","TOX","SH2D1A","CD3D",
                  "SYNE1","REST","KLRC2","PPP2R2B","SPON2",
                  "PLCL1","CEMIP2","RCAN3","SLC2A3","RHOH",
                  "RANBP9","PHKB","ATG13","COX19","DOP1B", #5
                  "SPEN","STX5","ADAM17","SAE1","LPCAT2", #6
                  "FABP5","CD63","AKAP13","CTSD","LGALS3" )  #定义基因。这里是找出logFC最大的5个基因和最小的5个基因进行绘制。
DotPlot(object = pbmc, features = cluster10Marker)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()





#save(pbmc,file="0.6步结束后工作区所有文件  cells500   features500  未进行第二次标准化.Rdata")    #保存目前的工作区所有文件
#load("0.6步结束后工作区所有文件  cells500   features500  未进行第二次标准化.Rdata")





###################################07.注释细胞类型###################################
#BiocManager::install ("SingleR")
#BiocManager::install("celldex")

#教程视频https://www.bilibili.com/video/BV1xF411f7jD/?spm_id_from=333.337.search-card.all.click&vd_source=09764cf484157501fb04568494820124
library(celldex)
library(SingleR)

##载入人类参考数据集
#BiocManager::install("ExperimentHub")
#   library(ExperimentHub)
#  ref_hs=HumanPrimaryCellAtlasData()   #因获取较难，因此获取后保存到本地
#   save(ref_hs,file="ref_hs.hs.Rdata")
load("ref_hs.hs.Rdata")
##载入小鼠参考数据集
#ref_mm <- MouseRNAseqData()
##singleR有7大数据集

###把rna的转录表达数据提取
testdata <- GetAssayData(pbmc, slot="data")   #提取表达矩阵，单细胞的表达矩阵
testdata[1:6,1:6]
clusters <- pbmc$seurat_clusters              #提取分组信息（前面通过第6步聚类分析得到的34组分组）  #注意第6步是用什么方法进行聚类
cellpred <- SingleR(test = testdata,   #表达矩阵
                    ref = ref_hs,   #参考矩阵。
                    labels = ref_hs$label.main,   #参考的细胞名称
                    clusters = clusters,assay.type.test = "logcounts",
                    assay.type.ref = "logcounts")
##添加到metadata当中
celltype = data.frame(ClusterID=rownames(cellpred), 
                      celltype=cellpred$labels, stringsAsFactors = FALSE)
pbmc@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  pbmc@meta.data[which(pbmc@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

levels(clusters)   #查看原始分得23类
unique(celltype$celltype)   #查看有多少种分类


####注意第6步是用什么方法进行聚类，是tsne聚类还是umap聚类，要对应好只运行下面其中一种即可

#tsne聚类
DimPlot(pbmc, reduction = "tsne",label = T)
DimPlot(pbmc, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（AS、HC全部）.pdf",width=8,height=6)
DimPlot(pbmc, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（AS、HC全部）.pdf",width=8,height=6)
DimPlot(pbmc, reduction = "tsne", group.by = "celltype",label = T)
dev.off()



library(ggsci)
library(ggplot2)
#差异比例箱线图
#比例图
#library(stringr)
#phe=str_split(rownames(pbmc@meta.data),'_',simplify = T)
#head(phe)
#table(phe[,2])
#table(phe[,3])
#pbmc$group4=phe[,2]
tb <- data.frame(table(pbmc@meta.data$celltype,pbmc@meta.data$orig.ident))
tb$Var3=tb$Var2
#tb$Var3=gsub("[-,C,E,-,1,2,3,4,5,6,7,8,9,0]", "", tb$Var3)
#tb$Total <- apply(tb,1,function(x)sum(tb[tb$Var1 == x[1],3]))
tb$Total <- apply(tb,1,function(x)sum(tb[tb$Var2 == x[2],3]))
tb<- tb %>% mutate(Percentage = round(Freq/Total,3) * 100)
table(tb$Var3,tb$Var1)

tb=tb[,c(1,4,6)]
tb$Var1=as.factor(tb$Var1)
tb$Var3=as.factor(tb$Var3)

ggplot(tb) +  
  geom_bar(aes(x =Percentage, y=Var3 , fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+  
  theme_classic() + 
  labs(x='Ratio',y = 'Sample')+ 
  coord_flip()+ 
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave(filename='percent_celltype.pdf',    
       height = 5,width = 5)


DimPlot(pbmc, reduction = "tsne", group.by = "celltype",    
        split.by = 'orig.ident',    
        label = T,pt.size = 0.1,label.size = 3,    
        repel = T,label.box = T) +  
  scale_colour_manual(values = pal_d3("category20")(20),           
                      aesthetics = c("colour", "fill"))
ggsave('umap-by-celltype-ggsci.pdf',height = 4,width=8)


















# 提取只含AS的pbmc
pbmc_AS <- pbmc[, grep("^AS", colnames(pbmc@assays$RNA), value = TRUE)]
# 提取只含HC的pbmc
pbmc_HC <- pbmc[, grep("^HC", colnames(pbmc@assays$RNA), value = TRUE)]

#save(pbmc_AS,file="0.7步后pbmc_AS  cells500   features500  未进行第二次标准化.Rdata")    #保存
#load("0.7步后pbmc_AS  cells500   features500  未进行第二次标准化.Rdata")
#save(pbmc_HC,file="0.7步后pbmc_HC  cells500   features500  未进行第二次标准化.Rdata")    #保存
#load("0.7步后pbmc_HC  cells500   features500  未进行第二次标准化.Rdata")



#AS的tsne聚类
DimPlot(pbmc_AS, reduction = "tsne",label = T)
DimPlot(pbmc_AS, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（AS）.pdf",width=8,height=6)
DimPlot(pbmc_AS, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（AS）.pdf",width=8,height=6)
DimPlot(pbmc_AS, reduction = "tsne", group.by = "celltype",label = T)
dev.off()

#HC的tsne聚类
DimPlot(pbmc_HC, reduction = "tsne",label = T)
DimPlot(pbmc_HC, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（HC）.pdf",width=8,height=6)
DimPlot(pbmc_HC, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（HC）.pdf",width=8,height=6)
DimPlot(pbmc_HC, reduction = "tsne", group.by = "celltype",label = T)
dev.off()



#读取细胞注释信息
cellAnn=read.table("07.all cellAnn cells500   features500  未进行第二次标准化.txt",sep="\t",header=F)
unique(cellAnn$V8)
#获取细胞注释后的细胞类型名称
# 分别提取各种细胞的所有细胞
T_cells <- cellAnn$V1[cellAnn$V8 == "T_cells"]

# 提取T_cells的pbmc
pbmc_T_cells<- pbmc[, grep("^AS", colnames(pbmc@assays$RNA), value = TRUE)]
pbmc_T_cells <- subset(pbmc, pbmc@meta.data$celltype== T_cells)
pbmc_T_cells <- pbmc[, grep(pbmc@meta.data$celltype== T_cells, value = TRUE)]


############### 提取T_cells的pbmc############
T_cells <- "T_cells"
# 使用 Seurat 的子集函数来提取特定细胞类型的子集
pbmc_T_cells <- subset(pbmc, subset = celltype == T_cells)

# 在AS中提取T_cells的pbmc
pbmc_AS_T_cells <- pbmc_T_cells[, grep("^AS", colnames(pbmc_T_cells@assays$RNA), value = TRUE)]
# 在HC中提取T_cells的pbmc
pbmc_HC_T_cells <- pbmc_T_cells[, grep("^HC", colnames(pbmc_T_cells@assays$RNA), value = TRUE)]

#AS的tsne聚类
DimPlot(pbmc_AS_T_cells, reduction = "tsne",label = T)
DimPlot(pbmc_AS_T_cells, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（AS_T_cells）.pdf",width=8,height=6)
DimPlot(pbmc_AS_T_cells, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（AS_T_cells）.pdf",width=8,height=6)
DimPlot(pbmc_AS_T_cells, reduction = "tsne", group.by = "celltype",label = T)
dev.off()

#HC的tsne聚类
DimPlot(pbmc_HC_T_cells, reduction = "tsne",label = T)
DimPlot(pbmc_HC_T_cells, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（HC_T_cells）.pdf",width=8,height=6)
DimPlot(pbmc_HC_T_cells, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（HC_T_cells）.pdf",width=8,height=6)
DimPlot(pbmc_HC_T_cells, reduction = "tsne", group.by = "celltype",label = T)
dev.off()




############### 提取B_cell的pbmc############
B_cell <- "B_cell"
# 使用 Seurat 的子集函数来提取特定细胞类型的子集
pbmc_B_cell <- subset(pbmc, subset = celltype == B_cell)

# 在AS中提取B_cell的pbmc
pbmc_AS_B_cell <- pbmc_B_cell[, grep("^AS", colnames(pbmc_B_cell@assays$RNA), value = TRUE)]
# 在HC中提取B_cell的pbmc
pbmc_HC_B_cell <- pbmc_B_cell[, grep("^HC", colnames(pbmc_B_cell@assays$RNA), value = TRUE)]

#AS的tsne聚类
DimPlot(pbmc_AS_B_cell, reduction = "tsne",label = T)
DimPlot(pbmc_AS_B_cell, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（AS_B_cell）.pdf",width=8,height=6)
DimPlot(pbmc_AS_B_cell, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（AS_B_cell）.pdf",width=8,height=6)
DimPlot(pbmc_AS_B_cell, reduction = "tsne", group.by = "celltype",label = T)
dev.off()

#HC的tsne聚类
DimPlot(pbmc_HC_B_cell, reduction = "tsne",label = T)
DimPlot(pbmc_HC_B_cell, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（HC_B_cell）.pdf",width=8,height=6)
DimPlot(pbmc_HC_B_cell, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（HC_B_cell）.pdf",width=8,height=6)
DimPlot(pbmc_HC_B_cell, reduction = "tsne", group.by = "celltype",label = T)
dev.off()



############### 提取NK_cell的pbmc############
NK_cell <- "NK_cell"
# 使用 Seurat 的子集函数来提取特定细胞类型的子集
pbmc_NK_cell <- subset(pbmc, subset = celltype == NK_cell)

# 在AS中提取NK_cell的pbmc
pbmc_AS_NK_cell <- pbmc_NK_cell[, grep("^AS", colnames(pbmc_NK_cell@assays$RNA), value = TRUE)]
# 在HC中提取NK_cell的pbmc
pbmc_HC_NK_cell <- pbmc_NK_cell[, grep("^HC", colnames(pbmc_NK_cell@assays$RNA), value = TRUE)]

#AS的tsne聚类
DimPlot(pbmc_AS_NK_cell, reduction = "tsne",label = T)
DimPlot(pbmc_AS_NK_cell, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（AS_NK_cell）.pdf",width=8,height=6)
DimPlot(pbmc_AS_NK_cell, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（AS_NK_cell）.pdf",width=8,height=6)
DimPlot(pbmc_AS_NK_cell, reduction = "tsne", group.by = "celltype",label = T)
dev.off()

#HC的tsne聚类
DimPlot(pbmc_HC_NK_cell, reduction = "tsne",label = T)
DimPlot(pbmc_HC_NK_cell, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（HC_NK_cell）.pdf",width=8,height=6)
DimPlot(pbmc_HC_NK_cell, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（HC_NK_cell）.pdf",width=8,height=6)
DimPlot(pbmc_HC_NK_cell, reduction = "tsne", group.by = "celltype",label = T)
dev.off()



############### 提取Monocyte的pbmc############
Monocyte <- "Monocyte"
# 使用 Seurat 的子集函数来提取特定细胞类型的子集
pbmc_Monocyte <- subset(pbmc, subset = celltype == Monocyte)

# 在AS中提取Monocyte的pbmc
pbmc_AS_Monocyte <- pbmc_Monocyte[, grep("^AS", colnames(pbmc_Monocyte@assays$RNA), value = TRUE)]
# 在HC中提取Monocyte的pbmc
pbmc_HC_Monocyte <- pbmc_Monocyte[, grep("^HC", colnames(pbmc_Monocyte@assays$RNA), value = TRUE)]

#AS的tsne聚类
DimPlot(pbmc_AS_Monocyte, reduction = "tsne",label = T)
DimPlot(pbmc_AS_Monocyte, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（AS_Monocyte）.pdf",width=8,height=6)
DimPlot(pbmc_AS_Monocyte, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（AS_Monocyte）.pdf",width=8,height=6)
DimPlot(pbmc_AS_Monocyte, reduction = "tsne", group.by = "celltype",label = T)
dev.off()

#HC的tsne聚类
DimPlot(pbmc_HC_Monocyte, reduction = "tsne",label = T)
DimPlot(pbmc_HC_Monocyte, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（HC_Monocyte）.pdf",width=8,height=6)
DimPlot(pbmc_HC_Monocyte, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（HC_Monocyte）.pdf",width=8,height=6)
DimPlot(pbmc_HC_Monocyte, reduction = "tsne", group.by = "celltype",label = T)
dev.off()



############### 提取GMP的pbmc############
GMP <- "GMP"
# 使用 Seurat 的子集函数来提取特定细胞类型的子集
pbmc_GMP <- subset(pbmc, subset = celltype == GMP)

# 在AS中提取GMP的pbmc
pbmc_AS_GMP <- pbmc_GMP[, grep("^AS", colnames(pbmc_GMP@assays$RNA), value = TRUE)]
# 在HC中提取GMP的pbmc
pbmc_HC_GMP <- pbmc_GMP[, grep("^HC", colnames(pbmc_GMP@assays$RNA), value = TRUE)]

#AS的tsne聚类
DimPlot(pbmc_AS_GMP, reduction = "tsne",label = T)
DimPlot(pbmc_AS_GMP, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（AS_GMP）.pdf",width=8,height=6)
DimPlot(pbmc_AS_GMP, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（AS_GMP）.pdf",width=8,height=6)
DimPlot(pbmc_AS_GMP, reduction = "tsne", group.by = "celltype",label = T)
dev.off()

#HC的tsne聚类
DimPlot(pbmc_HC_GMP, reduction = "tsne",label = T)
DimPlot(pbmc_HC_GMP, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（HC_GMP）.pdf",width=8,height=6)
DimPlot(pbmc_HC_GMP, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（HC_GMP）.pdf",width=8,height=6)
DimPlot(pbmc_HC_GMP, reduction = "tsne", group.by = "celltype",label = T)
dev.off()



############### 提取PRO-B_CELL_CD34+的pbmc############
PRO_B_CELL_CD34 <- "Pro-B_cell_CD34+"
# 使用 Seurat 的子集函数来提取特定细胞类型的子集
pbmc_PRO_B_CELL_CD34 <- subset(pbmc, subset = celltype == PRO_B_CELL_CD34)

# 在AS中提取PRO-B_CELL_CD34+的pbmc
pbmc_AS_PRO_B_CELL_CD34 <- pbmc_PRO_B_CELL_CD34[, grep("^AS", colnames(pbmc_PRO_B_CELL_CD34@assays$RNA), value = TRUE)]
# 在HC中提取PRO-B_CELL_CD34+的pbmc
pbmc_HC_PRO_B_CELL_CD34 <- pbmc_PRO_B_CELL_CD34[, grep("^HC", colnames(pbmc_PRO_B_CELL_CD34@assays$RNA), value = TRUE)]

#AS的tsne聚类
DimPlot(pbmc_AS_PRO_B_CELL_CD34, reduction = "tsne",label = T)
DimPlot(pbmc_AS_PRO_B_CELL_CD34, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（AS_PRO-B_CELL_CD34+）.pdf",width=8,height=6)
DimPlot(pbmc_AS_PRO_B_CELL_CD34, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（AS_PRO-B_CELL_CD34+）.pdf",width=8,height=6)
DimPlot(pbmc_AS_PRO_B_CELL_CD34, reduction = "tsne", group.by = "celltype",label = T)
dev.off()

#HC的tsne聚类
DimPlot(pbmc_HC_PRO_B_CELL_CD34, reduction = "tsne",label = T)
DimPlot(pbmc_HC_PRO_B_CELL_CD34, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（HC_PRO-B_CELL_CD34+）.pdf",width=8,height=6)
DimPlot(pbmc_HC_PRO_B_CELL_CD34, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（HC_PRO-B_CELL_CD34+）.pdf",width=8,height=6)
DimPlot(pbmc_HC_PRO_B_CELL_CD34, reduction = "tsne", group.by = "celltype",label = T)
dev.off()



############### 提取MEP的pbmc############
MEP <- "MEP"
# 使用 Seurat 的子集函数来提取特定细胞类型的子集
pbmc_MEP <- subset(pbmc, subset = celltype == MEP)

# 在AS中提取MEP的pbmc
pbmc_AS_MEP <- pbmc_MEP[, grep("^AS", colnames(pbmc_MEP@assays$RNA), value = TRUE)]
# 在HC中提取MEP的pbmc
pbmc_HC_MEP <- pbmc_MEP[, grep("^HC", colnames(pbmc_MEP@assays$RNA), value = TRUE)]   #MEP细胞在HC组中不存在

#AS的tsne聚类
DimPlot(pbmc_AS_MEP, reduction = "tsne",label = T)
DimPlot(pbmc_AS_MEP, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（AS_MEP）.pdf",width=8,height=6)
DimPlot(pbmc_AS_MEP, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（AS_MEP）.pdf",width=8,height=6)
DimPlot(pbmc_AS_MEP, reduction = "tsne", group.by = "celltype",label = T)
dev.off()

#HC的tsne聚类
DimPlot(pbmc_HC_MEP, reduction = "tsne",label = T)
DimPlot(pbmc_HC_MEP, reduction = "tsne", group.by = "celltype",label = T)
pdf(file="07.tsne细胞分类未注释（HC_MEP）.pdf",width=8,height=6)
DimPlot(pbmc_HC_MEP, reduction = "tsne",label = T)
dev.off()
pdf(file="07.tsne细胞分类细胞注释（HC_MEP）.pdf",width=8,height=6)
DimPlot(pbmc_HC_MEP, reduction = "tsne", group.by = "celltype",label = T)
dev.off()






#umap聚类
DimPlot(pbmc, reduction = "umap",label = T)
DimPlot(pbmc, reduction = "umap", group.by = "celltype",label = T)
pdf(file="07.umap细胞分类未注释.pdf",width=8,height=6)
DimPlot(pbmc, reduction = "umap",label = T)
dev.off()
pdf(file="07.umap细胞分类细胞注释.pdf",width=8,height=6)
DimPlot(pbmc, reduction = "umap", group.by = "celltype",label = T)
dev.off()

#将全部的细胞注释结果、分类后的细胞注释结果输出
write.table(clusters,file="07.clusters cells500   features500  未进行第二次标准化.txt",quote=F,sep="\t",col.names=F)
write.table(celltype, file="07.celltype cells500   features500  未进行第二次标准化.txt",quote=F,sep="\t",col.names=F)
write.table(pbmc@meta.data,file="07.all cellAnn cells500   features500  未进行第二次标准化.txt",quote=F,sep="\t",col.names=F)




#准备monocle分析需要的文件   准备细胞轨迹分析需要的文件
#行名为基因名、列名为样品名的matrix文件
monocle.matrix=as.matrix(pbmc@assays$RNA@data)    #此步耗内存
monocle.matrix[1:6,1:6 ]
dim(monocle.matrix)
class(monocle.matrix)
#write.table(monocle.matrix, file = "07.monocleMatrix.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)##这一步耗时较长、耗内存较大。
#因输出txt文件较大，输出及后面读取耗时较久，尝试保存R文件
save(monocle.matrix,file="07.monocleMatrix  cells500   features500  未进行第二次标准化.Rdata")
#load("07.monocleMatrix  cells500   features500  未进行第二次标准化.Rdata")

#样品文件，样品类型，哪个细胞来自哪个样品
monocle.sample=as.matrix(pbmc@meta.data)
monocle.sample[1:6,1:6 ]
write.table(monocle.sample, file = "07.monocleSample  cells500   features500  未进行第二次标准化.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

#基因文件
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
head(monocle.geneAnn)
write.table(monocle.geneAnn,file="07.monocleGene  cells500   features500  未进行第二次标准化.txt",sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)   #所有基因文件    #不用输出了，直接接着做细胞轨迹图

#Mark基因文件
sig.markers07=as.data.frame(sig.markers)
sig.markers07[1:6,1:7 ]
sig.markers07$gene=rownames(sig.markers07)  #新增一列，列名为gene，值等于行名
sig.markers07[1:6,1:7]
write.table(sig.markers07, file = "07.monocleMarkers  cells500   features500  未进行第二次标准化.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


class(celltype)
head(celltype)
write.table(celltype, file = "07.monocleClusterAnn  cells500   features500  未进行第二次标准化.txt", sep = "\t", quote = FALSE, row.names = F, col.names = F)


#以下几行代码用于细胞互作、细胞通讯分析
#save(pbmc,file="pbmc 细胞通讯使用cells500   features500  未进行第二次标准化.Rdata")    #保存，用于后续细胞通讯互作分析
#load("pbmc 细胞通讯使用cells500   features500  未进行第二次标准化.Rdata")
#head(pbmc@meta.data)
#meta.data需要准备的文件=data.frame(pbmc@meta.data)
#head(meta.data需要准备的文件)
#write.table(meta.data需要准备的文件, file = "meta.data需要准备的文件cells500   features500  未进行第二次标准化.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
