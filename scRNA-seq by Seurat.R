
rm(list=ls())
gc()
library(Seurat) #version 5
library(Seurat,lib.loc = "D:/R-4.3.1/library2/")
library(dplyr)
library(future)
library(future.apply)
library(msigdbr)
library(tidyverse)
library(patchwork)
library(monocle)
library(ggpubr)
library(devtools)
library(harmony)
library(SCpubr)
library(SeuratData)
library(batchelor)
library(SeuratWrappers)

setwd("E:\\科研\\24.7.6 IA\\19 scRNA")

cors<-ggsci::pal_igv(alpha = 0.5)(51)
cors2<-c(cors[2],cors[1],cors[3:51])
cors3<-ggsci::pal_frontiers()(10)

plan("multicore", workers = 3) ###set the compute core
options(future.globals.maxSize = 14000 * 1024^2)




# 1 读取数据 #######################################

## 循环读取h5文件 #######################
# 安装并加载包
#BiocManager::install("hdf5r")
library(hdf5r)
# 创建一个空的Seurat对象
seurat_list <- list()
# 设置文件夹路径
folder_list <- "GSE193533_RAW/"
# 获取文件夹中的h5文件列表
file_list <- list.files(path = folder_list, pattern = "\\.h5$", full.names = TRUE)
# 循环读取每个h5文件并转化为Seurat对象
for (file_path in file_list) {
  # 读取h5文件
  data <- Read10X_h5(file_path)
  # 提取文件名作为样本名称
  sample_name <- substr(basename(file_path), start = 1, stop = 10) ##提取文件名的第12到第15个字符作为seurat对象的名字
  # 创建seurat对象
  seurat_obj <- CreateSeuratObject(counts = data, project = sample_name)
  # 将 Seurat 对象添加到列表中
  seurat_list[[sample_name]] <- seurat_obj
}
# 合并Seurat对象
seurat_merge <- merge(seurat_list[[1]], y = c(seurat_list[[2]],seurat_list[[3]]), add.cell.ids = names(seurat_list))


saveRDS(seurat_merge,"最初的系数矩阵.rds")
rm(data,seurat_list,seurat_merge,seurat_obj,file_list,file_path,folder_list,sample_name)
gc()

# 4 标准化与PCA ################################
sqy<-readRDS("最初的系数矩阵.rds")

#使用PercentageFeatureSet函数计算线粒体基因的百分比
sqy[["percent.mt"]] <- PercentageFeatureSet(object = sqy, pattern = "^mt-",)#鼠源的换成mt

#计算核糖体基因百分比
sqy[["percent.rb"]] <- PercentageFeatureSet(sqy, pattern = "^Rp[Sl]")

#计算血红蛋白相关基因百分比
sqy[["percent.hb"]] <- PercentageFeatureSet(sqy, pattern = "^Hb[^(p)]")

#质控
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.hb")
pdf("质控.pdf",width = 12,height = 4)
VlnPlot(sqy, 
        fill.by = "feature", # "feature", "ident"
        features = feats,
        ncol = 5, pt.size = 0)
dev.off()
sqy=subset(x = sqy, 
           subset = percent.mt < 20 &
                    nFeature_RNA < 6000 &
                    nCount_RNA < 18000 &
                    percent.rb <30 &
                    percent.hb <1)    #对数据进行过滤
dim(sqy)

#测序深度的相关性图
plot1 <- FeatureScatter(object = sqy, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = sqy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
pdf("测序深度.pdf",width = 10,height = 4)
CombinePlots(plots = list(plot1, plot2))
dev.off()
rm(plot1,plot2)

#看细胞属于哪一期，并加到矩阵里
s.genes <- cc.genes$s.genes %>% str_to_title()
g2m.genes <- cc.genes$g2m.genes %>% str_to_title()
sqy <- CellCycleScoring(sqy, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
sqy@meta.data[1:5,]

#标准化
sqy<-NormalizeData(sqy,verbose = T,)   #标准化

#PCA
sqy<-FindVariableFeatures(sqy,selection.method = "vst", nfeatures = 2000)   #找前2000个差异显著的基因
sqy<-ScaleData(sqy,vars.to.regress = c("percent.mt","S.Score","G2M.Score"),verbose = T) #去除线粒体基因和分裂期的影响
sqy<-RunPCA(sqy,verbose = T,npcs = 50)  #pca降维
pdf("拐点pca.pdf")
ElbowPlot(sqy,ndims = 50)  #看拐点
dev.off()
#30







# 5 矫正和去批次 ##############################

## 2.5.1 Harmony ##############
#矫正前的降维图和vln
p1 <- DimPlot(object = sqy, reduction = "pca", pt.size = 1, group.by = "orig.ident",cols = cors2)
p2 <- VlnPlot(object = sqy, features = "PC_1", group.by = "orig.ident", pt.size = .1,cols = cors2)
pdf("矫正前pca.pdf",width = 10,height = 4)
p1|p2
dev.off()
#矫正
sqy<-RunHarmony(sqy,
                group.by.vars = c("orig.ident"),
                plot_convergence = TRUE,verbose = T)
harmony_embeddings <- Embeddings(sqy, 'harmony')
#矫正后的降维图和vln
p3 <- DimPlot(object = sqy, reduction = "harmony", pt.size = 1, group.by = "orig.ident",cols=cors2)
p4 <- VlnPlot(object = sqy, features = "harmony_1", group.by = "orig.ident", pt.size = .1,cols=cors2)
pdf("矫正后pca.pdf",width = 10,height = 4)
p3|p4
dev.off()
#清除环境
rm(p1,p2,p3,p4,harmony_embeddings,g2m.genes,s.genes,plot1,plot2)








# 6 添加临床信息 ###################################

Nomral_samples<-"GSM5813881"
IA_samples<-c("GSM5813883","GSM5813885")
sqy@meta.data$group[sqy$orig.ident%in%Nomral_samples]<-"Sham"
sqy@meta.data$group[sqy$orig.ident%in%IA_samples]<-"IA"

sqy@active.ident<-sqy$orig.ident %>% as.factor()
sqy<-RenameIdents(sqy,"GSM5813881"="Sham","GSM5813883"="Formed","GSM5813885"="Ruptured")
sqy$group2<-sqy@active.ident






# 8 去除双细胞 #####################################
library(DoubletFinder)

#UMAP/tsne
sqy <- FindNeighbors(object = sqy, dims = 1:30,reduction = "harmony")       #计算邻接距离
sqy <- FindClusters(object = sqy, resolution = 0.5,)         #对细胞分组,对细胞标准模块化,resolution在0.4-1.2，越大分的群越多
sqy <- RunTSNE(object = sqy, dims = 1:30,reduction = "harmony") 

## pK Identification (no ground-truth)
sweep.res.list <- paramSweep(sqy, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

## Homotypic Doublet Proportion Estimate
annotations <- sqy$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.075*nrow(sqy@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies
colnames(sqy@meta.data)[ncol(sqy@meta.data)]<-"pANN"
sqy <- doubletFinder(sqy, 
                     PCs = 1:15, 
                     pN = 0.25, 
                     pK = 0.09, 
                     nExp = nExp_poi, 
                     reuse.pANN = "pANN", 
                     sct = FALSE)

## Filtering
colnames(sqy@meta.data)[ncol(sqy@meta.data)]<-"DF.classifications"
pdf("Tsne_DF.classifications.pdf",width = 5,height = 4)
do_DimPlot(sqy,
           plot.title = "DF.classifications",
           reduction = "tsne",
           dims = c(1,2),
           group.by = "DF.classifications",
           pt.size = 0.5,
           legend.position = "top",
           label = F,
           label.box = F
)+
  theme_classic2()
dev.off()

table(sqy$DF.classifications)
sqy<-sqy[,sqy$DF.classifications%in%'Singlet']


## 9 UMAP降维 ################################

#UMAP/tsne
sqy <- FindNeighbors(object = sqy, dims = 1:30)      
sqy <- FindClusters(object = sqy, resolution = 0.5)         
sqy <- RunTSNE(object = sqy, dims = 1:30) 
#sqy<-RunTSNE(object = sqy, dims = 1:15)
#pdf(file = "tsne_cluster.pdf",width=5,height = 4)
#TSNEPlot(object = sqy, label = TRUE,cols=cors)#分组的umap 
#dev.off()

#各种UMAPPlot
pdf(file = "TSNE_cluster.pdf",width=5,height = 4)
do_DimPlot(sqy,
           plot.title = "Cluster",
           reduction = "tsne",
           dims = c(1,2),
           group.by = "seurat_clusters",
           pt.size = .5,
           legend.position = "top",
           label = T,
           label.box = T,
           repel = T
)+
  theme_classic2()
dev.off()

pdf(file = "TSNE_orig.ident.pdf",width=5,height = 4)
do_DimPlot(sqy,
           plot.title = "Sample",
           reduction = "tsne",
           dims = c(1,2),
           group.by = "orig.ident",
           pt.size = 0.5,
           legend.position = "top",
           label = T,
           label.box = F,
           repel = T
)+
  theme_classic2()
dev.off()

pdf(file = "UMAP_group2.pdf",width=5,height = 4)
do_DimPlot(sqy,
           plot.title = "Group",
           reduction = "tsne",
           dims = c(1,2),
           group.by = "group2",
           pt.size = 0.5,
           legend.position = "top",
           label = T,
           label.box = F,
           repel = T
)+
  theme_classic2()
dev.off()

colnames(sqy@meta.data)

pdf(file = "UMAP_phase.pdf",width=5,height = 4)
do_DimPlot(sqy,
           plot.title = "Phase",
           reduction = "tsne",
           dims = c(1,2),
           group.by = "Phase",
           pt.size = 0.5,
           legend.position = "top",
           label = T,
           label.box = F
)+
  theme_classic2()
dev.off()










# 10 注释细胞 #########################################

#查找每个聚类的差异基因(辅助注释)
logFCfilter=1
adjPvalFilter=0.05
sqy.markers <- FindAllMarkers(object = sqy,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter,)
sig.markers=sqy.markers[(abs(as.numeric(as.vector(sqy.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(sqy.markers$p_val_adj))<adjPvalFilter),]
write.csv(sig.markers,file="cluster_markers.csv")

top10 <- sqy.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file = "cluster_markers_top10.csv")

#第一次
genes <- list("immune"=c("Ptprc"),
              "epithelial"=c("Epcam"),
              "stromal"=c("Mme","Pecam1")
)

pdf(file = "Dot_第一次注释.pdf",width = 7,height = 7)
do_DotPlot(sample = sqy,features = genes,dot.scale = 12,
           colors.use = c("white","yellow","red"),
           legend.framewidth = 2, font.size =10,group.by = "seurat_clusters")
dev.off()

#第二次
genes <- list(
  'VSMC'=c('Myh11','Acta2'),
  'Mono/Macro'=c('Cd68','C1qa','C1qb'),
  'Neut'=c('S100a9','S100a8','Lcn2'),
  'Fibro'=c('Dcn'),
  'EC'=c('Cdh5','Pecam1'),
  'T/NK'=c('Cd3d','Nkg7')
)
sqy@active.ident<-sqy$seurat_clusters
pdf(file = "Dot_第二次注释.pdf",width = 10,height = 7)
do_DotPlot(sample = sqy,features = genes,dot.scale = 12,colors.use = c("white","yellow","red"),
           legend.framewidth = 2, font.size =10)
dev.off()

sqy@active.ident<-sqy$seurat_clusters
sqy<-RenameIdents(sqy,
                  '0'='VSMC',
                  '1'='Mono/Macro',
                  '2'='VSMC',
                  '3'='Neut',
                  '4'='VSMC',
                  '8'='Fibro',
                  '10'='EC',
                  '11'='T/NK',
                  '13'='Neut',
                  '5'='Mono/Macro',
                  '7'='VSMC',
                  '6'='Mono/Macro',
                  '9'='VSMC',
                  '12'='EC')
sqy$celltype<-sqy@active.ident%>%as.factor()
table(sqy$celltype)

##查找每种细胞的差异基因
logFCfilter=1
adjPvalFilter=0.05
sqy.markers <- FindAllMarkers(object = sqy,
                              only.pos = FALSE,
                              min.pct = 0.25,
                              logfc.threshold = logFCfilter,)
sig.markers=sqy.markers[(abs(as.numeric(as.vector(sqy.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(sqy.markers$p_val_adj))<adjPvalFilter),]
write.csv(sig.markers,file="cluster_markers_celltype.csv")

top10 <- sqy.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file = "cluster_markers_top10_celltype.csv")

#看各群的火山图
library(scRNAtoolVis)
pdf("jjvol_各种细胞的标志物.pdf",width = 9,height = 6)
jjVolcano(diffData = sqy.markers, topGeneN = 5,)+
  xlab("Celltype")+
  scale_color_manual(values = c('#2500ff','#dd5000'))+
  ggtitle("Top-5 markers of each cell-type")
dev.off()

#绘制marker在各个celltype的热图
pdf(file = "doheatmap_cluster_top10.pdf",width = 15, height = 15)
DoHeatmap(object = sqy, features = top10$gene) + NoLegend()
dev.off()

#celltype UMAP
pdf(file = "UMAP_celltype.pdf",width=5,height = 4)
do_DimPlot(sqy,
           plot.title = "Celltype",
           reduction = "tsne",
           dims = c(1,2),
           group.by = "celltype",
           pt.size = 0.3,
           legend.position = "top",
           label = T,
           label.box = F,
           repel = T
)+
  theme_classic2()
dev.off()
pdf(file = "UMAP_celltype_split.pdf",width=7,height = 6)
do_DimPlot(sqy,
           plot.title = "Celltype",
           reduction = "tsne",
           dims = c(1,2),
           group.by = "celltype",
           pt.size = 0.1,
           legend.position = "top",
           label = T,
           label.box = F,
           split.by = "group2"
)+
  theme_classic2()
dev.off()

# celltype Dotplot

pdf(file = "Dot_marker_cluster.pdf",width = 10,height = 5)
do_DotPlot(sample = sqy,features = genes,dot.scale = 12,colors.use = c("white","yellow","red"),
           legend.framewidth = 2, font.size =10,group.by = "celltype")
dev.off()









# 11 看每种细胞在每个样本中的比例 ################################

##### 1 堆砌柱状图
library(ggplot2)
tab<-table(Idents(sqy),sqy$orig.ident)%>%as.data.frame()
write.csv(tab,file = "每个样本各细胞数目.csv")

# 1.1绘制绝对数的柱状图
ggplot(tab,aes(x=Var2,y=Freq,fill=Var1))+ 
  geom_bar(stat='identity',position='stack',alpha=.5)+ 
  labs(title='',x='group',y='Cell Number')+ 
  theme(legend.justification = 'right', 
        legend.position = 'right', 
        legend.key.height = unit(0.1,'cm'),
        panel.background = element_blank(),
        axis.line=element_line(size=0.5,colour="black")
  )+ scale_fill_manual(values=cors)

# 1.2绘制百分比的柱状图
tab<-table(Idents(sqy),sqy$orig.ident)%>%
  prop.table(.,2)*100
tab<-tab%>%as.data.frame()


pdf("百分比堆砌柱状图，每个样本.pdf",width = 4,height = 4)
ggplot(tab,aes(x=Var2,y=Freq,fill=Var1))+ 
  geom_bar(stat='identity',position='stack',alpha=.5)+ 
  labs(title='',x='group',y='Cell Proportion (%)')+ 
  theme(legend.justification = 'right', 
        legend.position = 'right', 
        legend.key.height = unit(0.1,'cm'),
        panel.background = element_blank(),
        axis.line=element_line(size=0.5,colour="black")
  )+ scale_fill_manual(values=cors)+
  RotatedAxis()

dev.off()


pdf("Fea_exp.pdf",width = 12,height = 4)
FeaturePlot(sqy,features = c("Cyth4","Mill1","Mill2"),reduction = 'tsne',ncol = 3)
dev.off()

pdf("Vln_exp.pdf",width = 6,height = 12)
do_ViolinPlot(sqy,features = c("Cyth4","Mill1","Mill2"),
              group.by = "celltype",
              pt.size = 0,
              plot_boxplot = T,
              ncol = 1)
dev.off()

save.image("24.8.11.RData")
saveRDS(sqy,"sqy.rds")
rm(list=ls()) ; gc()


# 12 验证表达量 #######

P1<-VlnPlot(mmt,features = "Cyth4",group.by = "group",pt.size = -1)+
  stat_compare_means(method = "wilcox")+
  geom_boxplot(width=0.1,color="black",fill="white")+
  theme_classic2()+
  xlab(NULL)+
  NoLegend()
  
P2<-VlnPlot(mmt,features = "Mill1",group.by = "group",pt.size = -1)+
  stat_compare_means(method = "wilcox")+
  geom_boxplot(width=0.1,color="black",fill="white")+
  theme_classic2()+
  xlab(NULL)+
  NoLegend()

P3<-VlnPlot(mmt,features = "Mill2",group.by = "group",pt.size = -1)+
  stat_compare_means(method = "wilcox")+
  geom_boxplot(width=0.1,color="black",fill="white")+
  theme_classic2()+
  xlab(NULL)+
  NoLegend()

pdf("vln_验证基因表达.pdf",width = 8,height = 3)
P1|P2|P3
dev.off()
