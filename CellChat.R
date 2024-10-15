
library(Seurat,lib.loc = "D:/R-4.3.1/library2/")
library(dplyr)
library(future)
library(future.apply)
library(stringr)
library(CellChat)
library(patchwork)

plan("multicore", workers = 6) ###set the compute core
options(future.globals.maxSize = 160000 * 1024^2)

cors<-ggsci::pal_igv()(20)

#setwd("E:\\科研\\230712 糖足胞葬\\2.6 singlecell2\\1 foot\\融合\\cellchat2")


# 0 准备 ####

sce.mergeTEN<-readRDS("./sqy.rds")
#sce.mergeTEN<-sce.mergeTEN[,Idents(sce.mergeTEN)!="Doublets"] #去除双细胞
#>by the way，双细胞我不喜欢去计算，因为费时间，并且见仁见智，此双细胞并非彼双细胞
#>我见过奇葩的TCR阳性的macrophage，现实中竟然还真存在这种细胞

##把所有数据拆分成两个，normal和tumor
normaldata<-sce.mergeTEN[,sce.mergeTEN$group=="Sham"]
tumordata<-sce.mergeTEN[,sce.mergeTEN$group=="IA"]

#试运行时（比如学习时）可以只取前2000细胞
#normaldata<-normaldata[,1:2000]
#tumordata<-tumordata[,1:2000]



# 1读入数据，并保存rds文件######

#devtools::install_github("sqjin/CellChat")


##normal
table(sce.mergeTEN$celltype)#看各种细胞类型的名字和数目，以便后边改（要去掉没有的细胞）
#sqy<-sqy[,sqy$cell_type!=c("RBCs")]  #去掉没有的细胞，防止报错，慎重覆盖变量！
normal.input <- GetAssayData(normaldata, assay = "RNA", slot = "data") # normalized data matrix
labels <- factor(normaldata$celltype,levels = levels(Idents(normaldata))) #用哪一列当细胞类型
meta <- data.frame(group = labels, row.names = rownames(normaldata@meta.data)) # create a dataframe of the cell labels
meta$group<-str_replace_all(meta$group,pattern = "/",replacement = "_") #把斜杠替换成下划线，以免报错
#下句话记得把上面各种细胞类型粘上，然后把"/"改成"_"，并且不能放数目为0的细胞，防报错
meta$group <- factor(meta$group)
normal_cellchat <- createCellChat(object = normal.input, meta = meta, group.by = "group")
table(normal_cellchat@idents)  #看有没有缺失值NA和数目为0的细胞
saveRDS(normal_cellchat,file="normal_cellchat_original.rds")
rm(normaldata,normal.input,labels,meta)

##normal
table(sce.mergeTEN$celltype)#看各种细胞类型的名字和数目，以便后边改（要去掉没有的细胞）
#sqy<-sqy[,sqy$cell_type!=c("RBCs")]  #去掉没有的细胞，防止报错，慎重覆盖变量！
tumor.input <- GetAssayData(tumordata, assay = "RNA", slot = "data") # normalized data matrix
labels <- factor(tumordata$celltype,levels = levels(Idents(tumordata))) #用哪一列当细胞类型
meta <- data.frame(group = labels, row.names = rownames(tumordata@meta.data)) # create a dataframe of the cell labels
meta$group<-str_replace_all(meta$group,pattern = "/",replacement = "_") #把斜杠替换成下划线，以免报错
#下句话记得把上面各种细胞类型粘上，然后把"/"改成"_"，并且不能放数目为0的细胞，防报错
meta$group <- factor(meta$group)
tumor_cellchat <- createCellChat(object = tumor.input, meta = meta, group.by = "group")
table(tumor_cellchat@idents)  #看有没有缺失值NA和数目为0的细胞
saveRDS(tumor_cellchat,file="tumor_cellchat_original.rds")
rm(tumordata,normal.input,labels,meta)

####特别说明，cellchat支持seurat格式的文件，所以可以直接用createcellchat函数处理seurat文件
#例子：
#normal_cellchat<-createCellChat(sqy)




# 2 数据处理#####


rm(list=ls()) #清除缓存
gc()
normal_cellchat<-readRDS(file="normal_cellchat_original.rds")

#预处理要跑十分钟左右，去休息一下吧！
normal_cellchat@DB <- CellChatDB.mouse#物种：人类,可以改为“CellChatDB.mouse”
#View(CellChatDB.human$interaction)  #看有哪些互作
normal_cellchat <- subsetData(normal_cellchat) # subset the expression data of signaling genes for saving computation cost
normal_cellchat <- identifyOverExpressedGenes(normal_cellchat)  #确定每个细胞亚群中的过表达基因
normal_cellchat <- identifyOverExpressedInteractions(normal_cellchat) #寻找过表达的interaction
normal_cellchat <- projectData(normal_cellchat, PPI.human)  #向ppi投射
normal_cellchat <- computeCommunProb(normal_cellchat,raw.use = T) #算interaction的可能性
normal_cellchat <- filterCommunication(normal_cellchat, min.cells = 10)  #去除interaction很少的细胞
normal_cellchat <- computeCommunProbPathway(normal_cellchat)  #计算通路
normal_cellchat <- netAnalysis_computeCentrality(normal_cellchat, slot.name = "netP") #信号网络的一个拓扑学分析

##保存和读取
saveRDS(normal_cellchat,file="normal_cellchat.rds")
normal_cellchat<-readRDS("normal_cellchat.rds")

rm(list=ls()); gc()
tumor_cellchat@DB <- CellChatDB.mouse#物种：人类,可以改为“CellChatDB.mouse”
#View(CellChatDB.human$interaction)  #看有哪些互作
tumor_cellchat <- subsetData(tumor_cellchat) # subset the expression data of signaling genes for saving computation cost
tumor_cellchat <- identifyOverExpressedGenes(tumor_cellchat)  #确定每个细胞亚群中的过表达基因
tumor_cellchat <- identifyOverExpressedInteractions(tumor_cellchat) #寻找过表达的interaction
tumor_cellchat <- projectData(tumor_cellchat, PPI.human)  #向ppi投射
tumor_cellchat <- computeCommunProb(tumor_cellchat,raw.use = T) #算interaction的可能性
tumor_cellchat <- filterCommunication(tumor_cellchat, min.cells = 10)  #去除interaction很少的细胞
tumor_cellchat <- computeCommunProbPathway(tumor_cellchat)  #计算通路
tumor_cellchat <- netAnalysis_computeCentrality(tumor_cellchat, slot.name = "netP") #信号网络的一个拓扑学分析

##保存和读取
saveRDS(tumor_cellchat,file="tumor_cellchat.rds")
tumor_cellchat<-readRDS("tumor_cellchat.rds")

rm(list=ls()); gc()


# 3 看差异 #######

library(CellChat)
library(ComplexHeatmap)
library(patchwork)

cellchat.NL <- readRDS("normal_cellchat.rds")
cellchat.LS <- readRDS("tumor_cellchat.rds")

#查看两个数据集细胞分群，并检查是否一致：
levels(cellchat.NL@idents)
levels(cellchat.LS@idents)
identical(levels(cellchat.NL@idents),levels(cellchat.LS@idents))

#> 注意，本期cellchat组间差异分析流程仅适用于聚类后细胞群组成相同的不同数据集，且仅支持两组别的通讯水平差异比较。

#合并cellchat对象:
object.list <- list(NL = cellchat.NL,
                    LS = cellchat.LS) #对照组(NL)在前，比较组(LS)在后，注意顺序
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

dplyr::glimpse(cellchat)



## 1 总体比较：通讯数目与通讯强度差异 #####

p1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
p2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p1 + p2

## 2.细胞亚群水平的通讯差异 #####

par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T,)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# 3.传出信号 #####
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways,
                       object.list[[i+1]]@netP$pathways)
pathway.union

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]],
                                        pattern = "outgoing", #传出
                                        signaling = pathway.union,
                                        title = names(object.list)[i],
                                        width = 4,
                                        height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]],
                                        pattern = "outgoing", #传出
                                        signaling = pathway.union,
                                        title = names(object.list)[i+1],
                                        width = 4,
                                        height = 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

# 4.传入信号 #####
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]],
                                        pattern = "incoming", #传入
                                        signaling = pathway.union,
                                        title = names(object.list)[i],
                                        width = 4, height = 15,
                                        color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+1]],
                                        pattern = "incoming", #传入
                                        signaling = pathway.union,
                                        title = names(object.list)[i+1],
                                        width = 4, height = 15,
                                        color.heatmap = "GnBu")
draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))

(ht1+ht2)+(ht3+ht4)

# 4.1 总配受体对概率差异气泡图 ####
levels(cellchat@idents$joint) #查看细胞亚群
netVisual_bubble(cellchat,
                 sources.use = 3,
                 targets.use = c(1:6),
                 comparison = c(1, 2),
                 angle.x = 45)
# 4.2 区分上下调配体对 ####
p7 <- netVisual_bubble(cellchat,
                       sources.use = 3,
                       targets.use = c(1:6),
                       comparison = c(1,2),
                       max.dataset = 2,
                       title.name = "Increased signaling in LS",
                       angle.x = 45,
                       remove.isolate = T) #Increased为比较组通讯概率更强的配受体对信息
p8 <- netVisual_bubble(cellchat,
                       sources.use = 3,
                       targets.use = c(1:6),
                       comparison = c(1,2),
                       max.dataset = 1,
                       title.name = "Decreased signaling in LS",
                       angle.x = 45,
                       remove.isolate = T) #Decreased为对照组通讯概率更强的配受体对信息
p7 + p8

# 单个 ####
#使用网络图：
pathway.union
pathways.show <- c("CXCL") #选择目标信号通路
weight.max <- getMaxWeight(object.list,
                           slot.name = c("netP"),
                           attribute = pathways.show) #控制不同数据集的边权重

  par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]],
                      signaling = pathways.show,
                      layout = "circle",
                      edge.weight.max = weight.max[1],
                      edge.width.max = 10,
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}


# 看Cyth4和Ccr1 ####
sqy <- readRDS("E:/科研/24.7.6 IA/19.7 cellchat/sqy.rds")
sqy_Cyth4 <- sqy@assays$RNA@data["Cyth4",] %>% as.data.frame()
colnames(sqy_Cyth4) = "Cyth4"
sqy_Ccr1 <- sqy@assays$RNA@data["Tnf",] %>% as.data.frame()
colnames(sqy_Ccr1) = "Ccr1"

sqy_df<- cbind(sqy_Cyth4,sqy_Ccr1)
#sqy_df <- sqy_df %>% filter(Cyth4>0 & Ccr1>0)

library(ggplot2); library(ggpubr)
ggplot(sqy_df,aes(x=Cyth4,y=Ccr1))+
  geom_point()+
  stat_cor()+
  geom_smooth(method = lm)

FeatureScatter(sqy,feature1 = "Cyth4",feature2 = "")+
  theme_classic()+
  stat_cor()+
  geom_smooth(method = lm)+
  ggtitle("Correlation between Cyth4 and Tnf")

FeatureScatter(sqy[,sqy$celltype%in%"Mono/Macro"],feature1 = "Cyth4",feature2 = "Ccl5",slot = "count")+
  theme_classic()+
  stat_cor(method = "spearman")+
  geom_smooth(method = loess)+
  ggtitle("Correlation between Cyth4 and Tnf")

VlnPlot(sqy[,sqy$celltype%in%"Mono/Macro"],features = "Cyth4",group.by = "group",pt.size = -1)+
  stat_compare_means(method = "wilcox")+
  geom_boxplot(width=0.1,color="black",fill="white",outlier.size = 0)+
  theme_classic2()+
  xlab(NULL)+
  NoLegend()



# 5 出图 ####

# 1 总览细胞通讯
par(mfrow = c(1,2), xpd = TRUE)
normal_cellchat <- aggregateNet(cellchat.NL)
groupSize <- as.numeric(table(normal_cellchat@idents))
netVisual_circle(normal_cellchat@net$weight, arrow.size = 0.01,vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength") #看细胞通讯的weight权重

tumor_cellchat <- aggregateNet(cellchat.LS)
groupSize <- as.numeric(table(tumor_cellchat@idents))
netVisual_circle(tumor_cellchat@net$weight, arrow.size = 0.01,vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength") #看细胞通讯的weight权重

# 2 传出和传入信号差异
pdf("2.传出和传入信号差异.pdf",width = 12,height = 8)
(ht1+ht2)+(ht3+ht4)
dev.off()

# 3 总配受体对概率差异气泡图
pdf("3.总配受体对概率差异气泡图.pdf",width = 8,height = 8)
netVisual_bubble(cellchat,
                 sources.use = 3,
                 targets.use = c(1:6),
                 comparison = c(1, 2),
                 angle.x = 45)
dev.off()

# 4 上下调气泡图 
pdf("4.上下调气泡图.pdf",width = 8,height = 8)
p7+p8
dev.off()

# 5 Cyth4在巨噬细胞中和这几个ccl基因的关系 
sqy_Cyth4 <- sqy[,sqy$celltype%in%"Mono/Macro"]@assays$RNA@counts["Cyth4",] %>% as.data.frame()
colnames(sqy_Cyth4) = "Cyth4"
sqy_Ccl3 <- sqy[,sqy$celltype%in%"Mono/Macro"]@assays$RNA@counts["Ccl3",] %>% as.data.frame()
colnames(sqy_Ccl3) = "Ccl3"
sqy_Ccl5 <- sqy[,sqy$celltype%in%"Mono/Macro"]@assays$RNA@counts["Ccl5",] %>% as.data.frame()
colnames(sqy_Ccl5) = "Ccl5"
sqy_Ccl6 <- sqy[,sqy$celltype%in%"Mono/Macro"]@assays$RNA@counts["Ccl6",] %>% as.data.frame()
colnames(sqy_Ccl6) = "Ccl6"
sqy_Ccl7 <- sqy[,sqy$celltype%in%"Mono/Macro"]@assays$RNA@counts["Ccl7",] %>% as.data.frame()
colnames(sqy_Ccl7) = "Ccl7"
sqy_Ccl9 <- sqy[,sqy$celltype%in%"Mono/Macro"]@assays$RNA@counts["Ccl9",] %>% as.data.frame()
colnames(sqy_Ccl9) = "Ccl9"

df <- cbind(sqy_Cyth4,sqy_Ccl3,sqy_Ccl5,sqy_Ccl6,sqy_Ccl7,sqy_Ccl9)
write.csv(df,"Cyth4和ccl基因的count.csv",row.names = T,quote = F)

# 6 ccl基因在macro中上调

v3<-VlnPlot(sqy[,sqy$celltype%in%"Mono/Macro"],features = "Ccl3",group.by = "group",pt.size = -1)+
  stat_compare_means(method = "wilcox")+
  geom_boxplot(width=0.1,color="black",fill="white",outlier.size = 0)+
  theme_classic2()+
  xlab(NULL)+
  NoLegend()
v5<-VlnPlot(sqy[,sqy$celltype%in%"Mono/Macro"],features = "Ccl5",group.by = "group",pt.size = -1)+
  stat_compare_means(method = "wilcox")+
  geom_boxplot(width=0.1,color="black",fill="white",outlier.size = 0)+
  theme_classic2()+
  xlab(NULL)+
  NoLegend()
v6<-VlnPlot(sqy[,sqy$celltype%in%"Mono/Macro"],features = "Ccl6",group.by = "group",pt.size = -1)+
  stat_compare_means(method = "wilcox")+
  geom_boxplot(width=0.1,color="black",fill="white",outlier.size = 0)+
  theme_classic2()+
  xlab(NULL)+
  NoLegend()
v7<-VlnPlot(sqy[,sqy$celltype%in%"Mono/Macro"],features = "Ccl7",group.by = "group",pt.size = -1)+
  stat_compare_means(method = "wilcox")+
  geom_boxplot(width=0.1,color="black",fill="white",outlier.size = 0)+
  theme_classic2()+
  xlab(NULL)+
  NoLegend()
v9<-VlnPlot(sqy[,sqy$celltype%in%"Mono/Macro"],features = "Ccl9",group.by = "group",pt.size = -1)+
  stat_compare_means(method = "wilcox")+
  geom_boxplot(width=0.1,color="black",fill="white",outlier.size = 0)+
  theme_classic2()+
  xlab(NULL)+
  NoLegend()

v3|v5|v6|v7|v9
ggsave("6.ccl基因在macro中的上调.pdf",width = 12,height = 3)
