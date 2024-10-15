
rm(list=ls())
library(Seurat,lib.loc = "D:/R-4.3.1/library2/")
library(dplyr)
library(future)
library(future.apply)
library(monocle)
library(ggsci)
library(beepr)
library(ggplot2)
library(ggpubr)

plan("multicore", workers = 4) ###set the compute core
options(future.globals.maxSize = 15000 * 1024^2)
#getwd()

setwd("E:/科研/24.7.6 IA/19 scRNA")

###################### 用monocle做拟时序分析 #####################

sqy<-readRDS("sqyMM.rds")
#选定要拟时序分析的细胞种类
mmt<-sqy
rm(sqy)

gc()
#把TAM，monocyte等sub细胞分类信息给R
mmt$cell_type_val<-mmt$seurat_clusters
mmt@meta.data[1:5,]

#寻找高变基因，作为拟时序降维的基因
mmt<-FindVariableFeatures(mmt,nfeatures = 2000,)
#>这一步，不光是可以使用FindVariableFeatures寻找多变基因
#>也可以用前面的用于pca降维的2000个多变基因作为降维基因
#>二者都要试试，看哪个更符合生物学进程和假说，就用哪个
#>生信需要反复调整参数，而不是去跑一个流程就行
#>?FindVariableFeatures

#提取数据
matrix<-as.matrix(mmt@assays$RNA@counts)  #提取最原始的count数据
matrix<-rbind(matrix,"CytoTRACE2_Score" = mmt$CytoTRACE2_Score)
dim(matrix) #看多少行多少列
matrix[1:500,1:10]

#基因注释
gene_ann <- data.frame(
  gene_short_name = row.names(matrix), 
  row.names = row.names(matrix))

#细胞注释
#oct[["cell_group"]]<-Idents(oct)
mmt@meta.data[1:5,]
sample_ann <- mmt@meta.data

#建立monocle研究对象
fd <- new("AnnotatedDataFrame",data=gene_ann)
pd<-new("AnnotatedDataFrame",data=sample_ann)
sc_cds_2 <- newCellDataSet(matrix,  phenoData = pd,featureData =fd,expressionFamily = negbinomial.size(),lowerDetectionLimit=0.01)  #建立monocle研究对象，并进行归一化和质控
sc_cds_2 <- estimateSizeFactors(sc_cds_2)
sc_cds_2 <- estimateDispersions(sc_cds_2)

#质控
#?detectGenes：某基因必须在多少个细胞里有表达，才算有表达，并删掉没表达的，继续缩小数据
sc_cds_2 <- detectGenes(sc_cds_2, min_expr = 1)
expressed_genes <- row.names(subset(fData(sc_cds_2), num_cells_expressed >= 10)) #必须至少在10个细胞里（num_cells_expressed）有表达，才纳入该基因,可以自己调
fData(sc_cds_2)[1:5,]
#saveRDS(sc_cds_2,file = "9.06sc_cds_2.rds") #save一下sc_cds_2防止后边报错，闪退等
#sc_cds_2<-readRDS(file = "sc_cds_2.rds")

#设置高变基因
###以monocle下的FindVariableFeatures的高变基因为ordering
ordering_genes<-mmt@assays$RNA@var.features
###>或者
###>以组之间差异基因为ordering gene
#diff_test_res <- differentialGeneTest(sc_cds_2[expressed_genes,],
#                                      fullModelFormulaStr = "~cell_type_val",cores = 4,
#                                     verbose = T) #+num_genes_expressed+orig.ident
#ordering_genes <- row.names(subset(diff_test_res, qval < 1e-2)) #这个p值需要进行摸索，摸索到基因数在1000-2000
###>或者选logfc的top200作为高变基因
###>ordering_genes<-markers %>% group_by(cluster) %>% top_n(n=200,wt=avg_logFC)

#拟时序降维
sc_cds2 <- setOrderingFilter(sc_cds_2, ordering_genes)
#plot_ordering_genes(sc_cds2)  #看降维基因表达水平图？
sc_cds2<- reduceDimension(sc_cds2, max_components = 2, num_dim=6,reduction_method  = "DDRTree")#,residualModelFormulaStr = "~orig.ident")   #降到2维
sc_cds2 <- orderCells(sc_cds2,root_state = 3)  #把cell的顺序排出来、
saveRDS(sc_cds2,"monocle_for_plot.rds")
beepr::beep(1)



#看cytotrace score
plot_cell_trajectory(sc_cds2, color_by = "State",show_branch_points = T,cell_size = 1)#+facet_wrap(~group2) #拟时序，颜色表示细胞state
plot_cell_trajectory(sc_cds2, color_by = "CytoTRACE2_Score",show_branch_points = T,cell_size = 1) #拟时序，颜色表示假时间
phenoData <- sc_cds2@phenoData@data %>% as.data.frame()
ggplot(phenoData, aes(State, CytoTRACE2_Score,color=State))+
  geom_boxplot(width=.5)+
  stat_compare_means(label.x = 1.5)+
  theme_classic()+
  NoLegend()
ggsave("CytoTrace_Score.pdf",width = 4,height = 3)

plot_cell_trajectory(sc_cds2, color_by = "State",show_branch_points = T,cell_size = 1)#+facet_wrap(~group2) #拟时序，颜色表示细胞state
ggsave("traj_State.pdf",width = 4,height = 3)

plot_cell_trajectory(sc_cds2, color_by = "Pseudotime",show_branch_points = T)#拟时序，颜色表示细胞sub种类
ggsave("traj_time.pdf",width = 4,height = 3)

plot_genes_in_pseudotime(sc_cds2["Cyth4",],color_by = "Pseudotime",cell_size = .1)+
  theme_classic()
ggsave("cor_pse_cyth4.pdf",width = 4,height = 3)



#画想画的基因的热图
to_be_plot <- row.names(subset(fData(sc_cds2), gene_short_name %in% c("CANX","MKI67","PCNA","CCNB1")))
cds_subset1 <- sc_cds2[to_be_plot,]
plot_pseudotime_heatmap(cds_subset1,show_rownames = T,norm_method = c("log", "vstExprs"),num_clusters = 2)


#散点图（横轴，假时间；纵轴，自定义marker基因的表达量）
to_be_tested <- row.names(subset(fData(sc_cds2), gene_short_name %in% c("Cyth4","Mill1","Mill2")))
cds_subset2 <- sc_cds2[to_be_tested,]
plot_genes_in_pseudotime(cds_subset2,color_by = "Pseudotime")+
  stat_cor()+
  geom_smooth()



#寻找随着假时间变化的基因
diff_test_res <- differentialGeneTest(sc_cds2[expressed_genes,],fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 2)
sig_gene_names <- row.names(subset(diff_test_res, qval < 1))


#记得先看看想研究的基因是不是随时间变化的，再做后续分析
#输出随时间变化的基因
write.csv(sc_cds2[sig_gene_names,],"pseudotime_genes.csv")


#绘制随着假时间变化的基因表达的热图
load("24.8.30.RData")
pdf("pseudoheatmap.pdf",width = 5,height = 5)
pseudoplot<-plot_pseudotime_heatmap(sc_cds2[sig_gene_names,],
                                    num_clusters = 5,###进行调节分群，后面也要改
                                    cores = 5,
                                    show_rownames = F,return_heatmap = T)
dev.off()


#获得每个亚群的基因名称，以便于后续富集分析
clusters <- cutree(pseudoplot$tree_row, k = 5)###进行调节分群，前面也要改
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
write.csv(clustering,"cor_time_gene_cluster.csv")
beep(1)

