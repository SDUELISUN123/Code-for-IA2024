#> 5000个细胞我的电脑都受不了，得用服务器

#devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r") 
library(CytoTRACE2)
library(tidyverse)
library(Seurat,lib.loc = "D:/R-4.3.1/library2/")
library(patchwork)
library(ggpubr)

# 1 准备输入文件 ###############
setwd("E:/科研/24.7.6 IA/19 scRNA")
sqy<-readRDS("sqy.rds")
sqy<-sqy[,sqy$celltype%in%"Mono/Macro"]






# 2 运行cytoTRACE ###########

#输入seurat 对象
sqy <- cytotrace2(sqy,
                  is_seurat = TRUE, 
                  slot_type = "counts", 
                  species = 'human',
                  seed = 12345)








# 3 可视化 （原作者默认的）###################
# making an annotation dataframe that matches input requirements for plotData function
annotation <- data.frame(phenotype = sqy@meta.data$seurat_clusters) %>% 
  set_rownames(., colnames(sqy))

# plotting
plots <- plotData(cytotrace2_result = sqy, 
                  annotation = annotation, 
                  is_seurat = TRUE,pc_dims = 30,seed = 12345)
# 绘制CytoTRACE2_Potency的umap图
p1 <- plots$CytoTRACE2_UMAP
# 绘制CytoTRACE2_Potency的umap图
p2 <- plots$CytoTRACE2_Potency_UMAP
# 绘制CytoTRACE2_Relative的umap图 ，v1 
p3 <- plots$CytoTRACE2_Relative_UMAP 
# 绘制各细胞类型CytoTRACE2_Score的箱线图
p4 <- plots$CytoTRACE2_Boxplot_byPheno

#输出
pdf("CytoTRACCE可视化.pdf",width = 10,height = 10)
(p1+p2+p3+p4) + plot_layout(ncol = 2)
dev.off()






# 4 可视化（自己写的代码）##########################

#Featureplot
pdf("Fea_relative_score.pdf",width = 5,height = 5)
FeaturePlot(sqy, "CytoTRACE2_Relative",pt.size = 0.5,reduction = "tsne") + 
  scale_colour_gradientn(colours = 
                           (c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", 
                              "#66C2A5", "#5E4FA2")), 
                         na.value = "transparent", 
                         limits = c(0, 1), 
                         breaks = seq(0, 1, by = 0.2), 
                         labels = c("0.0 (More diff.)", 
                                    "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"), 
                         name = "Relative\norder \n", 
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black")) + 
  ggtitle("CytoTRACE 2") + 
  xlab("TSNE1") + ylab("TSNE2") + 
  theme(legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        plot.title = element_text(size = 12, 
                                  face = "bold", hjust = 0.5, 
                                  margin = margin(b = 20))) + 
  theme(aspect.ratio = 1)
dev.off()

#Boxplot
p1 <- ggboxplot(sqy@meta.data, 
                x="seurat_clusters", y="CytoTRACE2_Score", width = 0.6, 
                color = "black",#轮廓颜色
                fill="seurat_clusters",#填充
                palette = "npg",
                xlab = F, #不显示x轴的标签
                bxp.errorbar=T,#显示误差条
                bxp.errorbar.width=0.5, #误差条大小
                size=1, #箱型图边线的粗细
                outlier.shape=NA, #不显示outlier
                legend = "right",
                title =  "CytoTRACE2 Score",
                order = c('1','5','6')) +
  NoLegend()+
  RotatedAxis()
###指定组比较
p<-p1+stat_compare_means(label.x = 1)
pdf("Vln_relative_score.pdf",width = 5,height = 4)
p
dev.off()

saveRDS(sqy,"sqyMM.rds")
