library(circlize)      #引用包
geneFile="interGenes.List.txt"      #基因列表文件
posFile="geneREF.txt"               #基因位置信息文件
setwd("E:\\科研\\24.7.6 IA\\11 circlize")     #设置工作目录

#读取基因位置信息文件
genepos=read.table(posFile, header=T, sep="\t", check.names=F)
colnames(genepos)=c('genename','chr','start','end')
genepos=genepos[,c('chr','start','end','genename')]
row.names(genepos)=genepos[,'genename']

#读取基因列表文件
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
genepos=genepos[as.vector(geneRT[,1]),]
bed0=genepos

#绘制图形
pdf(file="circlize.pdf", width=6, height=6)
#初始化圈图
circos.clear()
circos.initializeWithIdeogram(species="hg38", plotType=NULL)
#展示每条染色体的注释信息
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col=rand_color(24))
  circos.text(mean(xlim), mean(ylim), chr, cex=0.6, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height=0.15, bg.border = NA)
#绘制基因组的图形
circos.genomicIdeogram(species = "hg38", track.height=mm_h(6))
#在染色体相应位置上标注基因的名称
circos.genomicLabels(bed0, labels.column=4, side = "inside", cex=0.8)
circos.clear()
dev.off()
