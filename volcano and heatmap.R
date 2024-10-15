#引用包
library(limma)
library(dplyr)
library(pheatmap)
library(ggplot2)

logFCfilter=0.5          #logFC的过滤条件
adj.P.Val.Filter=0.05      #矫正后p值的过滤条件
inputFile="merge.normalize.txt"      #表达数据文件
setwd("E:\\科研\\24.7.6 IA\\6 batchDiff")      #设置工作目录

#读取输入文件，并对输入文件整理
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#获取样品的分组信息
Type=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data=data[,order(Type)]       #根据样品的分组信息对样品进行排序
Project=gsub("(.+)\\_(.+)\\_(.+)", "\\1", colnames(data))    #获得GEO数据研究的id
Type=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
colnames(data)=gsub("(.+)\\_(.+)\\_(.+)", "\\2", colnames(data))

#差异分析
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("Control","Treat")
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(Treat-Control,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

#输出所有基因的差异情况
allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, file="all.txt", sep="\t", quote=F, col.names=F)

#输出显著的差异基因
diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="diff.txt", sep="\t", quote=F, col.names=F)

#输出差异基因表达量
diffGeneExp=data[row.names(diffSig),]
diffGeneExpOut=rbind(id=paste0(colnames(diffGeneExp),"_",Type), diffGeneExp)
write.table(diffGeneExpOut, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)

#绘制差异基因热图
geneNum=20     #定义展示基因的数目
diffUp=diffSig[diffSig$logFC>0,]
diffDown=diffSig[diffSig$logFC<0,]
geneUp=row.names(diffUp)
geneDown=row.names(diffDown)
if(nrow(diffUp)>geneNum){geneUp=row.names(diffUp)[1:geneNum]}
if(nrow(diffDown)>geneNum){geneDown=row.names(diffDown)[1:geneNum]}
hmExp=data[c(geneUp,geneDown),]
#准备注释文件
names(Type)=colnames(data)
Type=as.data.frame(Type)
Type=cbind(Project, Type)
#输出图形
pdf(file="heatmap.pdf", width=5, height=5)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue2", "white", "red2"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=5.5,
         fontsize_col=8)
dev.off()


#定义显著性
rt=read.table("all.txt", header=T, sep="\t", check.names=F)
Sig=ifelse((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")
#绘制火山图
rt = mutate(rt, Sig=Sig)
p = ggplot(rt, aes(logFC, -log10(adj.P.Val)))+
    geom_point(aes(col=Sig))+
    scale_color_manual(values=c("blue2", "grey","red2"))+
    labs(title = " ")+
    theme(plot.title = element_text(size=16, hjust=0.5, face = "bold"))
#输出图形
pdf(file="vol.pdf", width=5.5, height=4.5)
print(p)
dev.off()

