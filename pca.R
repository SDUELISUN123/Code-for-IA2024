#???ð?
library(ggplot2)
library(ggpubr)
setwd("E:\\????\\24.7.6 IA\\5 pca")    #???ù???Ŀ¼

#????PCA?????ĺ???
bioPCA=function(inputFile=null, outFile=null, titleName=null){
	#??ȡ?????ļ?,??ȡ????
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
	data=t(rt)
	Project=gsub("(.*?)\\_.*", "\\1", rownames(data))    #??ȡGEO???ݿ??о???id
	
	#PCA????
	data.pca=prcomp(data)
	pcaPredict=predict(data.pca)
	PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=Project)

	#????ͼ??
	pdf(file=outFile, width=5.5, height=4.25)
	p1=ggscatter(data=PCA, x="PC1", y="PC2", color="Type", shape="Type", 
	          ellipse=T, ellipse.type="norm", ellipse.border.remove=F, ellipse.alpha = 0.1,
	          size=2, main=titleName, legend="right")+
	          theme(plot.margin=unit(rep(1.5,4),'lines'), plot.title = element_text(hjust=0.5))
	print(p1)
	dev.off()
}

#???ú???, ???????ν???ǰ??ͼ??
bioPCA(inputFile="merge.preNorm.txt", outFile="PCA.preNorm.pdf", titleName="Before batch correction")
#???ú???, ???????ν???????ͼ??
bioPCA(inputFile="merge.normalize.txt", outFile="PCA.normalzie.pdf", titleName="After batch correction")

# boxplot
pre_matrix <- read.table("merge.preNorm.txt",header = T,sep = "\t",row.names = 1)
pdf("pre_normal.pdf",width = 15,height = 5)
boxplot(pre_matrix,outline=T, notch=T, las=2)
dev.off()


after_matrix <- read.table("merge.normalize.txt",header = T,sep = "\t",row.names = 1)
pdf("after_normal.pdf",width = 15,height = 5)
boxplot(after_matrix,outline=T, notch=T, las=2)
dev.off()
