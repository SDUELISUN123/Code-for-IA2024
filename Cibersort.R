inputFile="merge.normalize.txt"      #表达数据文件
library(CIBERSORT)      #引用包

#免疫细胞浸润分析
outTab=cibersort("ref.txt", inputFile, perm=1000,qn)

#对免疫浸润结果过滤，并且保存免疫细胞浸润结果
#outTab=outTab[outTab[,"P-value"]<0.05,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)
