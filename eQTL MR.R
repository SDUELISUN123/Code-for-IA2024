

###### 1 去除弱工具变量之后，做MR分析 #####



#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("VariantAnnotation")

#install.packages("devtools")
#devtools::install_github("mrcieu/gwasglue", force = TRUE)

#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")


#引用包
library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)
library(future)

plan("multicore", workers = 3) ###set the compute core
options(future.globals.maxSize = 14000 * 1024^2)

rm(list=ls())
gc()
exposureFile="~/sqy/24.7.14 IA/7 MR/exposure.F.csv"       #暴露数据文件
outcomeID="finn-b-I9_SAHANEUR"        #结局数据id(需修改)
outcomeName="Aneurysms, operations, SAH"       #图形中展示疾病的名称
path<-paste0("~/sqy/24.7.14 IA/7 MR/",outcomeID)
dir.create(path)
setwd(path)     #设置工作目录

#读取暴露数据
exposure_dat=read_exposure_data(filename=exposureFile,
                                sep = ",",
                                snp_col = "SNP",
                                beta_col = "beta.exposure",
                                se_col = "se.exposure",
                                pval_col = "pval.exposure",
                                effect_allele_col="effect_allele.exposure",
                                other_allele_col = "other_allele.exposure",
                                eaf_col = "eaf.exposure",
                                phenotype_col = "exposure",
                                id_col = "id.exposure",
                                samplesize_col = "samplesize.exposure",
                                chr_col="chr.exposure", pos_col = "pos.exposure",
                                clump=FALSE)

#读取结局数据
outcomeData=extract_outcome_data(snps=factor(exposure_dat$SNP), outcomes=outcomeID)
write.csv(outcomeData, file="outcome.csv", row.names=F)

#将暴露数据和结局数据合并
outcomeData$outcome=outcomeName
dat=harmonise_data(exposure_dat, outcomeData)

#
dat$samplesize.outcome=195203
dat2 <- subset(dat, dat$exposure %in% c("MICA"))
dat3 <- subset(dat, dat$exposure %in% c("CYTH4"))
mr_report(dat2)
mr_report(dat3)

#输出用于孟德尔随机化的工具变量
outTab=dat[dat$mr_keep=="TRUE",]
write.csv(outTab, file="table.SNP.csv", row.names=F)


#孟德尔随机化分析
start_time<-Sys.time()
mrResult=mr(dat)
end_time<-Sys.time()
MR_time<-end_time-start_time
#beepr::beep(1)
MR_time

#对结果进行OR值的计算
mrTab=generate_odds_ratios(mrResult)
write.csv(mrTab, file="table.MRresult.csv", row.names=F)

#异质性分析
heterTab=mr_heterogeneity(dat)
write.csv(heterTab, file="table.heterogeneity.csv", row.names=F)

#多效性检验
pleioTab=mr_pleiotropy_test(dat)
write.csv(pleioTab, file="table.pleiotropy.csv", row.names=F)

#绘制散点图
pdf(file="pic.scatter_plot.pdf", width=7.5, height=7)
mr_scatter_plot(mrResult, dat)
dev.off()

#森林图
res_single=mr_singlesnp(dat)      #得到每个工具变量对结局的影响
pdf(file="pic.forest.pdf", width=7, height=5.5)
mr_forest_plot(res_single)
dev.off()

#漏斗图
pdf(file="pic.funnel_plot.pdf", width=7, height=6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()

#留一法敏感性分析
pdf(file="pic.leaveoneout.pdf", width=7, height=5.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
dev.off()


##### 2 IVW算法过滤掉不相关和基因 ####

mrFile="table.MRresult.csv"        #孟德尔随机化分析的结果文件
pleFile="table.pleiotropy.csv"     #多效性的结果文件
#setwd("E:/科研/24.2.29 pan_neuro_MR/1 pd/13.IVWfilter")    #设置工作目录

#读取孟德尔随机化的结果文件
rt=read.csv(mrFile, header=T, sep=",", check.names=F)
#提取IVW方法pvalue<0.05的基因
ivw=data.frame()
for(geneName in unique(rt$exposure)){
  geneData=rt[rt$exposure==geneName,]
  #提取五种方法OR方向一致的基因
  if(nrow(geneData)==5){
    if(geneData[geneData$method=="Inverse variance weighted","pval"]<0.05){
      if(sum(geneData$or>1)==nrow(geneData) | sum(geneData$or<1)==nrow(geneData)){
        ivw=rbind(ivw, geneData)
      }
    }
  }
}

#读取多效性的结果文件
pleRT=read.csv(pleFile, header=T, sep=",", check.names=F)
#剔除多效性pvalue小于0.05的基因
pleRT=pleRT[pleRT$pval>0.05,]
immuneLists=as.vector(pleRT$exposure)
outTab=ivw[ivw$exposure %in% immuneLists,]
write.csv(outTab, file="IVW.filter.csv", row.names=F)

#输出最终的基因
MR_gene<-unique(outTab$exposure)
write.table(MR_gene,"MR_gene.txt",quote = F,row.names = F,col.names = F)

save.image(paste0("result_",outcomeID,".RData"))

