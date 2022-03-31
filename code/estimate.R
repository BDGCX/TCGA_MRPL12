##旧版timer数据库Gene里可以进行免疫浸润分析直接出图
###新版2.0使用详见https://zhuanlan.zhihu.com/p/447395807
###estimate包进行肿瘤免疫分析####
### estimate 包安装###
library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
###1.在UCSC Xena下载RNA-seq文件即exprSet
exprSet <- read.table(file = "KIRC-exprSet",sep = "\t",header = T)
rownames(exprSet) <- exprSet$GeneSymbol
exprSet <- exprSet[,-1]
###如果是RNA-SEQ的COUNTS矩阵要使用ESTIMATE的打包函数####
##只需要给一个RNA-seq的counts矩阵，如上所示，命名为exprSet，
##这样后续只需要使用我打包好的estimate函数即
dat=log2(edgeR::cpm(exprSet)+1)
library(estimate)
estimate <- function(dat,pro){
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(dat,file = input.f,sep = '\t',quote = F)
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina") ## 注意这个platform参数
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}
pro <- 'KIRC'
scores <- estimate(dat,pro)


###未使用
###2.进行Ensembl_id的转换，将Ensembl_id复制到excel里，
##使用分列功能分离小数点，然后再Ensembl数据库里进行转换详见
##https://www.bilibili.com/video/av806369373?from=search&seid=4814786577548386341&spm_id_from=333.337.0.0
mart_export <- read.table(file = "./mart_export.txt",sep =",",header = T )
write.csv(mart_export,file = "./genesymble.csv")
genesymbl <- read.csv(file = "./genesymble.csv",header = T)
colnames(genesymbl)[1] <-"Ensembl_ID" 
gtf_Ensembl_ID <- substr(exprSet[,1],1,15)
Ensembl_ID_To_Genename <- data.frame(Ensembl_ID=gtf_Ensembl_ID,exprSet)
Ensembl_ID_To_Genename <-Ensembl_ID_To_Genename[,-2]  
mergeRawCounts <- merge(genesymbl,Ensembl_ID_To_Genename,by="Ensembl_ID")
mergeRawCounts <- mergeRawCounts[,-c(1,3)]
colnames(mergeRawCounts)[1] <- "GeneSymbol"
rownames(mergeRawCounts) <- mergeRawCounts$GeneSymbol
mergeRawCounts <- mergeRawCounts[,-1]
# 5. ESTIMATE计算免疫得分
# 5.1 输入txt格式的表达矩阵，输出ESIMATE计算结果
filterCommonGenes(input.f= "./Rawdata/TCGA_HNSCpaired_Norexpr_data_paired.tsv", 
                  output.f="./Output/TCGA_estimate.gct", id="GeneSymbol")

## [1] "Merged dataset includes 9219 genes (1193 mismatched)."

estimateScore(input.ds = "./Output/TCGA_estimate.gct",
              output.ds = "./Output/TCGA_estimate_score.gct", 
              platform="illumina")

## [1] "1 gene set: StromalSignature  overlap= 135"
## [1] "2 gene set: ImmuneSignature  overlap= 139"

ESTI_score <- read.table("./Output/TCGA_estimate_score.gct",skip = 2,header = T,row.names = 1)
ESTI_score <- as.data.frame(t(ESTI_score[2:ncol(ESTI_score)]))
head(ESTI_score)
# 5.2 融合数据
table(row.names(ESTI_score) == rownames(phenotype))

## 
## TRUE 
##   86

ESTI_score$group <- phenotype$group
ESTI_score$sample <- rownames(ESTI_score)
ESTI_score_New =  melt(ESTI_score)

## Using group, sample as id variables

colnames(ESTI_score_New)=c("group","sample","status","score")  #设置行名
head(ESTI_score_New)

####多分组多平面带灰色框的图#####
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(ggpubr)
library(rstatix)
library(ggsignif)

iris %>% melt() %>%
  ggplot(aes(variable,value,fill=Species)) + 
  geom_boxplot()+
  scale_fill_jco()+
  stat_compare_means(label = "p.format")+
  facet_grid(.~variable,scales = "free",space="free_x")##产生标题灰框+
theme_bw()+theme(axis.text.x = element_blank())
