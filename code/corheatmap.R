IC <- read.csv(file = "Immune checkpoint.csv",header = F)
ic <- c(IC$V1)
IC <- exprSet[ic,]#没去除NA值
IC <- as.data.frame(t(IC))
MIC <- as.data.frame(t(exprSet["MRPL12",]))
# 构建相关关系矩阵
library(psych)
data.corr <- corr.test(IC, MIC, method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值

library(corrplot)
