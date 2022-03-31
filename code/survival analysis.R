library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
sa <- read.csv(file = "LIHC_6182_29_71.csv",header = T)
sa <- sa[-which(sa$Days>3000),]
sa$Status <- ifelse(sa$Status=='Alive',0,1)
# 利用ggsurvplot快速绘制漂亮的生存曲线图
sfit <- survfit(Surv(Days, Status)~Group, data=sa)
sfit
summary(sfit)
my.surv <- Surv(sa$Days,sa$Status=="1")
survdiff(my.surv~Group,data = sa)
coxmodel <- coxph(my.surv~Group,data = sa)
summary(coxmodel)###分别计算每种肿瘤的HR和P值
## more complicate figures.
ggsurvplot(sfit,
                 conf.int = TRUE,#增加置信区间
                 conf.int.style="step", # 设置置信区间的类型，有"ribbon"(默认),"step"两种。
                 conf.int.alpha=0, # 数值，指定置信区间填充颜色的透明度；
                 # 数值在0-1之间，0为完全透明，1为不透明。
                 risk.table =TRUE,
                 combine=T,
                 pval =TRUE,
                 xlab ="Follow up time(d)",
                 surv.median.line = "hv",# 增加中位生存时间
                 legend = c(0.75,0.86), # 指定图例位置
                 legend.labs = c("High", "Low"), # 指定图例分组标签
                 legend.title="MRPL12 in LIHC Exp",  # 图例标题
                 risk.table.title="" ,#  风险表的标题
                 ggtheme = theme_bw(),
                 surv.plot.height=1)

png(filename = "LIHC-survival.png",width = 1920,height = 1920,units = "px",res = 400)
Ps
dev.off()
