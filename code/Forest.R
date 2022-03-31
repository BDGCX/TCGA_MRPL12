####Forest森林图####
###library(survival)
##my.surv <- Surv(OS_MONTHS,OS_STATUS=='DECEASED')这个生存对象是看看病人的总生存期
##与死亡状态的关系,这个Surv函数第一个参数必须是数值型的时间，第二个参数是逻辑向量，1,0表示死亡与否
##kmfit1 <- survfit(my.surv~1)直接对生存对象拟合生存函数
##summary(kmfit1)这个函数会详细打印出所有的结果，值得仔细理解
##plot(kmfit1)画出生存曲线
##survdiff(my.surv~type, data=dat)### 根据生存对象再加上一个分组因子来拟合生存函数，并且比较不同因子分组的生存效果
###用my.surv <- surv(OS_MONTHS,OS_STATUS=='DECEASED')构建生存曲线。
###用kmfit2 <- survfit(my.surv~TUMOR_STAGE_2009)来做某一个因子的KM生存曲线。用 survdiff(my.surv~type, data=dat)来看看这个因子的不同水平是否有显著差异，其中默认用是的logrank test 方法。
###用coxph(Surv(time, status) ~ ph.ecog + tt(age), data=lung) 来检测自己感兴趣的因子是否受其它因子(age,gender等等)的影响。
###风险比(HR)-hazard ratio（输出里面，coef就是beta值，相应的exp(coef)就是HR了
###在oncolic下载各种肿瘤的生存期数据整合到一个csv里
survival <- read.csv(file="BRCA.csv",header = T)
library(survival)
my.surv <- Surv(survival$Days,survival$Status=="Dead")
##survfit(my.surv~1)
##kmfit <- survfit(my.surv~1)
##plot(kmfit)
##summary(kmfit)
survdiff(my.surv~Group,data = survival)
coxmodel <- coxph(my.surv~Group+Patient,data = survival)
summary(coxmodel)###分别计算每种肿瘤的HR和P值填入表中
##Concordance在0.5-1之间。0.5为完全不一致,说明该模型没有预测作用
# 1为完全一致,说明该模型预测结果与实际完全一致。
##1.使用surminer包中的ggforest画图，适用于做cox风险回归后
library(ggplot2)
library(ggpubr)
library(survminer)
ggforest(coxmodel, data =survival, 
         main = "Hazard ratio", ##设置标题
         cpositions = c(0.10, 0.22, 0.4), ###设置前三列的相对距离
         fontsize = 1.0, 
         refLabel = "1", noDigits = 4)
###2.使用forestplot包使用临床生信之家下载的数据forest，适用于直接从网站下载得到HR等信息
library(grid)
library(magrittr)
library(checkmate)
library(forestplot)
forest <- read.csv(file="forestplot.csv",header = T)
forest$Type
A <- c("",forest$Type)
###options(digits = 3)
B <- paste0(round(forest$HR,2),"(",round(forest$Lower,2),"~",round(forest$Upper,2),")")
B <- c("HR",B)
C <- c("P Value",format(forest$p.value,scientific = T))
tabletext <- cbind(A,B,C)
colnames(tabletext) <- NULL
b <- as.numeric(forest$Lower)
c <- as.numeric(forest$Upper)
m <- (b+c)/2
mean(as.numeric(forest$Lower))
mean(as.numeric(forest$Upper))
cochrane_from_rmeta <- data.frame(
  mean  = c(NA, m), 
  lower = c(NA, b),
  upper = c(NA ,c))
forestplot(tabletext, 
           mean = cochrane_from_rmeta$mean,
           lower = cochrane_from_rmeta$lower ,
           upper = cochrane_from_rmeta$upper)

Pf <- forestplot(tabletext, 
                 cochrane_from_rmeta,
                 # 添加水平线
                 align = "c", # 设置左边表格中字体的对齐方式
                 zero = 1, # 设置zero line的位置
                 new_page = TRUE,
                 clip=c(0.2,2.5), #Lower and upper limits for clipping confidence intervals to arrows
                 xlog=TRUE,
                 xticks.digits = 2,
                 fn.ci_norm = fpDrawCircleCI, # box的样式，默认为方块fpDrawCircleCI，summary默认为菱形fpDrawDiamondCI.可选项：fpDrawNormalCI、fpDrawDiamondCI、fpDrawCircleCI、fpDrawPointCI、fpDrawSummaryCI、fpDrawBarCI
                 lty.ci = 2, # HR线（穿过box的直线）的线型，默认为1（实直线）
                 col=fpColors(box="red",line="darkblue", 
                              summary="royalblue", hrz_lines = "#444444"),
                 vertices = TRUE,
                 boxsize=0.5,
                 txt_gp=fpTxtGp(label=gpar(cex=0.6), 
                                ticks=gpar(cex=0.6),
                                xlab=gpar(cex = 0.6), 
                                title=gpar(cex = 0.8)))
png(filename = "forest.png",width = 1920,height = 1480,units = "px",res = 400)
Pf
dev.off()

###forestplot参数####
### labeltext | 主要是以矩阵或者list形式将数据导入函数，最好以矩阵，因为数据一般都是矩阵的。 |
#  | mean | 误差条的均值 |
#  | upper | 误差条 95%置信区间上限 |
#  | align | 每列文字的对齐方式，偶尔会用到。如：align=c("l","c","c")l：左对齐r：右对齐c：居中对齐 |
# | is.summary | 主要的功能是让表格的每一行字体出现差异，从而区分表头。其值主要用TRUE/FALSE进行差异化分配。 |
#  | graph.pos | 定位森林图所在的位置。通过数字来确定为第几列。 |
#  | hrzl_lines | 以list形式设置表中线条的类型、影响范围。Eg：“3”=gpar(lwd=1,columns=1:4,col=’red’)意思就是第3行的线条，宽度为1，线段延伸至第四列。Col指的颜色。 |
#  | clip | x轴的最大最小范围 |
#  | xlab | x轴的标题 |
#  | zero | 森林图中基准线的位置（无效线的横坐标） |
#  | graphwidth | 森林图在表中的宽度如：graphwidth = unit(.4,"npc") |
#  | colgap |列与列之间的间隙宽度，默认是 6 mm，需要用 unit 的形式
#| lineheight | 行的高度，可以是数字，也可以是 unit 的形式 |
#  | line.margin | 行与行之间的间隙的宽度 |
#  | col | 森林图横线以及点的颜色。box：box（点估计值）的颜色line：穿过方块的横线的颜色zero：中间那条基准线的颜色summary：summary中菱形的颜色hrz_lines：表中第一条横线的颜色eg：col=fpcolors(box=’royblue’,line=’darkblue’, summary=’royblue’, hrz_lines=’red’) |
#  | txt_gp | 设置表格中文本的格式：用gpar进行赋值，其中cex为文本字体大小，ticks为坐标轴大小，xlab为坐标轴文字字体大小。label：表格主体文字的格式ticks：森林图下方的坐标轴的刻度文字格式xlab：定义的x轴标题格式title：标题文字的格式eg：txt_gp=fpTxtGp(label=gpar(cex=1.25), ticks=gpar(cex=1.1), xlab=gpar(cex = 1.2),title=gpar(cex = 1.2)) |
#  | xticks | 横坐标刻度根据需要可随意设置，如：xticks = c(0.5, 1,1.5, 2) |
#    | lwd.xaxis | X轴线宽 |
#    | lwd.zero | 无效线的宽度 |
#    | lwd.ci | 置信区间线条的宽度（粗细） |
#    | lty.ci | 置信区间的线条类型 |
#    | ci.vertices | 森林图可信区间两端添加小竖线（TRUE） |
#    | ci.vertices.height | 设置森林图可信区间两端的小竖线高度，默认是10%行高 |
#    | boxsize | box（点估计值）的大小 |
#    | mar | 图形页边距，如：mar=unit(rep(1.25, times = 4), "cm") |
#      | title | 添加标题 |
#      | legend | 当同时显示多个置信区间时，需要添加图例 |
#    new_page| 是否新页 |
#  | fn.ci_norm | box（点估计值）的形状，默认是方块。如：fn.ci_norm="fpDrawDiamondCI"：box 类型选择钻石 |
