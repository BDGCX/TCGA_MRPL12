###相关性分析散点图
mrpl <- exprSet["MRPL12",]
mrpl <- t(mrpl)
cor <- cbind(MRPL12=mrpl,scores)
head(cor)
cor <- as.data.frame(cor)
cor1 <- cor[,c(2,1)]
###第1种方法
library(ggstatsplot)
library(ggside)
ggscatterstats(data = cor1,                                          
  x = StromalScore,                                                  
  y = MRPL12,
  bf.message = FALSE,#去除贝叶斯相关统计值
  results.subtitle = FALSE,##去除统计结果小标题
  point.args = list(size = 2, alpha = 1, stroke = 0, na.rm = TRUE),
  xlab = "StromalScore",
  ylab = "log2(MRPL12 TPM+1)",
  marginal.type = "density", #类型可以换成density,boxplot,violin,densigram
  margins = "both",
  xfill = "#FFD39B",
  yfill = "#B0E2FF"
)
###2.cowplot包拼图(主要用这种)
library(cowplot)
library(ggpmisc)
corr <- round(cor(cor1$MRPL12,cor1$StromalScore),2)##计算相关系数
pmain <- ggplot(cor1, aes(x = StromalScore, y = MRPL12)) +
  ylab("log2(MRPL12 TPM+1)")+geom_point(size=1) +theme_bw()+
  geom_smooth(method="lm")+
  stat_poly_eq(aes(label=paste("italic(r)~`=`",corr,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=4)
##加文字annotate("text", x=4.5, y=4.25, parse=TRUE,label="r^2 == 0.0138 * ' p-value = 0.1529' ")
#label.x 和 label.y 的取值范围与X轴和Y轴的取值对应。
##也可以通过 label.x.npc 和 label.y.npc 来指定显示文本的位置，但 label.x.npc 和 label.y.npc的取值范围为0~1，如0.5表示居中
#stat_cor() 也可以自动添加相关系数和显著性水平(P值)。
# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = cor1, aes(x = StromalScore), color ="black",fill="#FFD39B" ,
               alpha = 0.7, size = 0.2)+theme_light()
# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(pmain, axis = "y") +
  geom_density(data = cor1, aes(y = MRPL12), color="black",fill = "#B0E2FF",
               alpha = 0.7, size = 0.2) +theme_light() 
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)
ggsave(p2,filename = 'correlaton-scatter.png')
####3.ggscatterhist()画边际图####
library(ggpubr)
library(ggsci)
ggscatterhist(cor1, x="StromalScore", y="MRPL12", group = NULL, color = "black",
              shape = 19, size = 1, linetype = "solid", bins = 30,
              margin.params = list(fill="#B0E2FF"),#适用于所有边际图
              margin.plot = c("density"), 
              margin.ggtheme = theme_light(), 
              margin.space = FALSE, main.plot.size = 2,
              margin.plot.size = 1, title = NULL, 
              xlab = "StromalScore", ylab = "log2(MRPL12 TPM+1)", 
              legend = "top",
              add.params = list(color = "blue", fill = "gray"),
              cor.coef=TRUE,
              add = "reg.line", conf.int = TRUE,
              ggtheme = theme_bw(), print = TRUE)
####4.ggExtra包画边际图####
library(ggplot2)
library(ggExtra)
sd <- ggplot(cor1, aes(x =StromalScore , y = MRPL12)) +
  # 散点图函数
  geom_point()+
  theme_bw() +
  labs(x = "StromalScore", y = "log2(MRPL12 TPM+1)")+# 添加x，y轴的名称
  geom_smooth()# 拟合回归线段以及置信域(默认0.95/通过level参数可自定义，method="lm"直线，loess为曲线)
library(ggpmisc)
corr=round(cor(cor1$MRPL12,cor1$StromalScore),2)##计算相关系数
sd <- sd+stat_poly_eq(aes(label=paste("italic(r)~`=`",corr,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=4)
ggMarginal(sd, type = "density",
           xparams = list(colour = "black", fill="#FFD39B",size = 0.1),
           yparams = list(colour = "black", fill="#B0E2FF",size = 0.1))

model1 <- lm(cor1$MRPL12 ~ cor1$StromalScore, data = cor1) #回归拟合
summary(model1) #回归分析表
