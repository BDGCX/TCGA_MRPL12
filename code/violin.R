rm(list=ls())
options(stringsAsFactors = F)
##加载需要的R包
library(ggplot2)
library(dplyr)
library(ggpubr)
library(hrbrthemes)
library(viridis)
library(RColorBrewer)
library(ggsci)
###1.取颜色方法RColorBrewer####
##myPalette <- colorRampPalette(brewer.pal(8,"Set1"))(26)
###2.取颜色方法RColorBrewer
##qual_col_pals <-  brewer.pal.info[brewer.pal.info$category == 'qual',]
#处理后有73种差异还比较明显的颜色，基本够用
##col_vector <-  unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#看下中间60种颜色的效果
###取颜色ggsci包（主要用这个）######
colpalettes<-unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))
###CCLE数据库-Dataset-CCLE data搜索-expression...下载数据
CCLE <- read.csv(file ="./MRPL12 Expression 22Q1 Public.csv",header = T)
ccle <- as.data.frame(cbind(Lineage=CCLE$Lineage,Expression.22Q1.Public=CCLE$Expression.22Q1.Public))
#####小提琴合并箱线图
sample_size = ccle %>% group_by(Lineage) %>% dplyr::summarize(num=n())
###将数量较少的样本去除,根据sample_size搜索行号
ccle <- ccle[-c(188,1123,1215,1254,1144,168,489,47,62,636,1030,1072,1379:1385),]
ccle$Expression.22Q1.Public <- as.numeric(ccle$Expression.22Q1.Public)
sample_size = ccle %>% group_by(Lineage) %>% dplyr::summarize(num=n())
######Kruskal-Wallis检验
kruskal.test(Expression.22Q1.Public~Lineage,data = ccle)

####"\n"相当于回车

####小提琴合并箱线图####
P1 <- ccle %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(Lineage," ","(N=",num, ")")) %>%
  ggplot( aes(x=myaxis, y=Expression.22Q1.Public)) +
  geom_violin(width=1.2) +
  geom_boxplot(aes(fill =Lineage),width=0.2, color="black", alpha=0.2,outlier.colour=NA)+#填充颜色，箱宽，不显示离群点
  scale_fill_manual(values= colpalettes)+##ggsci包###
  theme(legend.position="none",
      plot.title = element_text(size=14),
      panel.grid.major=element_line(colour="grey"),
      panel.grid.minor=element_line(colour="grey"),
      panel.background = element_rect(fill ="white", colour = "grey"),
      axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Italic",size=6))+ #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)字体簇为Times大小为20
  ggtitle("Kruskal-Wallis test p=6.461e-16") +
  xlab("")+ylab("Gene Expression")

png(filename = "ccle.png",width = 1920,height = 1480,units = "px",res = 400)
P1
dev.off()
####GTEx（UCSC Xena中下载GTEx数据）####
GTEx <- read.table(file = "denseDataOnlyDownload.tsv",sep = "\t",header = T)
View(GTEx)
gtex <- GTEx[,c(4,3)]
View(gtex)
colnames(gtex) <- c("Lineage","Values")
sample_size1 = gtex %>% group_by(Lineage) %>% dplyr::summarize(num=n())
gtex <- gtex[-c(758,1437,2730,4605,7895),]
sample_size1 = gtex %>% group_by(Lineage) %>% dplyr::summarize(num=n())
gtex1 <- filter(gtex,Values>0)###去除负值的样本
kruskal.test(Values~Lineage,data = gtex1)
####小提琴合并箱线图####
P2 <- gtex1 %>%
  left_join(sample_size1) %>%
  mutate(myaxis = paste0(Lineage," ","(N=",num, ")")) %>%
  ggplot( aes(x=myaxis, y=Values)) +
  geom_violin(width=1.8) +
  geom_boxplot(aes(fill =Lineage),width=0.2, color="black", alpha=0.2,outlier.colour=NA)+#填充颜色，箱宽，不显示离群点
  scale_fill_manual(values =colpalettes )+
  theme(legend.position="none",
        plot.title = element_text(size=14),
        panel.grid.major=element_line(colour="grey"),
        panel.grid.minor=element_line(colour="grey"),
        panel.background = element_rect(fill ="white", colour = "grey"),
        axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Italic",size=6),
        axis.title.y=element_text(size = 10))+ #设置x轴刻度标签的字体显示倾斜角度为45度，并向下调整1(hjust = 1)字体簇为Times大小为20
  ggtitle("Kruskal-Wallis test p=2.2e-16") +
  xlab("")+ylab("log2(TPM+1)")
####存储为高清图片
png(filename = "GTEx.png",width = 1980,height = 1280,units = "px",res = 400)
P2
dev.off()##即在此路径下产生一个高清图片

####UCSC Xena数据库下载 TCGA TARGET GTEx 数据绘制分半小提琴图,此路太难放弃了####
TG <- read.table(file = "TCGA-GTEx.tsv",sep = "\t",header = T)
View(TG)
TG <- TG[,c(7,5,3,4)]
colnames(TG) <- c("Lineage","Group","values","primary site")
TG <- TG[-c(6457,11606,12598,15589),]##去除组织部位未知的行
TG$Group <- gsub("Solid Tissue Normal","Normal Tissue",TG$Group)##替换组别名字
TG$`primary site`<- gsub("Brain","GBM",TG$`primary site`)##替换组织名字为TCGA肿瘤名称

#####临床生信之家下载泛癌TCGA+GTEx数据绘制分半小提琴图####
tg <- read.csv(file = "MRPL12-TCGA+GTEx.csv",header = T)
colnames(tg) <- c("X","Lineage","Group","values")
tg$Lineage <- substring(tg$Lineage,1,4)
X <-tg$Lineage 
SUB <- sub("[:(:]","\\",X)###用法见https://www.jb51.net/article/207342.htm
tg <- cbind(Lineage=SUB,tg)
tg <- tg[,-3]
SUB <- sub("OVT","OV",tg$Lineage)
tg <- cbind(Lineage=SUB,tg)
tg <- tg[,-2]
View(tg)
tg1 <- tg[,-2]
View(tg1)
###删除不想要的行
w <- which(tg1$Lineage=="UVM(T=80;N=0)")
TG <- tg1[-w,]
w <- which(TG$Lineage=="MESO(T=86;N=0)")
TG <- TG[-w,]
sample_size2 = TG %>% group_by(Lineage,Group) %>% dplyr::summarize(num=n())
###计算一下误差信息####
library(Rmisc)
Data_summary <- summarySE(TG, measurevar="values", groupvars=c("Lineage","Group"))
head(Data_summary)
####包装函数geom_split_violin####
GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, 
                      xminv = x - violinwidth * (x - xmin), 
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)
geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}
#####绘制分半小提琴图#####
myaxis <- c("ACC(N=258,T=79)","BLCA(N=21,T=408)","BRCA(N=459,T=1097)",
            "CESC(N=19,T=306)","CHOL(N=9,T=36)","COAD(N=779,T=455)",
            "DLBC(N=929,T=48)",
            "ESCA(N=1445,T=162)","GBM(N=2642,T=153)","HNSC(N=44,T=502)",
            "KICH(N=89,T=65)","KIRC(N=89,T=530)","KIRP(N=89,T=288)",
            "LGG(N=2642,T=510)",
            "LICH(N=226,T=371)","LUAD(N=578,T=513)",
            "LUSC(N=578,T=501)","OV(N=180,T=374)","PAAD(N=328,T=178)",
            "PRAD(N=245,T=496)","READ(N=779,T=165)","SKCM(N=1809,T=470)",
            "STAD(N=359,T=375)","TGCT(N=361,T=134)","THCA(N=653,T=510)",
            "UCEC(N=142,T=543)","UCSN=(142,T=56)")
P3 <- ggplot(data=TG,aes(x=Lineage, y=values,fill=Group)) + 
  geom_split_violin(trim=F,color="white",scale = "width") + #绘制分半的小提琴图
  geom_point(data = Data_summary,aes(x=Lineage, y=values),pch=19,position=position_dodge(0.3),size=0.5)+ #绘制均值为点图
  geom_errorbar(data = Data_summary,aes(ymin = values-ci, ymax=values+ci), #误差条表示95%的置信区间
                width=0.1, #误差条末端短横线的宽度
                position=position_dodge(0.3), 
                color="black",
                alpha = 0.7,
                size=0.5) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ #设置填充的颜色
  theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=5), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=6,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 9,face="plain"), #设置y轴标题的字体属性
        legend.text=element_text(face="plain", family="Arial", colour="black",  #设置图例的子标题的字体属性
                                 size=6),
        legend.title=element_text(face="plain", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=6),
        legend.position = c(0.1,0.1),##调整图例位置
        legend.background = element_blank(),##去掉图例背景
        ###legend.key = element_blank(),##去掉图例背景
        ##legend.background = element_rect(fill ="",color="",size= ),##调整图例参数
        legend.key.size = unit(10, "pt"),
        # 或legend.key.height = unit(35, "pt"),legend.key.width = unit(55, "pt")
        panel.grid.major=element_line(colour="grey"),
        panel.grid.minor=element_line(colour="grey"),
        panel.background = element_rect(fill ="white", colour = "grey"))+  #不显示网格线
  ylab("log2(TPM+1)")+xlab("")+ #设置x轴和y轴的标题
  scale_x_discrete(labels= myaxis)+###自定义X轴刻度标签
  stat_compare_means(aes(group = Group),
                     label = "p.signif",##显示显著性水平,p.format显示四舍五入的p-value，p.adj为调整后的p-value
                     method = "anova",##选择统计检验方法
                     label.y = max(TG$values),###调整标签位置
                     size=2,
                     symnum.args=list(cutpoints = c(0,0.001, 0.01, 0.05, 1), symbols = c( "***", "**", "*", "ns")),
                     hide.ns = T)###T为隐藏不显著标记(ggpubr包)

png(filename = "TG.png",width = 2000,height = 1280,units = "px",res = 400)
P3
dev.off()
###单纯TCGA泛癌分析(TCGA里的正常太少了）####
a <- grep("GTEX",tg$X)###在长字符串种模糊查找某短字符串，返回元素下标，grepl()则返回逻辑值
tcga <- tg[-a,]
tcga <- tcga[,-2]
colnames(tcga) <- c("Lineage","Group","values")
v <- which(tcga$Lineage=="SKCM")
tcga<- tcga[-v,]
Data_summary1 <- summarySE(tcga, measurevar="values", groupvars=c("Lineage","Group"))
myaxis1 <- c("BLCA(N=19,T=408)","BRCA(N=113,T=1097)",
             "CESC(N=3,T=306)","CHOL(N=9,T=36)","COAD(N=41,T=455)",
             "ESCA(N=11,T=162)","GBM(N=5,T=153)","HNSC(N=44,T=502)",
             "KICH(N=24,T=65)","KIRC(N=72,T=530)","KIRP(N=32,T=288)",
             "LICH(N=50,T=371)","LUAD(N=59,T=513)",
             "LUSC(N=49,T=501)","PAAD(N=4,T=178)",
             "PRAD(N=52,T=496)","READ(N=10,T=165)",
             "STAD(N=32,T=375)","THCA(N=58,T=510)",
             "UCEC(N=35,T=543)")
P4 <- ggplot(data=tcga,aes(x=Lineage, y=values,fill=Group)) + 
  geom_split_violin(trim=F,color="white",scale = "width") + #绘制分半的小提琴图
  geom_point(data = Data_summary1,aes(x=Lineage, y=values),pch=19,position=position_dodge(0.3),size=0.5)+ #绘制均值为点图
  geom_errorbar(data = Data_summary1,aes(ymin = values-ci, ymax=values+ci), #误差条表示95%的置信区间
                width=0.1, #误差条末端短横线的宽度
                position=position_dodge(0.3), 
                color="black",
                alpha = 0.7,
                size=0.5) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ #设置填充的颜色
  theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=5), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=6,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 9,face="plain"), #设置y轴标题的字体属性
        legend.text=element_text(face="plain", family="Arial", colour="black",  #设置图例的子标题的字体属性
                                 size=6),
        legend.title=element_text(face="plain", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=6),
        legend.position = c(0.1,0.1),##调整图例位置
        legend.background = element_blank(),##去掉图例背景
        ###legend.key = element_blank(),##去掉图例背景
        ##legend.background = element_rect(fill ="",color="",size= ),##调整图例参数
        legend.key.size = unit(10, "pt"),
        # 或legend.key.height = unit(35, "pt"),legend.key.width = unit(55, "pt")
        panel.grid.major=element_line(colour="grey"),
        panel.grid.minor=element_line(colour="grey"),
        panel.background = element_rect(fill ="white", colour = "grey"))+  #不显示网格线
  ylab("log2(TPM+1)")+xlab("")+ #设置x轴和y轴的标题
  scale_x_discrete(labels= myaxis1)+###自定义X轴刻度标签
  stat_compare_means(aes(group = Group),
                     label = "p.signif",##显示显著性水平,p.format显示四舍五入的p-value，p.adj为调整后的p-value
                     method = "anova",##选择统计检验方法
                     label.y = max(TG$values),###调整标签位置
                     size=2,
                     symnum.args=list(cutpoints = c(0,0.001, 0.01, 0.05, 1), symbols = c( "***", "**", "*", "ns")),
                     hide.ns = T)###T为隐藏不显著标记(ggpubr包)
png(filename = "TCGA.png",width = 2000,height = 1280,units = "px",res = 400)
P4
dev.off()
