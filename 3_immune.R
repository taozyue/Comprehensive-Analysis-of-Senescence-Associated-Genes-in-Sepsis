
# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

depens<-c('tibble', 'survival', 'survminer', 'sva', 'limma', "DESeq2","devtools",
          'limSolve', 'GSVA', 'e1071', 'preprocessCore', 'ggplot2', "biomaRt",
          'ggpubr', "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor", "timeROC","pracma")
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))
    BiocManager::install(depen,update = FALSE)
}


#devtools::install_github("IOBR/IOBR",ref="master")
#推荐本地安装
#devtools::install_local('IOBR-master.zip')
library("IOBR")
load('key_train_exprSet.Rdata')
library(IOBR)
library(limma)
library(tidyr)
mypal <- c("#db6968","#4d97cd","#99cbeb","#459943",
           "#fdc58f","#e8c559","#a3d393","#f8984e","#93cc82","#4d97cd","#f6f5ee","#ea9c9d",
           "#c74546","#88c4e8","#606f8a","#e8c559","#ea9c9d","#005496","#88c4e8","#db6968","#982b2b","#0074b3",
           "#e5ce81","#f47720","#459943","#bdc3d2")
rt=as.matrix(rt)
exp=rt
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=normalizeBetweenArrays(mat)
mat=mat[rowMeans(mat)>0,]
cibersort<-deconvo_tme(mat, method = "cibersort",perm = 1000 )
write.table(cibersort,"CIBERSORT.txt",sep="\t",quote = F,row.names = F)
cibersort$ID=seq(1,ncol(mat),1)
res<-cell_bar_plot(input = cibersort, title = "CIBERSORT Cell Fraction",coord_filp = F)

pdf(file="CIBERSORT.pdf",  width=10, height=7)
print(res)
dev.off()

ciber=IOBR::deconvo_cibersort(eset=rt,arrays = T,perm = 100)
identical(ciber$ID,rownames(anno))

rownames(ciber)=ciber$ID
ciber=ciber[,-1]
ciber=ciber[,-(23:25)]

identical(rownames(ciber),rownames(anno))

sig_sense=read.csv('diffSig_sence.csv')
sig_sense=sig_sense$X

rt=rt[rownames(rt) %in% sig_sense,]
rt=as.data.frame(t(rt))

#install.packages("psych")
library(psych)
pp <- corr.test(rt,ciber,method="spearman",adjust = "fdr")

cor <- pp$r
pvalue <- pp$p

# 瞄一眼，有NA的列去除,这些列全是0!!!!!!!!!!!!!!!!!!!!!
colnames(ciber)
ciber=ciber[, ! colnames(ciber) %in% c("NK_cells_activated_CIBERSORT" , 
                                       "Dendritic_cells_resting_CIBERSORT" )]

library(psych)
library(reshape2)
pp <- corr.test(rt,ciber,method="spearman",adjust = "fdr")

cor <- pp$r
pvalue <- pp$p

myfun <- function(pval) {
  stars = ""
  if(pval <= 0.001)
    stars = "***"
  if(pval > 0.001 & pval <= 0.01)
    stars = "**"
  if(pval > 0.01 & pval <= 0.05)
    stars = "*"
  if(pval > 0.05 & pval <= 0.1)
    stars = ""
  stars
}


heatmap <- melt(cor)
colnames(heatmap)=c('sample','gene','cor')  

heatmap=mutate(heatmap,pvalue=melt(pvalue)[,3]) %>%
  mutate(signif = sapply(pvalue, function(x) myfun(x)))


ggplot(heatmap,aes(sample,gene,col=cor))+
  geom_tile(color="grey70",fill="white",size=1)+
  geom_point(aes(size = abs(cor)),shape=15) + 
  geom_text(aes(label=signif),size=6,color="white",
            hjust=0.5,vjust=0.7)+
  labs(x = NULL,y = NULL,color=NULL) + 
  scale_color_viridis_c()+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(text=element_text(),
        axis.ticks.x = element_blank(),axis.text.x = element_text( angle =45),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(fill=NA,color="grey70",
                                    size=2, linetype="solid")) +
  scale_size(range=c(1,10),guide=NULL)+
  guides(color = guide_colorbar(direction = "vertical",
                                reverse = F,barwidth = unit(.5, "cm"),
                                barheight = unit(15, "cm")))

#### 和诊断基因

load('gene_rfe_rf.Rdata')

rt=rt[,colnames(rt) %in% gene]

library(psych)
pp <- corr.test(rt,ciber,method="spearman",adjust = "fdr")

cor <- pp$r
pvalue <- pp$p

# 瞄一眼，有NA的列去除,这些列全是0
colnames(ciber)
ciber=ciber[, ! colnames(ciber) %in% c("Plasma_cells_CIBERSORT" , "T_cells_follicular_helper_CIBERSORT",
                                       "Dendritic_cells_resting_CIBERSORT")]

library(psych)
pp <- corr.test(rt,ciber,method="spearman",adjust = "fdr")

cor <- pp$r
pvalue <- pp$p

myfun <- function(pval) {
  stars = ""
  if(pval <= 0.001)
    stars = "***"
  if(pval > 0.001 & pval <= 0.01)
    stars = "**"
  if(pval > 0.01 & pval <= 0.05)
    stars = "*"
  if(pval > 0.05 & pval <= 0.1)
    stars = ""
  stars
}



heatmap <- melt(cor)
colnames(heatmap)=c('sample','gene','cor')  

heatmap=mutate(heatmap,pvalue=melt(pvalue)[,3]) %>%
  mutate(signif = sapply(pvalue, function(x) myfun(x)))


ggplot(heatmap,aes(sample,gene,col=cor))+
  geom_tile(color="grey70",fill="white",size=1)+
  geom_point(aes(size = abs(cor)),shape=15) + 
  geom_text(aes(label=signif),size=6,color="white",
            hjust=0.5,vjust=0.7)+
  labs(x = NULL,y = NULL,color=NULL) + 
  scale_color_viridis_c()+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(text=element_text(),
        axis.ticks.x = element_blank(),axis.text.x = element_text( angle =45),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(fill=NA,color="grey70",
                                    size=2, linetype="solid")) +
  scale_size(range=c(1,10),guide=NULL)+
  guides(color = guide_colorbar(direction = "vertical",
                                reverse = F,barwidth = unit(.5, "cm"),
                                barheight = unit(15, "cm")))



######################mcpcounter(绝对)
load('key_train_exprSet.Rdata')

library(IOBR)
mcp=IOBR::deconvo_mcpcounter(eset=rt)
rownames(mcp)=mcp$ID
mcp=mcp[,-1]

identical(rownames(mcp),rownames(anno))


sig_sense=read.csv('diffSig_sence.csv')
sig_sense=sig_sense$X


rt=rt[rownames(rt) %in% sig_sense,]
rt=as.data.frame(t(rt))

colnames(mcp)

mcp=mcp[,! colnames(mcp) %in% c("Endothelial_cells_MCPcounter","Fibroblasts_MCPcounter" )]
library(psych)
pp <- corr.test(rt,mcp,method="spearman",adjust = "fdr")

cor <- pp$r
pvalue <- pp$p

myfun <- function(pval) {
  stars = ""
  if(pval <= 0.001)
    stars = "***"
  if(pval > 0.001 & pval <= 0.01)
    stars = "**"
  if(pval > 0.01 & pval <= 0.05)
    stars = "*"
  if(pval > 0.05 & pval <= 0.1)
    stars = ""
  stars
}


heatmap <- melt(cor)
colnames(heatmap)=c('sample','gene','cor')  

heatmap=mutate(heatmap,pvalue=melt(pvalue)[,3]) %>%
  mutate(signif = sapply(pvalue, function(x) myfun(x)))


ggplot(heatmap,aes(sample,gene,col=cor))+
  geom_tile(color="grey70",fill="white",size=1)+
  geom_point(aes(size = abs(cor)),shape=15) + 
  geom_text(aes(label=signif),size=6,color="white",
            hjust=0.5,vjust=0.7)+
  labs(x = NULL,y = NULL,color=NULL) + 
  scale_color_viridis_c()+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(text=element_text(),
        axis.ticks.x = element_blank(),axis.text.x = element_text( angle =45),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(fill=NA,color="grey70",
                                    size=2, linetype="solid")) +
  scale_size(range=c(1,10),guide=NULL)+
  guides(color = guide_colorbar(direction = "vertical",
                                reverse = F,barwidth = unit(.5, "cm"),
                                barheight = unit(15, "cm")))

#### 和诊断基因

load('gene_rfe_rf.Rdata')

rt=rt[,colnames(rt) %in% gene]

library(psych)
pp <- corr.test(rt,mcp,method="spearman",adjust = "fdr")

cor <- pp$r
pvalue <- pp$p


myfun <- function(pval) {
  stars = ""
  if(pval <= 0.001)
    stars = "***"
  if(pval > 0.001 & pval <= 0.01)
    stars = "**"
  if(pval > 0.01 & pval <= 0.05)
    stars = "*"
  if(pval > 0.05 & pval <= 0.1)
    stars = ""
  stars
}


heatmap <- melt(cor)
colnames(heatmap)=c('sample','gene','cor')  

heatmap=mutate(heatmap,pvalue=melt(pvalue)[,3]) %>%
  mutate(signif = sapply(pvalue, function(x) myfun(x)))


ggplot(heatmap,aes(sample,gene,col=cor))+
  geom_tile(color="grey70",fill="white",size=1)+
  geom_point(aes(size = abs(cor)),shape=15) + 
  geom_text(aes(label=signif),size=6,color="white",
            hjust=0.5,vjust=0.7)+
  labs(x = NULL,y = NULL,color=NULL) + 
  scale_color_viridis_c()+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(text=element_text(),
        axis.ticks.x = element_blank(),axis.text.x = element_text( angle =45),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(fill=NA,color="grey70",
                                    size=2, linetype="solid")) +
  scale_size(range=c(1,10),guide=NULL)+
  guides(color = guide_colorbar(direction = "vertical",
                                reverse = F,barwidth = unit(.5, "cm"),
                                barheight = unit(15, "cm")))



library(ggpubr)
library(ggstatsplot)
library(ggside)
data=cbind(rt,mcp)
gene

ggscatterstats(data,x = "ETS2",y="Neutrophils_MCPcounter")

p <- ggplot(data = data, aes(x = ETS2, y = Neutrophils_MCPcounter)) +
  geom_point() + geom_smooth(method = lm) +
  labs(title = 'EST2 and neutrophil correlation') +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  stat_cor(method = 'pearson', aes(x = ETS2, y = Neutrophils_MCPcounter))
p


p1 <- ggplot(data = data, aes(x = MMP9, y = Neutrophils_MCPcounter)) +
  geom_point() + geom_smooth(method = lm) +
  labs(title = 'MMP9 and neutrophil correlation') +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  stat_cor(method = 'pearson', aes(x = MMP9, y = Neutrophils_MCPcounter))
p1
##利用该包进行绘制相关性的散点图和条形图
set.seed(123)


####炎症因子
load('key_train_exprSet.Rdata')
gene=read.table('BIOCARTA_INFLAM_PATHWAY.v7.5.1.gmt')[3:31]
gene=as.data.frame(t(gene))
gene=gene$V1
rt1=rt[rownames(rt)%in% gene,]
identical(colnames(rt1),rownames(anno))

load('gene_rfe_rf.Rdata')

rt2=rt[rownames(rt)%in% gene,]
identical(colnames(rt2),rownames(anno))

rt1=as.data.frame(t(rt1))
rt2=as.data.frame(t(rt2))
library(psych)
pp <- corr.test(rt2,rt1,method="spearman",adjust = "fdr")

cor <- pp$r
pvalue <- pp$p


myfun <- function(pval) {
  stars = ""
  if(pval <= 0.001)
    stars = "***"
  if(pval > 0.001 & pval <= 0.01)
    stars = "**"
  if(pval > 0.01 & pval <= 0.05)
    stars = "*"
  if(pval > 0.05 & pval <= 0.1)
    stars = ""
  stars
}

library(tidyverse)

library(psych)
library(reshape2)
heatmap <- melt(cor)
colnames(heatmap)=c('sample','gene','cor')  

heatmap=mutate(heatmap,pvalue=melt(pvalue)[,3]) %>%
 mutate(signif = sapply(pvalue, function(x) myfun(x)))

library(ggsci)
ggplot(heatmap,aes(sample,gene,col=cor))+
  geom_tile(color="grey70",fill="white",size=1)+
  geom_point(aes(size = abs(cor)),shape=15) + 
  geom_text(aes(label=signif),size=6,color="white",
            hjust=0.5,vjust=0.7)+
  labs(x = NULL,y = NULL,color=NULL) + 
  scale_color_viridis_c()+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(text=element_text(),
        axis.ticks.x = element_blank(),axis.text.x = element_text( angle =45),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(fill=NA,color="grey70",
                                    size=2, linetype="solid")) +
  scale_size(range=c(1,10),guide=NULL)+
  guides(color = guide_colorbar(direction = "vertical",
                                reverse = F,barwidth = unit(.5, "cm"),
                                barheight = unit(15, "cm")))
