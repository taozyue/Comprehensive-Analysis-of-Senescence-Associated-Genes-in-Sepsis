load('key_train_exprSet.Rdata')
rt=rt[,anno$group=='sepsis']

sig_sense=read.csv('diffSig_sence.csv')
sig_sense=sig_sense$X

rt=rt[rownames(rt) %in% sig_sense,]

#聚类
library(ConsensusClusterPlus)
rt=as.matrix(rt)
maxK=10
results=ConsensusClusterPlus(rt,
                             maxK=maxK,
                             reps=100,
                             pItem=0.8,
                             pFeature=1,
                             title='final_cluster',
                             clusterAlg="km",
                             distance="euclidean",
                             seed=123456,
                             plot="pdf",)

#一致性打分
dir.create('calscore')
setwd('calscore/')
calcICL(results, title="consensusScore", plot="pdf")
#输出结果
clusterNum=2     #分几类，根据判断标准判断
cluster=results[[clusterNum]][["consensusClass"]]
#setwd('./stroke_mito/')
setwd('..')
write.table(cluster, file="cluster.txt", sep="\t", quote=F, col.names=F)



### 作热图
load('key_train_exprSet.Rdata')
rt=as.data.frame(rt)
sig_sense=read.csv('diffSig_sence.csv')
sig_sense=sig_sense$X

rt=rt[rownames(rt) %in% sig_sense,]

cluster=read.table('cluster.txt',row.names = 1)
cluster=cluster[order(cluster$V2),,drop=F]
colnames(cluster)='Cluster'
rt=rt[,rownames(cluster)]
cluster$Cluster=factor(cluster$Cluster)
pheatmap::pheatmap(rt,scale = 'row',cluster_cols = F,border_color = NA,
                   breaks = seq(-2,2,length.out = 100),annotation_col = cluster,
                   show_colnames = F)

###加临床
##

pdata=read.table("GSE57065/clinical.txt",row.names = 1,header = 1)
rownames(pdata)=paste0('GSE57065_',rownames(pdata))


colnames(pdata)

## 尽可能多的挑选临床信息!!!!!!!
anno=pdata[,c('age:ch1','gender:ch1')]
## 只选一列要deop=F!!!!!!!!!
#anno2=pdata2[,c("time after stroke (h):ch1"),drop=F]
anno=pdata
library(tidyverse)
cluster$sample=rownames(cluster)
anno$sample=rownames(anno)


anno=merge(cluster,anno,all=T)
anno=anno[!is.na(anno$Cluster),]
anno=anno[!is.na(anno$age),]
rownames(anno)=anno$sample
anno$sample=NULL
colnames(anno)

load('pred_logistic.Rdata')
pred_df=pred_df[rownames(anno),,drop=F]
anno$Probability=pred_df$pred


anno=anno[order(anno$Cluster),]
rt=rt[,rownames(anno)]
anno$age=as.numeric(anno$age)
#有临床信息的
pheatmap::pheatmap(rt,scale = 'row',cluster_cols = F,border_color = NA,
                   breaks = seq(-2,2,length.out = 100),annotation_col = anno,
                   show_colnames = F)

pdf(file="4.cluster年龄差异.pdf", width=8, height=6)
library(ggpubr)
ggboxplot(
 anno,
  x = "Cluster",
  y = "age",
  color = "black",
  fill = "Cluster",
  xlab = "Cluster",
  ylab = "age", palette=c("#20854e","#E18727",'#bc3c29')
)+stat_compare_means()
dev.off()

colnames(anno)

pdf(file="4.cluster诊断效能差异.pdf", width=8, height=6)
library(ggpubr)
ggboxplot(
  anno,
  x = "Cluster",
  y = "Probability",
  color = "black",
  fill = "Cluster",
  xlab = "Cluster",
  ylab = "Prob", palette=c("#20854e","#E18727",'#bc3c29')
)+stat_compare_means()
dev.off()


colnames(anno)

library(ggstatsplot)

ggstatsplot::ggbarstats(anno,x='gender',y='Cluster')

colnames(anno)

pdf(file="4.gender_cluster.pdf", width=8, height=6)
ggstatsplot::ggbarstats(anno,x='gender',y='Cluster')
dev.off()

pdf(file="4.state_cluster.pdf", width=8, height=6)
ggstatsplot::ggbarstats(anno,x='state',y='Cluster')
dev.off()

pdf(file="4.time_cluster.pdf", width=8, height=6)
ggstatsplot::ggbarstats(anno,x='time',y='Cluster')
dev.off()

### 箱线图
load('key_train_exprSet.Rdata')
rt=rt[,anno$group=='sepsis']

sig_sense=read.csv('diffSig_sence.csv')
sig_sense=sig_sense$X

rt=rt[rownames(rt) %in% sig_sense,]

cluster=read.table('cluster.txt',row.names = 1)
cluster=cluster[order(cluster$V2),,drop=F]
colnames(cluster)='Cluster'
rt=rt[,rownames(cluster)]
cluster$Cluster=factor(cluster$Cluster)

rt2=as.data.frame(t(rt))
rt2$group=cluster$Cluster

library(tidyverse)
rt2=tidyr::pivot_longer(rt2,cols = -c('group'),names_to = "gene",values_to = 'expression')


library(ggpubr)
pdf(file="4.cluster中核心衰老基因的差异.pdf", width=8, height=6)
ggboxplot(
  rt2,
  x = "gene",
  y = "expression",
  color = "black",
  fill = "group",
  xlab = "group",
  ylab = "expression", palette=c("#20854e","#E18727",'#bc3c29')
) +
  stat_compare_means(
    aes(group = group),
    label = "p.signif", 
    method = "wilcox.test",
    hide.ns = T,
    size = 4.5
  ) +
  theme(axis.text.x = element_text(
    angle =45,
    hjust = 1,
    vjust = 1
  ))
dev.off()

####MCP
library("IOBR")
load('key_train_exprSet.Rdata')

mcp=IOBR::deconvo_mcpcounter(eset=rt)
identical(mcp$ID,rownames(anno))

rownames(mcp)=mcp$ID
mcp=mcp[,-1]
# 外周血去除内皮和纤维
mcp=mcp[,-(9:10)]

identical(rownames(mcp),rownames(anno))

cluster=read.table('cluster.txt')
colnames(cluster)=c('Sample','Cluster')

mcp=mcp[cluster$Sample,]

mcp$group=cluster$Cluster


library(tidyverse)
rt2=tidyr::pivot_longer(mcp,cols = -c('group'),names_to = "cell",values_to = 'score')

rt2$group=factor(rt2$group)
library(ggpubr)

pdf(file="4.cluster中MCP免疫细胞的差异.pdf", width=8, height=6)
ggboxplot(
  rt2,
  x = "cell",
  y = "score",
  color = "black",
  fill = "group",
  xlab = "group",
  ylab = "score", palette=c("#20854e","#E18727",'#bc3c29')
) +
  stat_compare_means(
    aes(group = group),
    label = "p.signif", 
    method = "wilcox.test",
    hide.ns = T,
    size = 4.5
  ) +
  theme(axis.text.x = element_text(
    angle =45,
    hjust = 1,
    vjust = 1
  ))
dev.off()

####炎症因子
load('key_train_exprSet.Rdata')
gene=read.table('BIOCARTA_INFLAM_PATHWAY.v7.5.1.gmt')[3:31]
gene=as.data.frame(t(gene))
gene=gene$V1
rt1=rt[rownames(rt)%in% gene,]
identical(colnames(rt1),rownames(anno))


cluster=read.table('cluster.txt')
colnames(cluster)=c('Sample','Cluster')


rt1=as.data.frame(t(rt1))
rt1=rt1[cluster$Sample,]

rt1$group=cluster$Cluster


library(tidyverse)
rt2=tidyr::pivot_longer(rt1,cols = -c('group'),names_to = "cell",values_to = 'expression')

rt2$group=factor(rt2$group)
library(ggpubr)

pdf(file="4.cluster中炎症因子的差异.pdf", width=8, height=6)
ggboxplot(
  rt2,
  x = "cell",
  y = "expression",
  color = "black",
  fill = "group",
  xlab = "group",
  ylab = "expression", palette=c("#20854e","#E18727",'#bc3c29')
) +
  stat_compare_means(
    aes(group = group),
    label = "p.signif", 
    method = "wilcox.test",
    hide.ns = T,
    size = 4.5
  ) +
  theme(axis.text.x = element_text(
    angle =45,
    hjust = 1,
    vjust = 1
  ))
dev.off()