#GSE16561
##——————————————————————————数据下载和整理——————————————————————————————————
if(! require("GEOquery")) BiocManager::install("GEOquery",update=F,ask=F)

library(GEOquery)
##下载并且了解你的数据
gset=getGEO('GSE16561',destdir = '.',getGPL = F)
#
exprSet=exprs(gset[[1]])
pdata=pData(gset[[1]])

library(data.table)
exprSet=fread('GSE16561_RAW.txt',data.table = F)
rownames(exprSet)=exprSet$ID_REF
exprSet=exprSet[,-1]

#整理文件最重要！
#问题解决：如何选择需要的样本？如何生成group_list告诉系统？
#group_list=
#exprSet=

group_list=c(rep('IS',39),rep('CT',24))
group_list=factor(group_list,levels = c('CT','IS'))

##————————————标准化和log2转化——————————————————#
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
max(exprSet)
exprSet=log2(exprSet+1)
exprSet=as.data.frame(exprSet)

##——————————————————ID转换——————————————————————————————
#开始转换
#载入注释文件，在你自己的文件夹里
#没data.table,geoquery啥的先装一下
#install.package('data.table')
#BiocManager::install('GEOquery',update = F,ask = F)
library(data.table)
anno = fread("GPL6883-11606.txt",data.table = F,sep = '\t') 
#看一下列名
colnames(anno) 
#一定选ID探针列和gene symbol列，！！！此处是可变的！！！
anno = anno[,c(1,6)] 
#!!!!
anno = anno[!anno$ILMN_Gene== "",] 
#无脑跑代码即可
tmp = rownames(exprSet)%in% anno[,1] 
exprSet = exprSet[tmp,] 
dim(exprSet) 
match(rownames(exprSet),anno$ID) 
anno = anno[match(rownames(exprSet),anno$ID),] 
match(rownames(exprSet),anno$ID) 
dim(exprSet) 
dim(anno) 
tail(sort(table(anno[,2])), n = 12L) 
#注释到相同基因的探针，保留表达量最大的那个探针 
{ 
  MAX = by(exprSet, anno[,2],  
           function(x) rownames(x)[ which.max(rowMeans(x))]) 
  MAX = as.character(MAX) 
  exprSet = exprSet[rownames(exprSet) %in% MAX,] 
  rownames(exprSet) = anno[ match(rownames(exprSet), anno[,1] ),2] 
} 
dim(exprSet) 
colnames(exprSet)=rownames(pdata)
save(exprSet,pdata,file ='GSE16561.Rdata')
write.table(exprSet,file ='GSE16561.txt',sep = '\t',col.names = NA,quote = F)


design=model.matrix(~ group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1) 

library(readr)
gene=read.table('REACTOME_MITOPHAGY.v7.5.1.gmt')
gene=as.data.frame(t(gene))
gene=gene[3:31,]

alldiff_mito=allDiff[rownames(allDiff) %in% gene,]

rt=exprSet[rownames(alldiff_mito),]
pheatmap::pheatmap(rt,scale = 'row',cluster_cols = F,border_color = NA,
                   breaks = seq(-2,2,length.out = 100))



#GSE58294--------------------------
############
setwd('d:/work/stroke_mito/')
library(affy)
# cel文件存放的路径（下载所有cel文件后，不需要彻底解压）
dir_cels <- "GSE58294_RAW/"  
# 读取readaffy
affy <- ReadAffy(celfile.path = dir_cels)
# rma标准化（推荐）
gset <- rma(affy)
exprSet=exprs(gset)
max(exprSet)
# 不需要log2转化

library(stringr)
expr_col= str_split(colnames(exprSet),pattern = '_',simplify = T)[,1]

colnames(exprSet)=expr_col

## 完成
exprSet=as.data.frame(exprSet)


## ID转换
library(limma)
## 可以看到均值是有轻微的不均的
boxplot(exprSet,outline=FALSE, notch=T, las=2)
library(limma) 
## 分位数校正
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,las=2)
max(exprSet)
exprSet=as.data.frame(exprSet)


library(data.table)
##下载https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16043
anno = fread("d:/Rlearning/platform/GPL570-55999.txt",data.table = F,sep = '\t') 
#看一下列名
colnames(anno) 
#一定选ID探针列和gene symbol列，！！！此处是可变的！！！
anno = anno[,c(1,11)] 
#若有不同，记得修改！
anno = anno[!anno$`Gene Symbol`== "",] 
#跑代码即可
tmp = rownames(exprSet)%in% anno[,1] 
exprSet = exprSet[tmp,] 
dim(exprSet) 
match(rownames(exprSet),anno$ID) 
anno = anno[match(rownames(exprSet),anno$ID),] 
match(rownames(exprSet),anno$ID) 
dim(exprSet) 
dim(anno) 
tail(sort(table(anno[,2])), n = 12L) 
#注释到相同基因的探针，保留表达量最大的那个探针 
{ 
  MAX = by(exprSet, anno[,2],  
           function(x) rownames(x)[ which.max(rowMeans(x))]) 
  MAX = as.character(MAX) 
  exprSet = exprSet[rownames(exprSet) %in% MAX,] 
  rownames(exprSet) = anno[ match(rownames(exprSet), anno[,1] ),2] 
} 
# 完美解决
dim(exprSet) 

#顺便再看看分组
library(GEOquery)
gset=getGEO('GSE58294',destdir = '.',getGPL = F)
pdata=pData(gset[[1]])

identical(colnames(exprSet),rownames(pdata))

save(exprSet,pdata,file ='GSE58294.Rdata')
write.table(exprSet,file ='GSE58294.txt',sep = '\t',col.names = NA,quote = F)






#######去批次

library(limma)
# 没有先安装
#BiocManager::install('sva',ask = F,update = F)
library(sva)

# 运行以下代码去批次
mergeFile="merge.preNorm.txt"            #合并后的文件名称
normalizeFile="merge.normalzie.txt"      #矫正后的文件名称
# 一定是在昨天输出这两个txt的文件夹下
files=c('GSE28750.txt','GSE57065.txt')       #输入文件名称
#获取交集基因
length(files)
library(data.table)

geneList=list()
for(i in 1:length(files)){
  fileName=files[i]
  rt=fread(fileName,data.table = F)
  header=unlist(strsplit(fileName, "\\.|\\-"))
  geneList[[header[1]]]=as.vector(rt[,1])
}
# 共同基因
intersectGenes=Reduce(intersect, geneList)
library(limma)
#数据合并
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
  fileName=files[i]
  header=unlist(strsplit(fileName, "\\.|\\-"))
  #读取输入文件，并对输入文件进行整理
  rt=fread(fileName,data.table = F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
  rt=avereps(data)
  colnames(rt)=paste0(header[1], "_", colnames(rt))
  #对数值大的数据取log2
  #qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  #LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
  #if(LogC){
  #rt[rt<0]=0
  #rt=log2(rt+1)}
  #rt=normalizeBetweenArrays(rt)
  #数据合并
  if(i==1){
    allTab=rt[intersectGenes,]
  }else{
    allTab=cbind(allTab, rt[intersectGenes,])
  }
  batchType=c(batchType, rep(header[1],ncol(rt)))
}

allTabOut=rbind(geneNames=colnames(allTab), allTab)
write.table(allTabOut, file=mergeFile, sep="\t", quote=F, col.names=F)

#对数据进行批次矫正，输出矫正后的结果
library(sva)
normalizeTab=ComBat(allTab, batchType, par.prior=TRUE)
normalizeTab=rbind(geneNames=colnames(normalizeTab), normalizeTab)
write.table(normalizeTab, file=normalizeFile, sep="\t", quote=F, col.names=F)

#合并前的PCA###################
#install.packages('ggplot2')
library(ggplot2)        #引用包
#读取输入文件,提取数据
rt=read.table('merge.preNorm.txt', header=T, sep="\t", check.names=F, row.names=1)
data=t(rt)
Project=gsub("(.*?)\\_.*", "\\1", rownames(data))
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))

#PCA分析
data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=Project)

#定义颜色
bioCol=c("#bc2c29",'#0072b5',"#20854e","#ef7c1c","#EE4C97","#FF9900","#20854E","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121")
bioCol=bioCol[1:length(levels(factor(Project)))]

#绘制PCA图

p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Type)) +
  scale_colour_manual(name="",  values=bioCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)

#合并后PCA###################
#读取输入文件,提取数据
rt=read.table('merge.normalzie.txt', header=T, sep="\t", check.names=F, row.names=1)
data=t(rt)
Project=gsub("(.*?)\\_.*", "\\1", rownames(data))
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))

#PCA分析
data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=Project)

#定义颜色
bioCol=c("#bc2c29",'#0072b5',"#20854e","#ef7c1c","#EE4C97","#FF9900","#20854E","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121")
bioCol=bioCol[1:length(levels(factor(Project)))]

#绘制PCA图

p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Type)) +
  scale_colour_manual(name="",  values=bioCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)




### 富集（运气做）
rt=read.table('merge.normalzie.txt', header=T, sep="\t", check.names=F, row.names=1)
# 牢记样本顺序

group_list1=c(rep('sepsis',10),rep('healthy',20))

group_list2=c(rep('sepsis',82),rep('healthy',25))

group_list=c(group_list1,group_list2)
group_list=factor(group_list,levels=c('healthy','sepsis'))
#3的数据集，提前2个
a=grep('GSE28750',colnames(rt))
b=grep('GSE57065',colnames(rt))

exprSet=rt[,c(a,b)]



gene=read.table('Senescence.txt',header = F,)

Senescence=list(gene$V1)
names(Senescence)='Senescence'


gene=read.table('NETs.txt',header = F,)

NETs=list(gene$V1)
names(NETs)='NETs'

library(GSVA)
exprSet=as.matrix(exprSet)
score=gsva(exprSet,NETs,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
score
max(score)
min(score)


score=as.data.frame(t(score))

score$group=group_list
library(ggpubr)

ggboxplot(score,x='group',y='NETs',
          color = "black",
          fill = "group",
          palette=c("#20854e","#E18727"))+stat_compare_means()




