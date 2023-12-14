setwd('./scRNA/')
gc()
install.packages("Matrix")
#BiocManager::install('Seurat',update = F,ask = F)
library(Seurat)
#加载R包
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(clustree)
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(CellChat)
dir_name=c('N1','N2','S1','S2')
datalist=list()
for (i in 1:length(dir_name)){
  dir.10x = paste0("",dir_name[i])
  my.data <- Read10X(data.dir = dir.10x) 
  datalist[[i]]=CreateSeuratObject(counts = my.data, project = dir_name[i], 
                                   #每个基因至少在3个细胞中表达，每一个细胞至少有250个基因表达
                                   min.cells = 3, min.features = 250)
}

#修改名称
names(datalist)=dir_name
#添加分组信息
datalist[[1]]@meta.data$group <- 'N'
datalist[[2]]@meta.data$group <- 'N'
datalist[[3]]@meta.data$group <- 'S'
datalist[[4]]@meta.data$group <- 'S'
#2、细胞质控####
# 批量计算线粒体和rRNA占比
for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 计算线粒体占比
  datalist[[i]] <- sce
  rm(sce)
}
#循环绘制每个样本的QC图
qc_feature <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
for(i in 1:length(datalist)){
  pdf(file=paste0("1_",i,"_quality_control.pdf"),width = 15,height=7)
  print(VlnPlot(datalist[[i]], features = qc_feature, ncol = 5, pt.size = 0.5))
  dev.off()
}
#质控前的
violin=list()
for (i in 1:length(datalist)){
  violin[[i]] <- VlnPlot(datalist[[i]],
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                         pt.size = 0.1, 
                         ncol = 4)
}
pearplot_befor <- CombinePlots(plots = violin , nrow=length(datalist), legend="none")
pearplot_befor
ggsave(filename = 'QC_before.pdf',plot = pearplot_befor,he=15,wi=15)
#样本合并
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
#统计每一个样本的个数
raw_count <- table(sce@meta.data$orig.ident)
table(sce@meta.data$orig.ident)

#各个样本的细胞数目
#N1   N2   S1   S2 
#7593 8269 6750 4401 
pearplot_befor1<-VlnPlot(sce,
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                         pt.size = 0.1, 
                         ncol = 4)
pearplot_befor1
ggsave(filename = 'QC_before1.pdf',plot = pearplot_befor1,he=7,wi=15)

#循环绘制特征之间的相互关系图
for(i in 1:length(datalist)){
  pdf(file=paste0("1_",i,"_feature_relationship.pdf"),width = 12,height=10)
  p1 <- FeatureScatter(datalist[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p2 <- FeatureScatter(datalist[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  print(p1 + p2)
  dev.off()
}
#过滤细胞
datalist <- lapply(X = datalist, FUN = function(x) {
  x<-subset(x,subset = nFeature_RNA > 200 & 
              nFeature_RNA < 6000 & 
              percent.mt <= 10)
})#过滤掉具有 1.<200 个基因、>6000 个基因的细胞。2.过滤线粒体 >20% 线粒体基因的细胞。3.过滤核糖体基因大于1%的细胞

#再次整合
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
clean_count <- table(sce@meta.data$orig.ident)
table(sce@meta.data$orig.ident)

N1   N2   S1   S2 
6129 7465 6283 3561 
#过滤前后样本细胞数据的统计
summary_cells <- as.data.frame(cbind(raw_count,clean_count))
counts <- rbind(as.data.frame(cbind(summary_cells[,1],rep("raw",each = length(summary_cells[,1])))),
                as.data.frame(cbind(summary_cells[,2],rep("clean",each = length(summary_cells[,2])))))
counts$sample <- rep(rownames(summary_cells),times =2)
colnames(counts)<- c("count","Stat","sample")
counts[,1] <- as.numeric(counts[,1])
counts$Stat <- factor(counts$Stat, levels=c("raw", "clean"), ordered=TRUE)
fit_cell_count <- ggplot(data =counts, mapping = aes(x = sample, y=count))+ 
  geom_bar(aes(fill = Stat),stat = 'identity', position = 'dodge') + scale_fill_brewer(palette = "Set1") +
  theme(text=element_text(size=10),legend.title=element_blank(),
        panel.background = element_rect(fill = "white", colour = "black",size = 0.2),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.background = (element_rect(colour= "white",fill = "white")))

fit_cell_count
ggsave(filename = 'fit_cell_count.pdf',plot = fit_cell_count,width = 9,height = 9)
#质控后的小提琴图
violin_after=list()
for (i in 1:length(datalist)){
  violin_after[[i]] <- VlnPlot(datalist[[i]],
                               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                               pt.size = 0.1,
                               ncol = 4)
}
pearplot_after <- CombinePlots(plots = violin_after , nrow=length(datalist), legend="none")
pearplot_after
ggsave(filename = 'QC_after.pdf',plot = pearplot_after,he=15,wi=15)
pearplot_after1 <- VlnPlot(sce,
                           features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                           pt.size = 0.1,
                           ncol = 4)
pearplot_after1
ggsave(filename = 'QC_after1.pdf',plot = pearplot_after1,he=7,wi=15)
#质控前后图片的合并
pearplot_befor1
pearplot_after1
qc_merge<- CombinePlots(plots = list(pearplot_befor1,pearplot_after1) , 
                        nrow=2, legend='none')
qc_merge
ggsave(filename = 'qc_merge.pdf',plot = qc_merge,he=9,wi=15)

#3、数据预处理####
#合并数据,多样本整合
#归一化
datalist <- lapply(datalist,function(x) {
  NormalizeData(x)
})
#寻找高变异度基因
datalist <- lapply(datalist, function(x) {
  FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
datalist
#整合成一个对象,寻找锚点
sce.anchors <- FindIntegrationAnchors(object.list = datalist, dims = 1:20)
#根据锚点来整合
sce.combined <- IntegrateData(anchorset = sce.anchors, dims = 1:20)
DefaultAssay(sce.combined) <- "integrated" #更改默认数组

#保存datalist文件
save(sce.combined,file = 'sce.combined')
load("sce.combined")

#对整合后的数据进行尺度变换
#ScaleData
all.genes <- rownames(sce.combined[["RNA"]]@data) #针对所有基因
length(all.genes)
sce.combined <- ScaleData(sce.combined, features = all.genes)
###########对整合之后的数据进行降维
#PCA降维，选择合适的拐点
sce.combined <- RunPCA(sce.combined, features = VariableFeatures(sce.combined)) 
#可视化PCA降维之后的结果
pdf('2_PCA.pdf', width = 10, height = 6)
DimPlot(sce.combined, reduction = "pca")
dev.off()
#查看拐点
elbowplot1 <- ElbowPlot(sce.combined, ndims=50, reduction="pca") 
elbowplot1
ggsave(filename = 'sc_pca.pdf',plot = elbowplot1,he=10,wi=15)
#可视化前4个PC的top20个基因
VizDimLoadings(sce.combined, dims = 1:4, nfeatures = 20, reduction = "pca")
#前4个PC
DimHeatmap(sce.combined, dims = 1:4, cells = 500, balanced = TRUE)

#对整合之后的数据进行聚类

sce.combined <- RunUMAP(sce.combined, reduction = "pca", dims = 1:20)
sce.combined <- RunTSNE(sce.combined , reduction = "pca", dims = 1:20)
sce.combined <- FindNeighbors(sce.combined, reduction = "pca", dims = 1:20)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce.combined <- FindClusters(sce.combined,resolution = res)
} #设置不同的分辨率，观察分群效果(选择哪一个？)

#可视化不同分辨率，分群效果
apply(sce.combined @meta.data[,grep("integrated_snn",colnames(sce.combined @meta.data))],2,table)
colnames(sce.combined @meta.data)
pdf('2_umap_resolution_high.pdf', width =20)
plot_grid(ncol = 3, DimPlot(sce.combined , reduction = "umap", group.by = "integrated_snn_res.0.8") + 
            ggtitle("louvain_0.8"), DimPlot(sce.combined , reduction = "umap", group.by = "integrated_snn_res.1") + 
            ggtitle("louvain_1"), DimPlot(sce.combined , reduction = "umap", group.by = "integrated_snn_res.0.3") + 
            ggtitle("louvain_0.3"))
dev.off()
pdf('2_umap_resolution_low.pdf', width = 18)
plot_grid(ncol = 3, DimPlot(sce.combined , reduction = "umap", group.by = "integrated_snn_res.0.01") + 
            ggtitle("louvain_0.01"), DimPlot(sce.combined , reduction = "umap", group.by = "integrated_snn_res.0.1") + 
            ggtitle("louvain_0.1"), DimPlot(sce.combined , reduction = "umap", group.by = "integrated_snn_res.0.2") + 
            ggtitle("louvain_0.2"))
dev.off()
###聚类树
pdf('2_Tree_diff_resolution.pdf', width = 10,height = 10)
clustree(sce.combined @meta.data, prefix = "integrated_snn_res.",layout = "sugiyama")
dev.off()

#查看聚类图之后，那么接下来分析，按照分辨率为0.3进行
sel.clust = "integrated_snn_res.0.5"
sce.combined <- SetIdent(sce.combined, value = sel.clust)
table(sce.combined@active.ident)


#绘制UMAP聚类图
pdf('3_UMAP_multi_samples_combined.pdf',width = 11,height = 6)
p1 <- DimPlot( sce.combined , reduction = "umap", group.by = "group")
p2 <- DimPlot(sce.combined , reduction = "umap", group.by = "integrated_snn_res.0.5",label = TRUE)
p1 + p2
dev.off()

#绘制tsne聚类图
pdf('3_tsne_multi_samples_combined.pdf',width = 11,height = 6)
p1 <- DimPlot( sce.combined , reduction = "tsne", group.by = "group")
p2 <- DimPlot(sce.combined , reduction = "tsne", label = TRUE)
p1 + p2
dev.off()
#分别展示两组的UMAP聚类图：split.by
pdf('3_UMAP_multi_samples_split.pdf',width = 10,height = 6)
DimPlot(sce.combined, reduction = "umap", split.by = "group",label = TRUE)+ggsci::scale_color_igv()
dev.off()
#DefaultAssay(sce.combined) <- "integrated" #更改默认数组
DefaultAssay(sce.combined) <- "RNA"
pdf('3_FeaturePlot_cluster_markers.pdf',width = 10,height = 6)
FeaturePlot(sce.combined,features = sig_sense,order = T,label = T,pt.size = 1,split.by = 'group')
dev.off()

pdf('3_VlnPlot_cluster_markers.pdf',width = 10,height = 6)
VlnPlot(sce.combined,features = sig_sense)
dev.off()

VlnPlot(sce.combined,features = 'MMP9')


save(sce.combined,file ='sce.combined.Rdata')
load('sce.combined.Rdata')

DefaultAssay(sce.combined) <- "integrated"
library(SingleR)
refdata <- SingleR::HumanPrimaryCellAtlasData()
library(Seurat)
testdata <- GetAssayData(sce.combined, slot="data")
clusters <- sce.combined@meta.data$integrated_snn_res.0.5
cellpred <- SingleR(test = testdata, ref = refdata,
                    labels =refdata$label.main,
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
sce.combined@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.combined@meta.data[which(sce.combined@meta.data$integrated_snn_res.0.5 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.combined$celltype)
Idents(sce.combined)=sce.combined$celltype

save(sce.combined,file ='sce.combined.Rdata')
load('sce.combined.Rdata')

library(RColorBrewer) 
library(viridis)
library(wesanderson)
####COLOR PALETTE FOR PLOTS ####
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pie(rep(6,n), col=sample(color, n))
col_vector

col_vector =c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest1"), wes_palette("Cavalcanti1"), wes_palette("GrandBudapest2"), wes_palette("FantasticFox1"))
pal <- wes_palette("Zissou1", 10, type = "continuous")
pal2 <- wes_palette("Zissou1", 5, type = "continuous")
pal[3:10]

sig_sense=read.csv('diffSig_sence.csv')
sig_sense=sig_sense$X
load("sce.combined.Rdata")
## 17基因在各细胞内
DefaultAssay(sce.combined) <- "RNA"
pdf('vln_17gene_anno.pdf',width =10,height = 8)
VlnPlot(sce.combined, features = sig_sense, group.by = 'celltype',pt.size = 0.3,ncol = 5)
dev.off()

## 17基因在不同疾病程度内
pdf('vln_17gene_PATIENT.pdf',width =10,height = 8)
VlnPlot(sce.combined, features = sig_sense, group.by = 'orig.ident',pt.size = 0.3,ncol = 5)
dev.off()



DefaultAssay(sce.combined) <- "RNA"
sensescore=list(sig_sense)
names(sensescore)='sensescore'

library(Seurat)
sce.combined <- AddModuleScore(object = sce.combined, features = sensescore,
                               name ='sensescore',replace = T)

VlnPlot(sce.combined,features = 'sensescore1', 
        pt.size = 0, adjust = 2,group.by = "orig.ident")

library(ggplot2)
pdf('4_VlnPlot_celltype_markers.pdf',width = 20,height = 10)
VlnPlot(sce.combined,features ='sensescore1', pt.size = 0,split.by = 'celltype') 
dev.off()

pdf('4_FeaturePlot_celltype_markers.pdf',width = 20,height = 10)
FeaturePlot(sce.combined,features = 'sensescore1',order = F,split.by = 'celltype',label = T,cols =viridis(10))
dev.off()

table(sce.combined$celltype)
table(sce.combined$group)
## 显著性
NK<- subset(sce.combined, celltype %in% c("NK_cell"))
NK_N<- subset(NK, group %in%  c("N"))
NK_S<- subset(NK, group %in%  c("S"))
wilcox.test(NK_N$sensescore1, NK_S$sensescore1, alternative = "two.sided") #p-value < 2.2e-16

## 显著性
Monocyte<- subset(sce.combined, celltype %in% c("Monocyte"))
Monocyte_N<- subset(Monocyte, group %in%  c("N"))
Monocyte_S<- subset(Monocyte, group %in%  c("S"))
wilcox.test(Monocyte_N$sensescore1, Monocyte_S$sensescore1, alternative = "two.sided") #p-value < 2.2e-16

##################
##markergene
library(data.table)
library(tidyverse)
library(Seurat)

#gene=read.csv('genes.csv')
colnames(sce.combined@meta.data)
####
df <- as.data.frame(table(sce.combined$celltype, sce.combined$group))
colnames(df) <- c("celltype","sample","value")
p <- ggplot(df, aes(x=sample, y=value, fill=celltype)) + 
  geom_bar(stat= "identity", position = "fill", width = 0.9) + 
  ggsci::scale_fill_futurama() +
  scale_y_continuous(labels = scales::percent) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle = 90), axis.title = element_blank())
ggsave('celltype_proportion_1.pdf', p, width = 6, height = 4)
p <- ggplot(df, aes(x=sample, y=value, fill=celltype)) + 
  geom_bar(stat= "identity", position = "fill", width = 0.9) + 
  ggsci::scale_fill_futurama() +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle = 90), axis.title = element_blank())
ggsave('celltype_proportion_2.pdf', p, width = 6, height = 3.5)
#################美化单细胞绘图

pdf('5celltype.pdf',width = 20,height = 10)
DimPlot(sce.combined, reduction = "tsne",group.by = "celltype",
        label = T,label.box = T,label.size = 3,repel = T) + theme_bw() +
  labs( x= "tSNE 1",y= "tSNE 2",title = "cell type") +
  theme(panel.grid=element_blank(), # 去网格线
        plot.title = element_text(size = 15,color="black",hjust = 0.5),
        axis.text.x = element_text(size = 12, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        axis.title.x = element_text(size = 12, color = 'black'),
        axis.title.y = element_text(size = 12, color = 'black'),
        axis.ticks = element_line(color = 'black', lineend = 'round'),
        legend.position = 'bottom',
        legend.text = element_text(size = 12, color = 'black'),
        legend.title = element_text(size = 12, color = 'black'),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.65)) +
  scale_color_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.65))
dev.off()

#DimPlot(sce.combined, order=T, group.by = "celltype",label=T, cols = col_vector)

##4. 拟时序分析
library(Seurat)
#没有monocle要先安装 BiocManager::install
library(monocle)
library(tidyverse)
library(patchwork)

#挑出NK cells，这里你可以换

dir.create("pseudotime")
##提取细胞子集，NK cells可以换！！！
scRNAsub <- subset(sce.combined, celltype %in% c('NK_cell'))
dim(scRNAsub)
set.seed(123456)
#a=sample(1:ncol(scRNAsub),2000,replace = F)
#scRNAsub=scRNAsub[,a]
#dim(scRNAsub)
##重新降维聚类,提取子集后必须重新降维
scRNAsub <- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub <- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
ElbowPlot(scRNAsub, ndims=20, reduction="pca")
pc.num=1:10
##细胞聚类
scRNAsub <- FindNeighbors(scRNAsub, dims = pc.num) 
scRNAsub <- FindClusters(scRNAsub )

##都不用更改代码
#UMAP
scRNAsub <- RunUMAP(scRNAsub, dims = pc.num)
scRNAsub <- RunTSNE(scRNAsub, dims = pc.num)
DimPlot(scRNAsub, reduction = "tsne",label = T,cols = col_vector) 

VlnPlot(scRNAsub,features = 'sensescore1',cols = col_vector)

scRNAsub$group
VlnPlot(scRNAsub,features = 'sensescore1',group.by = 'orig.ident')

VlnPlot(scRNAsub,features =sig_sense,group.by = 'group')

FeaturePlot(scRNAsub,features = 'sensescore1',split.by = 'orig.ident',label = T,cols = c('#dadada','#bc3c29'))


setwd('./scRNA/')
scRNAsub=readRDS('./pseudotime/scRNA_NK.RDS')
dim(scRNAsub)
## 拟时序分析 先关掉再开
library(Seurat)
#没有monocle要先安装 BiocManager::install
library(monocle)
library(tidyverse)
library(patchwork)
data=as.matrix(scRNAsub@assays$RNA@counts)
data <- as(data, 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
## 以下代码一律不得修改
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())

mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)

##使用monocle选择的高变基因，不修改
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
plot_ordering_genes(mycds)


#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序    trace('project2MST', edit = T, where = asNamespace("monocle"))
mycds <- orderCells(mycds)

#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
plot1

##Pseudotime轨迹图
plot2 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
plot2

plot3 <-  plot_cell_trajectory(mycds, color_by = "orig.ident")
plot3

plot4 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
plot4

##合并出图
plotc <- plot1|plot2|plot3|plot4
pdf('6_monocle.pdf', width = 10, height = 6)
plotc
dev.off

sig_sense=read.csv('diffSig_sence.csv')
sig_sense=sig_sense$X

# 展现趋势统一的基因
### 指定基因的可视化
s.genes <- c("CD55","GMFG","IGFBP7","CCL5","CCL4","IL32")
p1 <- plot_genes_jitter(mycds[s.genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(mycds[s.genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "State")
plotc <- p1|p2|p3
ggsave("6.Genes_Jitterplot_NK.pdf", plot = plotc, width = 16, height = 8)



#######################
## 细胞通讯分析
# 先定义细胞
save(scRNAsub,file ='scRNAsub.Rdata')
#scRNAsub=subset(scRNA,celltype %in% 'NK')
summary(scRNAsub$sensescore1)
scRNAsub$sensegroup='sense_median'
save
for (i in 1:nrow(scRNAsub@meta.data)) {
  if (scRNAsub$sensescore1[i]>0.26529) {
    scRNAsub$sensegroup[i]='sense_high'
  }
  if (scRNAsub$sensescore1[i]<0.02596) {
    scRNAsub$sensegroup[i]='sense_low'
  }
}

table(sce.combined$celltype) 
sce.combined$sensegroup=''
# 排除内皮的其他细胞
scRNAother=subset(sce.combined, celltype != 'NK_cell')

scRNA_chat=merge(scRNAsub,c(scRNAother))

scRNA_chat$sensecell=paste0(scRNA_chat$sensegroup,scRNA_chat$celltype)

## 用最新的R
# BiocManager::install("sqjin/CellChat")
library(CellChat)

# 选取病鼠，如果内存不够，再进一步减少细胞数，例如随机抽2000个
#scRNA_chat <- subset(scRNA_chat, orig.ident=='MCAO1')

meta =scRNA_chat@meta.data # a dataframe with rownames containing cell mata data

data_input <- as.matrix(scRNA_chat@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "sensecell")

CellChatDB <- CellChatDB.human
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 

dplyr::glimpse(CellChatDB$interaction)
##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)

unique(cellchat@idents)

# 等待 计算通讯概率，推断细胞通讯网络
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)#获取所有的配受体对以及其通讯概率

write.csv(df.net,file ='cellchat.csv',quote=F)
#returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways

#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) 

cellchat <- aggregateNet(cellchat)#通过统计互作树或者讯通概率计算聚合网络
dev.off()

groupSize <- as.numeric(table(cellchat@idents))

#可视化加和的细胞间通讯网络。例如，使用circle plot显示任意两个细胞亚群之间的通讯次数或总通讯强度(权重)：

pdf('7_通讯1.pdf', width = 10, height = 6)par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
table(scRNA_chat$sensecell)
p_bubble= netVisual_bubble(cellchat,
                           sources.use = c('sense_highNK_cell','sense_lowNK_cell'),
                           targets.use = c('B_cell','Monocyte', 'Neutrophils','Platelets','T_cells'),
                           
                           remove.isolate = FALSE)+coord_flip()
pdf('7_通讯2.pdf', width = 10, height = 6)
p_bubble
dev.off()

#必须把上一个图关掉
cellchat@netP[["pathways"]]
pdf('7_通讯3.pdf', width = 10, height = 6)
netVisual_aggregate(cellchat, signaling = 'MIF',  
                    sources.use = c('sense_highNK_cell','sense_lowNK_cell'),
                    targets.use = c('B_cell','Monocyte', 'Neutrophils','Platelets','T_cells'))
dev.off()                    




#必须把上一个图关掉
dev.off()
netVisual_aggregate(cellchat, signaling = 'MIF',  
                    sources.use = c('Mitophagy_highEndothelial cells','Mitophagy_lowEndothelial cells'),
                    targets.use = c('Astrocytes','Epithelial cells', 'Microglia','Monocytes','Oligodendrocytes'))


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = 'MIF', width = 8, height = 2.5, font.size = 10)


h1=netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
h2=netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
h1 + h2


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg1

