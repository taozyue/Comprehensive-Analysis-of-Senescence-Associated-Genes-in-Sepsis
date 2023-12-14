
rt=read.table('merge.normalzie.txt', header=T, sep="\t", check.names=F, row.names=1)
# 牢记样本顺序

group_list1=c(rep('sepsis',10),rep('healthy',20))

group_list2=c(rep('sepsis',82),rep('healthy',25))


group_list=c(group_list1,group_list2)
group_list=factor(group_list,levels=c('healthy','sepsis'))

a=grep('GSE28750',colnames(rt))
b=grep('GSE57065',colnames(rt))

exprSet=rt[,c(a,b)]

library(limma)
design=model.matrix(~ group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
#tumor处需修改case
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1) 
write.table(allDiff,file=paste0("2.","DIFF_all.xls"),sep="\t",quote=F,col.names=F)
diffSig=allDiff[with(allDiff, (abs(logFC)>0.5 &P.Value < 0.05)), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file=paste0("2.","DIFF_less.xls"),sep="\t",quote=F,col.names=F)



gene=read.table('Senescence.txt',header = F,)
gene=gene$V1


alldiff_sence=allDiff[rownames(allDiff) %in% gene,]
write.csv(alldiff_sence,file ='alldiff_sence.csv',quote = F)
diffSig_sence=diffSig[rownames(diffSig) %in% gene,]
#alldiff_sig=alldiff_sence[alldiff_sence$adj.P.Val<0.05,]
write.csv(diffSig_sence,file ='diffSig_sence.csv',quote = F)

### 作热图
rt1=exprSet[,group_list=='healthy']
rt2=exprSet[,group_list=='sepsis']

rt=cbind(rt1,rt2)

anno=data.frame(row.names = colnames(rt),group=c(rep('healthy',45),rep('sepsis',92)))
save(rt,anno,file ='key_train_exprSet.Rdata')


rt=rt[rownames(rt) %in% rownames(diffSig_sence),]

pheatmap::pheatmap(rt,scale = 'row',cluster_cols = F,border_color = NA, fontsize = 9,
                   fontsize_row=10,
                   fontsize_col=14,
                   breaks = seq(-1.6,1.6,length.out = 90),show_colnames = F,annotation_col = anno)

## 保存整理好的这些文件
### 火山图

library(ggplot2)

diffSig_sence$label=ifelse(diffSig_sence$adj.P.Val<0.05,rownames(diffSig_sence),'')


ggplot(diffSig_sence,aes(logFC, -log10(adj.P.Val)))+
  # 横向水平参考线：
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  # 纵向垂直参考线：
  geom_vline(xintercept = c(-0.2,0.2,0.5,-0.5), linetype = "dashed", color = "#999999")+
  # 散点图:
  geom_point(aes(size=-log10(adj.P.Val), color= -log10(adj.P.Val)))+
  # 指定颜色渐变模式：
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  # 指定散点大小渐变模式：
  #scale_size_continuous(range = c(1,3))+
  # 主题调整：
  theme_bw()+# 调整主题和图例位置：
theme(panel.grid = element_blank()
)+
  # 设置部分图例不显示：
  guides(col = guide_colourbar(title = "-Log10(adj.P.Val)"),
         size = "none")+
  # 添加标签：
  geom_text(aes(label=label, color = -log10(adj.P.Val)), size = 3, vjust = 1.5, hjust=1)+
  # 修改坐标轴：
  xlab("Log2FC")+
  ylab("-Log10(adj.P.Val)")

### 箱线图
rt=rt[rownames(diffSig_sence),]
rt=as.data.frame(t(rt))
rt$group=anno$group
ss=intersect(gene,colnames(rt))
rt=rt[,c(ss,'group')]

library(tidyverse)
rt2=tidyr::pivot_longer(rt,cols = -c('group'),names_to = "gene",values_to = 'expression',
                        )


library(ggpubr)

ggboxplot(
  rt2,
  x = "gene",
  y = "expression",
  color = "black",
  fill = "group",
  xlab = "group",
  ylab = "expression", palette=c("#20854e","#E18727")
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

### 相关性
load('key_train_exprSet.Rdata')

#gene=read.table('Senescence.txt',header = F,)
#gene=gene$V1

rt=rt[rownames(rt) %in% rownames(diffSig_sence),]
rt=as.data.frame(t(rt))

library(corrplot)
M <- cor(rt)
res1 <- cor.mtest(rt, conf.level = .95) 
library(ggcorrplot)
ggcorrplot(
  M,
  hc.order = F,
  #type = "lower",
  outline.color = NA,
  ggtheme = ggplot2::theme_gray,
  colors = c("blue", "white", "red")
)

rt=rt[anno$group=='sepsis',]
library(corrplot)
M <- cor(rt)
res1 <- cor.mtest(rt, conf.level = .95) 
library(ggcorrplot)
ggcorrplot(
  M,
  hc.order = F,
  #type = "lower",
  outline.color = NA,
  ggtheme = ggplot2::theme_gray,
  colors = c("blue", "white", "red")
)

sig_sense=read.csv('diffSig_sence.csv')
sig_sense=sig_sense$X

rt=rt[,sig_sense]
library(corrplot)
M <- cor(rt)
res1 <- cor.mtest(rt, conf.level = .95) 
library(ggcorrplot)
ggcorrplot(
  M,
  hc.order = F,
  #type = "lower",
  outline.color = NA,
  ggtheme = ggplot2::theme_gray,
  colors = c("blue", "white", "red")
)

library(ggstatsplot)
install.packages("ggstatsplot")
ggstatsplot::ggscatterstats(rt,x='MMP9','CD55')
ggstatsplot::ggscatterstats(rt,x='RPS27A','TOMM20')


###


# 加载数据集

load('key_train_exprSet.Rdata')
set.seed(12345)


library(caret)
sig_sense=read.csv('diffSig_sence.csv')
sig_sense=sig_sense$X

rt=rt[rownames(rt) %in%sig_sense,]
rt=as.data.frame(t(rt))

rt$group=anno$group
rt$group=factor(rt$group)
control <- rfeControl(functions = caretFuncs, method = "cv", number = 5)


set.seed(245083)
# 执行SVM-RFE算法
results <- rfe(rt[,1:17], 
               rt[,18],
               rfeControl = control,
               method = "svmRadial")

# 结果分析
print(results)
# 列出选择的变量集
rfe_gene=predictors(results)[1:10]


######################################
#############随机森林（减少变量）#################
####################################
#install.packages("randomForest")
library(randomForest)
set.seed(245083)
#!!!
x=rt[,1:17]
rf=randomForest(x = x,y = rt$group,ntree = 1000)

#rf=randomForest(group~., data=data, ntree=500)

plot(rf, main="Random forest", lwd=2)
#找出误差最小的点
optionTrees=which.min(rf$err.rate[,1])
optionTrees
set.seed(245083)
rf2=randomForest(x = x,y = rt$group, ntree=optionTrees)
#查看基因的重要性
importance=importance(x=rf2)
dev.off()
#绘制基因的重要性图
varImpPlot(rf2, main="")

#挑选疾病特征基因
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
# 可以调节！
rfGenes=names(rfGenes[1:10])     #挑选重要性评分大于0的基因

gene=intersect(rfe_gene,rfGenes)
gene

save(gene,file ='gene_rfe_rf.Rdata')




### 多因素逻辑回归
load('gene_rfe_rf.Rdata')
rt=rt[,c(gene,'group')]
rt$group=ifelse(rt$group=='healthy',0,1)


glm=glm(group ~.,data = rt,family= binomial(link='logit'))

OR=exp(glm$coefficients)[2:(length(gene)+1)]
OR_CI=exp(confint(glm,level = 0.95))
OR.95L=OR_CI[2:(length(gene)+1),1]
OR.95H=OR_CI[2:(length(gene)+1),2]


df=data.frame(gene=gene,OR=OR,OR.95L=OR.95L,OR.95H=OR.95H)


#install.packages("forestploter")
library(forestploter)

df$`OR (95% CI)` <- ifelse(is.na(df$OR), "",
                           sprintf("%.2f (%.2f to %.2f)",
                                   df$OR,df$OR.95L,df$OR.95H))


df$"" = paste(rep(" ",57), collapse = " ")
plot <- forest(df[, c(1, 6, 5)],
               est = df$OR,
               lower =df$OR.95L,
               upper =df$OR.95H,ref_line = 1,
               ci_column = 2)
plot

pred=predict(glm,rt,type = 'response')
pred_df=as.data.frame(pred)
save(pred_df,file ='pred_logistic.Rdata')



library(pROC)
plot.roc(rt$group,pred,print.auc=T,print.thres=T)

###bootstrap验证
rsq <- function(formula, data, indices) { 
  # 通过data得到重抽样样本
  d <- data[indices,] 
  # 通过formula得到模型
  fit.glm <- glm(formula, data=d)
  # 用boot的模型在原数据集上预测
  glm.pred <- predict(fit.glm, data) 
  # 得到预测的ROC
  glmROC=roc(rt$group,as.numeric(glm.pred))
  plot.roc(rt$group,glm.pred,add=T,col='#0072b5')
  # 返回ROC的AUC
  return(glmROC$auc)
} 

library(boot) 
set.seed(123456) 


# 运行bootstrap
results <- boot(data=rt, statistic=rsq, 
                R=1000, formula=group~.)

#
print(results)
plot(results)

mean(results$t)

# 获取置信区间
boot.ci(results, type=c("perc", "bca"))

rm(outtab)

###boot 灵敏度
rsq <- function(formula, data, indices) { 
  # 通过data得到重抽样样本
  d <- data[indices,] 
  # 通过formula得到模型
  fit.glm <- glm(formula, data=d)
  # 用boot的模型在原数据集上预测
  glm.pred <- predict(fit.glm, data,type='response') 
  print(glm.pred)
  b=ifelse(glm.pred>0.523,1,0)
  outtab=table(b,rt$group)
  sen=outtab[2,2]/(outtab[1,2]+outtab[2,2])
  spe=outtab[1,1]/(outtab[1,1]+outtab[2,1])
  return(sen)
} 

library(boot) 
set.seed(123456) 

# 运行bootstrap
results <- boot(data=rt, statistic=rsq, 
                R=1000, formula=group~.)

#
print(results)
plot(results)

## 和我们一样0.96
mean(results$t)

# 获取置信区间
boot.ci(results, type=c("perc", "bca"))



###boot 特异度
rsq <- function(formula, data, indices) { 
  # 通过data得到重抽样样本
  d <- data[indices,] 
  # 通过formula得到模型
  fit.glm <- glm(formula, data=d)
  # 用boot的模型在原数据集上预测
  glm.pred <- predict(fit.glm, data,type='response') 
  print(glm.pred)
  b=ifelse(glm.pred>0.523,1,0)
  outtab=table(b,rt$group)
  sen=outtab[2,2]/(outtab[1,2]+outtab[2,2])
  spe=outtab[1,1]/(outtab[1,1]+outtab[2,1])
  return(spe)
} 

library(boot) 
set.seed(123456) 

# 运行bootstrap
results <- boot(data=rt, statistic=rsq, 
                R=1000, formula=group~.)

#
print(results)
plot(results)

## 和我们一样0.96
mean(results$t)

# 获取置信区间
boot.ci(results, type=c("perc", "bca"))




#install.packages('regplot')
library(regplot)
regplot(glm,observation = rt[2,],interval = 'confidence',clickable = T,title = 'Nomogram')


library(rms)
glm2=lrm(group~., data = rt,x=T,y=T)

cal<-calibrate(glm2, method = 'boot', B=1000, data = rt)
par(oma=c(3,3,3,3)) 
par(mar=c(6,5,4,3) + 0.1) 
plot(cal,
     xlim=c(0,1.0),ylim=c(0,1.0),
     xlab = "Predicted Probability",
     ylab = "Observed Probability",main = "Calibration Curve")

#install.packages('set')
#devtools::install_github("yikeshu0611/do")
#if (!require(ggDCA)) {
#  devtools::install_github("yikeshu0611/ggDCA")
#}
#install.packages("rmda")
library(rmda)
library(ggDCA)
library(ggplot2)
library(rms)
library(caret)


set.seed(123)
colnames(rt)

model1 <- decision_curve(group ~ MMP9, data = rt,
                         thresholds = seq(0, 0.4, by = 0.005), bootstraps = 10)
model2 <- decision_curve(group ~ CD55, data = rt,
                         thresholds = seq(0, 0.4, by = 0.005), bootstraps = 10)
model3 <- decision_curve(group ~ HGF, data = rt,
                         thresholds = seq(0, 0.4, by = 0.005), bootstraps = 10)
model4 <- decision_curve(group ~ ETS2, data = rt,
                         thresholds = seq(0, 0.4, by = 0.005), bootstraps = 10)

colnames(rt)
#!!!!!!!!!!!!!!!!!!
model0 <- decision_curve(group ~ MMP9 +CD55+ HGF +ETS2+ IL18+IL10+IGFBP7,
                         data = rt,thresholds = seq(0, 0.4, by = 0.005), bootstraps = 10)

cols<-c("#88c4e8","#db6968","#982b2b","#0074b3",
        "#e5ce81","#f47720","#459943","#bdc3d2")
par(oma=c(1,1,1,1)) 
par(mar=c(6,5,4,3) + 0.1) 
plot_decision_curve( list(model1,model2,model3,model4,model0),
                     # curve.names = c("Baseline model", "Full model"),
                      col = cols,
                     confidence.intervals = F,  #remove confidence intervals
                     cost.benefit.axis = FALSE, #remove cost benefit axis
                     legend.position = "bottomleft") #add the legend
