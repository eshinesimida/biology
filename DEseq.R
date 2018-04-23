options(warn=-1)
library(DESeq2)
library(limma)
library(pasilla)
data(pasillaGenes)
exprSet <- counts(pasillaGenes)
head(exprSet)  ##表达矩阵如下

group_list <- pasillaGenes$condition

colData <- data.frame(row.names=colnames(exprSet), 
                      group_list=group_list
)
## 这是一个复杂的方法构造这个对象！
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ group_list)



library(airway)
data(airway)
library(DESeq2)
dds1  <- DESeqDataSet(airway, design = ~ cell+ dex)
design(dds1) <- ~ dex  



dds2 <- DESeq(dds1)

##直接用DESeq函数即可
## 下面的代码如果你不感兴趣不需要运行，免得误导你
## 就是看看normalization前面的数据分布差异
rld <- rlogTransformation(dds2)  ## 得到经过DESeq2软件normlization的表达矩阵！
exprSet_new <- assay(rld)
par(cex = 0.7)
n.sample <- ncol(exprSet)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))
boxplot(exprSet, col = cols,main="expression value",las=2)
boxplot(exprSet_new, col = cols,main="expression value",las=2)
hist(exprSet)
hist(exprSet_new)

