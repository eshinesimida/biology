options(warn=-1)
library(DESeq2)
library(limma)
library(pasilla)
data(pasillaGenes)
exprSet <- counts(pasillaGenes)
head(exprSet)  
group_list <- pasillaGenes$condition

colData <- data.frame(row.names=colnames(exprSet), 
                      group_list=group_list
)

dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ group_list)

#normalization

dds2 <- DESeq(dds)
rld <- rlogTransformation(dds2)  
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

#
resultsNames(dds2)

res <-  results(dds2, contrast=c("group_list","treated","untreated")) 


resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
head(resOrdered)
