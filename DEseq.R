options(warn=-1)
library(DESeq2)
library(limma)
library(pasilla)
data(pasillaGenes)
exprSet <- counts(pasillaGenes)
head(exprSet)  ##�����������

group_list <- pasillaGenes$condition

colData <- data.frame(row.names=colnames(exprSet), 
                      group_list=group_list
)
## ����һ�����ӵķ��������������
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ group_list)



library(airway)
data(airway)
library(DESeq2)
dds1  <- DESeqDataSet(airway, design = ~ cell+ dex)
design(dds1) <- ~ dex  



dds2 <- DESeq(dds1)

##ֱ����DESeq��������
## ����Ĵ�������㲻����Ȥ����Ҫ���У��������
## ���ǿ���normalizationǰ������ݷֲ�����
rld <- rlogTransformation(dds2)  ## �õ�����DESeq2����normlization�ı������
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
