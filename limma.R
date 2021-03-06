library(CLL)
data(sCLLex)
exprSet <- exprs(sCLLex)   
samples <- sampleNames(sCLLex)
pdata <- pData(sCLLex)
group_list <- as.character(pdata[,2])
dim(exprSet)

par(cex = 0.7)
n.sample <- ncol(exprSet)
if(n.sample > 40){
  par(cex = 0.5)
}

cols <- rainbow(n.sample*1.2)

boxplot(exprSet, col = cols, main = 'expression value', las = 2)


library(limma)
design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(exprSet)
design

contrast.matrix <- makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
contrast.matrix 
##step1
fit <- lmFit(exprSet,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput <- topTable(fit2, coef=1, n=Inf)
nrDEG <- na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(nrDEG)

