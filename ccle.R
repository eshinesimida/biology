library(GEOquery)
ccleFromGEO <- getGEO("GSE36133")
annotBlock1 <- pData(phenoData(ccleFromGEO[[1]]))
exprSet <- exprs(ccleFromGEO[[1]])

keyColumns <- c('title', 'source_name_ch1','characteristics_ch1','characteristics_ch1.1',
                'characteristics_ch1.2')

options(stringsAsFactors = F)
allAnnot <- annotBlock1[,keyColumns]

