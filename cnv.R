library(scRNAseq) 
data(fluidigm)
exprSet=assays(fluidigm)$rsem_counts
library(edgeR)
exprSet=log(cpm(exprSet)+1)
boxplot(exprSet[,1:5],las=2)
exprSet[1:4,1:4]
pos = read.table('../airway/human.gene.positions')
head(pos)
exprSet=exprSet[rownames(exprSet) %in% pos[,7],]
pos=pos[match(rownames(exprSet),pos[,7]),c(7,2:4)]
table(pos[,2])
dim(pos)
dim(exprSet)

new_chr=gsub('chr','',pos[,2])
table(new_chr)
new_chr[new_chr=='X']=23
new_chr[new_chr=='Y']=24
new_chr=as.numeric(new_chr)
pos[,2]=new_chr
pos=pos[order(pos[,2],pos[,3]),] 
write.table(pos,'scRNAseq_pos.txt',row.names = F,col.names = F,sep = '\t',quote = F)
write.table(exprSet,'scRNAseq_exprSet.txt',quote = F,sep = '\t')



rm(list = ls())
load(file = 'GBM_for_CNV_input.Rdata')
colnames(pos)=c('gene','chr','start','end')
exprSet[1:4,1:4]
head(pos)
exprSet <- exprSet[match(pos$gene,rownames(exprSet)),]
res <- cbind(pos,exprSet)
head(res[,1:10])
table(res$chr)

all_cnv <- lapply(split(res,res$chr), function(x){
  # x=split(res,res$chr)[[1]]
  anno=x[,1:4]
  dat=x[,5:434]
  if(nrow(dat)>100){
    cnv <- lapply(51:(nrow(dat)-50), function(i){
      this_cnv <- unlist( lapply(1:ncol(dat), function(j){
        sum(dat[(i-50):(i+50),j])/101
      }))
      return(this_cnv)
    })
    cnv=do.call(rbind,cnv)
    cnv=cbind(anno[51:(nrow(x)-50),],cnv)
    # cnv[1:4,1:8]
  }else{
    return(NULL)
  }
})
all_cnv=do.call(rbind,all_cnv) 
head(all_cnv[1:4,1:8])
table(all_cnv$chr)

library(pheatmap)
## 按照 列（样本）进行归一化
D=t(scale(all_cnv[,5:ncol(all_cnv)] ))
## 归一化之后进行转置，因为我们的热图列是基因，行是样本。

## 根据阈值进行截断
D[D>2]=2
D[D< -2] = -2

dim(D)
colnames(D)=paste0('genes_',1:ncol(D))
rownames(D)=colnames(exprSet)

require(RColorBrewer)
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

library(stringr)
annotation_row = data.frame(
  patients=str_split(rownames(D),'_',simplify = T)[,1]
)
rownames(annotation_row) = rownames(D)

annotation_col = data.frame(
  chr= factor(all_cnv$chr,levels = unique(all_cnv$chr))
)
rownames(annotation_col) = colnames(D)
png('GBM_CNV_by_jimmy.png')
pdf('GBM_CNV_by_jimmy.pdf')
pheatmap(D,cluster_rows = T,col=rev(cols),
         annotation_col=annotation_col,
         annotation_row = annotation_row,
         cluster_cols = F,show_rownames=F,show_colnames=F)

dev.off()
