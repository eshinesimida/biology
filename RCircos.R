library(RCircos)

data(UCSC.HG38.Human.CytoBandIdeogram)
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info, chr.exclude=NULL,tracks.inside=10, tracks.outside=0 )

RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot() 

#导入内建人类染色体数据
a <- function(...){ 
  data(UCSC.HG19.Human.CytoBandIdeogram);
  head(UCSC.HG19.Human.CytoBandIdeogram);
  ## 这里换了个参考基因组版本，请注意
  chr.exclude <- NULL;  
  #设置不显示的染色体，如 c(1,3)          
  cyto.info <- UCSC.HG19.Human.CytoBandIdeogram; 
  #设置染色体数据
  tracks.inside <- 10;  
  #设置内部track 个数
  tracks.outside <- 0;  
  #设置外部track 个数
  
  #导入上面四个基本参数
  RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside);
  
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()
}
a()

b <- function(...){
  a()
  #Gene Labels and connectors on RCircos Plot
  #RCircos.Gene.Connector.Plot   绘制染色体表意文字和基因名称之间的连接
  #RCircos.Gene.Name.Plot 在数据轨道上绘制基因名称
  data(RCircos.Gene.Label.Data);
  name.col <- 4;
  side <- "in";
  track.num <- 1;
  RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side);
  track.num <- 2;
  RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col,track.num, side);
}
b()


data(RCircos.Heatmap.Data);
head(RCircos.Heatmap.Data);

c <- function(...){
  b()
  data.col <- 6;
  track.num <- 5;
  side <- "in";
  RCircos.Heatmap.Plot(RCircos.Heatmap.Data, data.col, track.num, side);
  
}
c()


data(RCircos.Scatter.Data);
head(RCircos.Scatter.Data);


d <- function(...){
  c()
  RCircos.Scatter.Data$chromosome=paste0('chr',RCircos.Scatter.Data$chromosome)
  head(RCircos.Scatter.Data);
  data.col <- 5;
  track.num <- 6;
  side <- "in";
  by.fold <- 1;
  RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col,track.num, side, by.fold);
  
  
}
d()

data(RCircos.Line.Data);
head(RCircos.Line.Data);

e <- function(...){
  d()
  RCircos.Line.Data$chromosome=paste0('chr',RCircos.Line.Data$chromosome)
  head(RCircos.Line.Data);
  data.col <- 5;
  track.num <- 7;
  side <- "in";
  RCircos.Line.Plot(RCircos.Line.Data, data.col, track.num, side);
  
}
e()


data(RCircos.Histogram.Data);
head(RCircos.Histogram.Data);
f <- function(...){
  e()
  
  data.col <- 4;
  track.num <- 8;
  side <- "in";
  RCircos.Histogram.Plot(RCircos.Histogram.Data, data.col, track.num, side);
  
}
f()

data(RCircos.Tile.Data);
head(RCircos.Tile.Data)



g <- function(...){
  f()
  
  
  track.num <- 9;
  side <- "in";
  RCircos.Tile.Plot(RCircos.Tile.Data, track.num, side);
}
g()

data(RCircos.Link.Data);
head(RCircos.Link.Data);
h <- function(...){
  g()
  
  
  track.num <- 11;
  RCircos.Link.Plot(RCircos.Link.Data, track.num, TRUE);
  
  data(RCircos.Ribbon.Data);
  RCircos.Ribbon.Plot(ribbon.data=RCircos.Ribbon.Data, track.num=11, by.chromosome=FALSE, twist=FALSE)
  
}
h()
