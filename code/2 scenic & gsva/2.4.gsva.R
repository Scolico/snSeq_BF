setwd('G:\\paper1\\fig1.2\\gsva')

sp<-read.csv('sp.csv',row.names = 1)
age.genes.freq<-read.csv('age.genes.freq.0.6.csv',row.names = 1)

age.genes.freq<-data.frame(table(
  sp$gene[abs(sp$spearman.cor)>0.6]
  ))
age.genes.unique<-age.genes.freq$Var1
#[age.genes.freq$Freq==1]


sp.unique<-sp[sp$gene %in% age.genes.unique,]

library(reshape2)
top3<-sp.unique %>% group_by(cluster) %>% top_n(n=5,wt=abs(spearman.cor))
nrow=length(unique(top3$gene))
pheat.dat<-matrix(nrow=nrow,ncol=36,0)
colnames(pheat.dat)<-unique(sp.unique$cluster)
rownames(pheat.dat)<-unique(top3$gene)
library(dplyr)

dat<-sp[sp$gene %in% unique(top3$gene),]
for( i in 1:length(rownames(dat))){
  sp.cor=dat[i,1]
  cluster=dat[i,2]
  gene=dat[i,3]
  pheat.dat[gene,cluster]<-sp.cor
}

library(pheatmap)
library(viridis)
pheatmap(pheat.dat,
         cluster_cols = F,
         cluster_rows = T,
         scale = 'none',
         clustering_method = 'median',
         display_numbers = matrix(ifelse(abs(pheat.dat) > .6, "*", ""), nrow(pheat.dat)),
         number_color = 'white',
         color = viridis(100))

         