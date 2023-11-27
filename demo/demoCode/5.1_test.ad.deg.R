library(Seurat)
library(dplyr)
library(ggplot2)
library(Libra)
library(ggrepel)

##load data
pathway<-"E:/paper1/fig6/" #Input file path
setwd(pathway)
#create a folder to store output files
dir.create('demo/output') 
data0<-readRDS(paste0(pathway,'integrated_ad_clean.rds') )

##randomly select cells to construct a subset for test

set.seed(666)
cells<-rownames(sample_frac(data0@meta.data,0.1))
data<-subset(data0,cells = cells)

## draw fig5a
DimPlot(data,group.by = 'group',label = T,reduction = 'umap')+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())+
  NoLegend()+
  ggtitle(NULL)
ggsave('test.ad.umap.pdf',width =6 ,height = 6)

## draw fig5b
DimPlot(data,group.by = 'phetp',reduction = 'umap',
        cols = c('#9467bd','#F2BF4A'),shuffle = T)+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())+
  NoLegend()+
  ggtitle(NULL)
ggsave('test.ad.umap.phetp.pdf',width =6 ,height = 6)


### DEG analysis 
# time ~ 8 min
wilcox<-run_de(data,meta=data@meta.data,cell_type_col = 'group',label_col = 'phetp',
               de_family = 'singlecell',de_method = 'wilcox',replicate_col = NULL)
write.csv(wilcox,file = 'demo/output/wilcox.csv')

#draw fig5d
# choose a cell_type
cell_type<-'GABA_Sox6'
de<-wilcox[wilcox$cell_type=='GABA_Sox6',]

vlcon<-de
vlcon<-subset(vlcon,avg_logFC!="NA" & abs(avg_logFC)!=Inf) #删除log2FoldChange中的空值
vlcon<-subset(vlcon,p_val!="NA" &p_val_adj!=0) #删除pvalue中的空值
#vlcon<-subset(vlcon,entrez!="NA")
vlcon<-vlcon[order(vlcon$p_val_adj,abs(vlcon$avg_logFC), decreasing = c(T,F)), ]
#对log2foldchange进行上调下调的判断
threshold <- factor(ifelse(vlcon$p_val_adj < 0.05 &
                             
                             abs(vlcon$avg_logFC) >= 0.1 ,
                           
                           ifelse(vlcon$avg_logFC >= 0.1 ,'Up','Down'),'Not'),
                    levels = c('Up','Down','Not'))
ggplot(vlcon,aes(x=avg_logFC,y=-log2(p_val_adj),colour=threshold)) +
  xlab("log2FC")+ylab("-log2(p.adj)") +
  geom_point(size = 2,alpha=0.8) +
  scale_color_manual(values=c("#BF242A","#1685A9","lightgrey"))+
  geom_vline(xintercept = c(-0.1, 0.1), lty = 2,colour="#000000")+ 
  geom_hline(yintercept = c(-log2(0.05)), lty = 2,colour="#000000")+
  geom_text_repel(
    data = vlcon[vlcon$p_val_adj<0.05&abs(vlcon$avg_logFC)>0.2,],
    aes(label = gene),
    size = 2,max.overlaps = 15,
    color = "black",
    segment.color = "black", show.legend = FALSE )+
  theme(
    axis.text=element_text(size=10),
    axis.title=element_text(size=10),
    axis.line = element_line(),
    panel.background = element_blank()
  )+
  NoLegend()
ggsave('test.de.pdf',width =6 ,height = 6)



