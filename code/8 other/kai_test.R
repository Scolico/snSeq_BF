setwd('C:\\Users\\AScol\\Desktop\\文章1\\data.part3')
library(Seurat)
VlnPlot(Chat.ad,'Ackr4',split.by = 'phetp')

VlnPlot(Chat.ad,features = c('Nrxn3','Nlgn3','Nlgn2','Nrxn2','Nlgn1','Nrxn1'),split.by = 'phetp')
VlnPlot(Chat.ad,'Flt3',split.by = 'phetp')

genes<-c('Nrxn3','Nlgn3','Nlgn2','Nrxn2','Nlgn1','Nrxn1','Flt3','Ackr4')
meta<-Chat.ad@meta.data

matrix<-Chat.ad@assays$RNA@counts
matrix<-NormalizeData(matrix)
#matrix<-scale(matrix)
#matrix<-NormalizeData(matrix,normalization.method='RC',scale.factor=1e6)
#matrix<-log2(matrix+1)

mat<-matrix[rownames(matrix) %in% genes,]
mat<-data.frame(t(data.frame(mat)))

mat$seurat_clusters<-meta$seurat_clusters
mat$phetp<-meta$phetp

library(ggunchained)
library(ggplot2)
term<-'Ackr4'
ggplot(mat,aes(x=seurat_clusters,
                      y=mat$Ackr4,
                      fill=phetp))+
  geom_violin(size = 0.2,trim = F,scale = 'area',adjust=1)+
  #geom_split_violin(size = 0.3,trim = F,adjust = 1.5)+
  #geom_boxplot(width=.5,
   #            position = position_dodge(.5),outlier.alpha = 1,size=0.3)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x =element_blank(), 
        axis.title = element_blank())+
  scale_fill_manual(values = c("#FFC843","#9242C1"))+
  #ylim(c(-1,1))+
  #facet_wrap('~seurat_clusters')+
  ggtitle(term)

ggsave(paste0(str_remove(term,':'),'.pdf'),width =3,height = 1.8)

###0值占比的卡方检验
kai.res<-data.frame(NULL)
mat2<-mat[,1:8]
mat2[mat2!=0]<-1
mat2<-cbind(mat2,mat[,9:10])

gene<-'Nrxn3'
test<-xtabs(~Nrxn3+phetp,data = mat2[mat2$seurat_clusters=='0',])

kai<-chisq.test(test)

res<-data.frame(gene=gene,`X-squared`=kai$statistic,p_value=kai$p.value)

kai.res<-rbind(kai.res,res)

###非0值细胞基因表达的表达差异
library(dplyr)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggtext)

gene<-'Flt3'
mat.n<-mat[,c(gene,'seurat_clusters','phetp')]
mat.n<-mat.n[mat.n[1] !=0,]
df_p_val1.Flt3 <- mat.n %>% 
  group_by(seurat_clusters) %>% 
  wilcox_test(formula = Flt3 ~ phetp) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='seurat_clusters')


pval<-df_p_val1.Flt3
ggplot(data = mat.n,mapping = aes(x=seurat_clusters,y=gene))+
  #geom_violin(data = mat.n,aes(x=seurat_clusters,y=.data[[gene]],fill=phetp))+
  stat_boxplot(data = mat.n,aes(x=seurat_clusters,y=.data[[gene]],fill=phetp),
               geom = "errorbar", width = 0.2,size=.4,
               position = position_dodge(.8)) + 
  geom_boxplot(data = mat.n,aes(x=seurat_clusters,y=.data[[gene]],fill=phetp),
               width=0.5,size=.4,
               outlier.size = .3,outlier.colour = 'black',
               position = position_dodge(.8),coef=1.5)+
  stat_summary(data = mat.n,aes(x=seurat_clusters,y=.data[[gene]],fill=phetp),
               fun=mean,geom="point",colour="magenta",alpha=1,size=1,shape=18,
               position = position_dodge(.8))+
  scale_fill_manual(values = c("#9242C1","#FFC843"))+
  stat_pvalue_manual(pval,label = '{p.signif}',tip.length = 0)+
  labs(x='seurat_clusters',y=gene)+
  guides(fill=guide_legend(title = 'group'))+
  theme_test()+
  theme(axis.text = element_text(color = 'black'),
        axis.title.x = element_blank())

ggsave(paste0(gene,'.pdf'),width =4,height = 2.38)