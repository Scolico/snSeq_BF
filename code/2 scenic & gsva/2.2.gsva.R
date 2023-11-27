library(Seurat)
library(ggplot2)

data<-test.clean
meta.data<-data@meta.data
matrix<-data@assays$RNA@counts
matrix<-NormalizeData(matrix,normalization.method = 'RC')
##与年龄相关的基因
data$time.sp<-0
data$time.sp[data$orig.ident %in% c('ctp4')]<-'P4'
data$time.sp[data$orig.ident %in% c('ctp14')]<-'P14'
data$time.sp[data$orig.ident %in% c('ctm3')]<-'3M'
data$time.sp[data$orig.ident %in% c('ctm9')]<-'9M'
data$time.sp[data$orig.ident %in% c('ctm15')]<-'15M'


meta.data$time.sp[meta.data$orig.ident %in% c('ctp4')]<-1
meta.data$time.sp[meta.data$orig.ident %in% c('ctp14')]<-2
meta.data$time.sp[meta.data$orig.ident %in% c('ctm3')]<-3
meta.data$time.sp[meta.data$orig.ident %in% c('ctm9')]<-4
meta.data$time.sp[meta.data$orig.ident %in% c('ctm15')]<-5
mmeta<-meta.data

sp<-do.call(rbind,tapply(rownames(mmeta),mmeta$subtype,function(i){
  cell<-i
  #cell<-rownames(mmeta)[mmeta$subclass=='MGL']
  cluster<-unique(mmeta$subtype[rownames(mmeta) %in% cell])
  time.sp<-mmeta$time.sp[rownames(mmeta) %in% cell]
  mt<-matrix[,colnames(matrix) %in% cell]
  mt<-t(as.matrix(mt))
  print(paste0('calulate:',cluster))
  pb<-cor(mt,as.numeric(time.sp),method = 'spearman')
  colnames(pb)<-'spearman.cor'
  pb<-data.frame(pb,cluster=cluster)
  pb<-na.omit(pb)
}))

ggplot(sp,aes(x=cluster,y=spearman.cor))+
  stat_boxplot(geom = 'errorbar',width=0.3,position = position_dodge(1),size=0.2)+
  geom_boxplot(aes(fill=cluster),position = position_dodge(1),outlier.size = 0.1,size=0.2)+
  #geom_hline(yintercept = c(0), lty = 1,size=0.3,colour="blue")+
  geom_hline(yintercept = c(0.6,-0.6), lty = 2,size=0.3,colour="red")+
  theme(axis.title.x = element_blank(),
        axis.line.y = element_line(linetype = 1),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1))+
  ylim(-.9,.9)+
  #coord_flip()+
  ggtitle('Distribution of Sperman Correlation Coefficient')+
  scale_color_manual(values = rainbow(36))+
  NoLegend()




library(stringr)
a<-str_split(rownames(sp),'.')
b<-str_remove(rownames(sp),pattern = '[A-z]+\\.')
c<-str_remove(b,pattern = '[0-9]+\\.')
sp$gene<-c

#age.genes<-unique(sp$gene[abs(sp$spearman.cor)>0.8])
age.genes.freq<-data.frame(table(sp$gene[abs(sp$spearman.cor)>0.6]))
age.genes.unique<-age.genes.freq$Var1
#[age.genes.freq$Freq==1]

sp.unique<-sp[sp$gene %in% age.genes.unique,]


library(msigdbr)
library(GSVA)
library(Seurat)
library(stringr)
library(Matrix)
library(limma)
###选择参考基因集###
msig<-msigdbr_collections()
msgdC2 = msigdbr(species = "Mus musculus", category = "C5",subcategory = "BP")
geneSet1 = msgdC2 %>% split(x = .$gene_symbol, f = .$gs_name)

msgdC2 = msigdbr(species = "Mus musculus", category = "C5",subcategory = "CC")
geneSet2 = msgdC2 %>% split(x = .$gene_symbol, f = .$gs_name)

msgdC2 = msigdbr(species = "Mus musculus",category = "C5",subcategory = "MF")
geneSet3 = msgdC2 %>% split(x = .$gene_symbol, f = .$gs_name)

msgdC2 = msigdbr(species = "Mus musculus",category = "C2", subcategory = "CP:KEGG")
geneSet4 = msgdC2 %>% split(x = .$gene_symbol, f = .$gs_name)

msgdC2 = msigdbr(species = "Mus musculus",category = "C2", subcategory = "CP:REACTOME")
geneSet5 = msgdC2 %>% split(x = .$gene_symbol, f = .$gs_name)

my_genesets<-list()
for(i in list(geneSet1,geneSet2,geneSet3,geneSet4,geneSet5)){
  my_genesets<-append(my_genesets,i)
}

###整体分析得到的差异基因

my_genelist<-as.character(age.genes.unique)

###差异基因所涉及的通路
library(tidyverse)
a<-unlist(lapply(my_genesets,function(a){
  isFALSE(is_empty(intersect(a,my_genelist)))
}))
targetSets<-my_genesets[a]

###计算gsva分数

matrix<-data@assays$RNA@counts
matrix<-NormalizeData(matrix,normalization.method='RC',scale.factor=1e6)
matrix<-Matrix(log2(matrix+1),sparse=T)

aa<-tapply(rownames(mmeta),mmeta$subtype,list)  
aa<-tapply(rownames(mmeta),mmeta$subtype,function(a){
  m<-matrix[,colnames(matrix) %in% a]
  gsva(m,targetSets,min.sz=5,kcdf='Gaussian') ###用counts kcdf='Poisson',else'Gaussian'
})

gsva.score.mean<-lapply(gsva.score, function(a){
  a.mmeta<-mmeta[rownames(mmeta) %in% colnames(a),]
  a<-t(a)
  by(a,a.mmeta$time.sp,colMeans)
})

gsva.score.mean.matrix<-lapply(gsva.score.mean,function(a){
  states<-names(a)
  enrichTerms<-names(a[[1]])
  x<-matrix(unlist(a),nrow=lengths(a),byrow=F)
  x<-data.frame(x)
  rownames(x)<-enrichTerms
  x$diif<-apply(x,1,function(c){
    max(c)-min(c)
  })
  colnames(x)<-c(paste0('states_',states),'diff')
  return(x)
})

b<-gsva.score.mean.matrix$GABA.15

diff.top10<-lapply(gsva.score.mean.matrix,function(a){
  a<-a[order(a$diff,decreasing = T),]
  b<-a[1:10,]
})

library(pheatmap)
library(viridis)


for(i in names(diff.top10)){
  d<-diff.top10[[i]]
  ncol<-length(colnames(d))
  pdf(paste0('gsva.top10.',i,'.pdf'),width =12 ,height =3 )
  dat<-d[,1:(ncol-1)]
  pheatmap(dat,
           cluster_cols = F,
           scale = 'row',fontsize = 8,
           color = inferno(100),
           cellwidth = 8,cellheight = 8,
           legend_breaks = c(min(dat),max(dat)),
           legend_labels = c('min','max'),
           main =paste0('Enrichment Terms in Developmental States for ',i))
  dev.off()
}

diff.top100.terms<-lapply(diff.top100,rownames)

diff.terms<-data.frame(term=unlist(diff.top100.terms),cluster=names(diff.top100.terms))

length(unique(diff.terms$term))

a<-data.frame(table(diff.terms$term))

my_col<-c('#D52A29','#ff7f0e','#176CD8','#2FA148','#8D574C')
data$time.sp<-factor(data$time.sp,levels = c('P4','P14','3M','9M','15M'))

setwd('G:\\paper1\\fig1.2\\gsva')
gene='Lncpint'
subtype='ODC.1'
VlnPlot(data,gene,idents = subtype,group.by = 'time.sp',
        pt.size = 0,cols = my_col[2:5])+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0,hjust = 0.5))+NoLegend()+
  ggtitle(subtype)+
  ylab(gene)
ggsave(filename=paste0('vln.',subtype,'.',gene,'.pdf'),width = 5,height = 3)

