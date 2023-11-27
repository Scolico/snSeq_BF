library(msigdbr)
library(GSVA)
library(Seurat)
library(stringr)
library(Matrix)
library(limma)
library(pheatmap)
data<-demo.ad
###construct a  reference enrichment term set ###
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

###construct a DEGs list
my_genelist<-read.table('chat.degs.txt')
colnames(my_genelist)<-'gene'
my_genelist<-c(my_genelist$gene)

###select DEGs-related terms
library(tidyverse)
a<-unlist(lapply(my_genesets,function(a){
  isFALSE(is_empty(intersect(a,my_genelist)))
}))
targetSets<-my_genesets[a]

###select 100 terms for test
targetSets<-targetSets[1:100]

###caculat gsva score
# time ~ 15 min
matrix<-data@assays$RNA@counts
matrix<-NormalizeData(matrix,normalization.method='RC',scale.factor=1e6)
matrix<-Matrix(log2(matrix+1),sparse=T)
gsva.score<-gsva(matrix,targetSets,min.sz=5,kcdf='Gaussian') ###用counts kcdf='Poisson',else'Gaussian'

meta<-data@meta.data
meta<-meta[,c("orig.ident","group")]
meta$phetp<-'AD'
meta$phetp[meta$orig.ident %in% c('ctm15','ctm9')]<-'Norm'


meta.1<-meta[meta$group=='GABA_Sox6',]

gsva.score.1<-gsva.score[,rownames(meta.1)]

####difference analysis via limma
logFCcutoff<-log2(0.1)
adjPvalueCutoff<-0.05

phetp<-factor(meta.1$phetp,levels = c('AD','Norm'))
design<-model.matrix(~0+phetp)
colnames(design)<-c('AD','Norm')

fit<-lmFit(gsva.score.1,design)
fit<-eBayes(fit)
allGeneSets<-topTable(fit,coef = 'AD',number = Inf)
DEgeneSets<-topTable(fit,coef = 'AD',number = Inf, p.value = adjPvalueCutoff,adjust.method = 'BH')
gsva.res<-decideTests(fit,p.value = adjPvalueCutoff)
res<-as.data.frame(summary(gsva.res))

###calculate average GSVA score
cell<-tapply(rownames(meta.1),list(meta.1$phetp),list)
cell.n<-data.frame(table(meta.1$phetp))
cell.n$name<-paste0(cell.n$Var1,'/',cell.n$Var2)
names(cell)<-cell.n$name
spc.gsva<-gsva.score[rownames(gsva.score) %in% DEgs,]
tt <- sapply(cell,function(b) rowMeans(gsva.score.1[,b])) ###按分组信息求平均值

tt<-data.frame(tt)
tt$diff<-tt$AD-tt$Norm
tt<-tt[order(abs(tt$diff),decreasing = T),]


pheatmap(tt[1:20,1:2],cluster_cols = F)