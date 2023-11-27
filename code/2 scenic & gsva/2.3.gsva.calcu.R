library(stringr)
t<-paste0(unlist(regulons.enrich[1]),collapse =",")
names(t)<-NULL

names(t)<-'xx'

t<-lapply(regulons.enrich,function(a){
  paste0(unlist(a),collapse =",")
})
tt<-t(as.data.frame(t))
ttt<-data.frame(names=names(t),genes=tt)
ttt<-ttt[ttt$names %in% unique(regulons.top3$Regulon),]
write.table(ttt[1:35,],'regulons.enrich1_35.txt',sep = "\t",
            col.names = F,row.names = F,quote = F)


library(ggplot2)
library(Seurat)

meta.data<-data@meta.data
matrix<-data@assays$RNA@counts
matrix<-NormalizeData(matrix,normalization.method = 'RC')

meta.data$time.sp[meta.data$orig.ident %in% c('ctp4')]<-1
meta.data$time.sp[meta.data$orig.ident %in% c('ctp14')]<-2
meta.data$time.sp[meta.data$orig.ident %in% c('ctm3')]<-3
meta.data$time.sp[meta.data$orig.ident %in% c('ctm9')]<-4
meta.data$time.sp[meta.data$orig.ident %in% c('ctm15')]<-5
mmeta<-meta.data

sp<-do.call(rbind,tapply(rownames(mmeta),mmeta.data$subclass,function(i){
  cell<-i
  #cell<-rownames(mmeta)[mmeta$subclass=='MGL']
  cluster<-unique(mmeta$subclass[rownames(mmeta) %in% cell])
  time.sp<-mmeta$time.sp[rownames(mmeta) %in% cell]
  mt<-matrix[,colnames(matrix) %in% cell]
  mt<-t(as.matrix(mt))
  print(paste0('calulate:',cluster))
  pb<-cor(mt,as.numeric(time.sp),method = 'spearman')
  colnames(pb)<-'spearman.cor'
  pb<-data.frame(pb,cluster=cluster)
  pb<-na.omit(pb)
}))
