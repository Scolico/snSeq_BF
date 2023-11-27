DotPlot(test, features = c('Rbfox3','C1qb','Slc6a13','Slc1a3','Aqp4','Mbp','Pdgfra','Cldn5','Tmem212','Slc17a6','Gad2','Snap25'), assay = 'RNA')

mainclass <- setNames(as.character(test$seurat_clusters), colnames(test))
mainclass[mainclass %in% c(10,24,29,36)] <- 'Glutamatergic neurons'
mainclass[mainclass %in% c(3:5,8,9,12:15,17,19:21,23,25:28,30)] <- 'GABAergic neurons'
mainclass[mainclass %in% c(1,2,6,7,11,16,18,22,31:35,37:39)] <- 'Non_neu'

test$mainclass <- mainclass
DimPlot(test, group.by = 'mainclass')

test.neuron.ct <- subset(test, mainclass %in% c('GABAergic neurons', 'Glutamatergic neurons'))

subtype <- setNames(as.character(test$seurat_clusters), colnames(test))
subtype[subtype %in% c(14,18,22)] <- 'Glutamatergic neurons'
subtype[subtype %in% c(3:8,11:12,15,17,19,20)] <- 'GABAergic neurons'
subtype[subtype == 1] <- 'ASC.1'
subtype[subtype == 11] <- 'ASC.2'
subtype[subtype == 37] <- 'ASC.3'
subtype[subtype == 7] <- 'OPC.1'
subtype[subtype == 36] <- 'OPC.2'
subtype[subtype == 2] <- 'ODC.1'
subtype[subtype == 6] <- 'ODC.2'
subtype[subtype == 33] <- 'COP' # committed oligodendrocyte precurser
subtype[subtype == 31] <- 'ENDO'
subtype[subtype == 32] <- 'FB'
subtype[subtype == 35] <- 'PER'
subtype[subtype == 18] <- 'MGL'
subtype[subtype == 10] <- 'GLUT.1'
subtype[subtype == 24] <- 'GLUT.2'
subtype[subtype == 29] <- 'GLUT.3'
subtype[subtype == 36] <- 'GLUT.4'
subtype[subtype == 3] <- 'GLUT.1'
subtype[subtype == 26] <- 'Ependymal'
subtype[subtype %in% c(16,22)] <- 'UkN'

library(plyr)


curid <- c(1:39)
test@meta.data$subtype <- plyr::mapvalues(test$seurat_clusters, from = curid, to = c('ASC.1',
                                                                             'ODC.1',
                                                                             'GABA.1','GABA.2','GABA.3',
                                                                             'ODC.2',
                                                                             'OPC.1',
                                                                             'GABA.4','GABA.5',
                                                                             'GLUT.1',
                                                                             'ASC.2',
                                                                             'GABA.6','GABA.7','GABA.8','GABA.9',
                                                                             'UKN',
                                                                             'GABA.10',
                                                                             'MGL',
                                                                             'GABA.11','GABA.12','GABA.13',
                                                                             'UKN',
                                                                             'GABA.14',
                                                                             'GLUT.2',
                                                                             'GABA.15','GABA.16','GABA.17','GABA.18',
                                                                             'GLUT.3',
                                                                             'GABA.19',
                                                                             'VEC',
                                                                             'VLMC',
                                                                             'COP',
                                                                             'UKN',
                                                                             'PER',
                                                                             'GLUT.4',
                                                                             'EPEN',
                                                                             'OPC.2',
                                                                             'ASC.3'))



library(scrattch.hicat)
library(dendextend)
library(dplyr)
library(matrixStats)
library(Matrix)

norm.dat <- test@assays$integrated@data
norm.dat <- Matrix(norm.dat, sparse = T)
Idents(test) <- test$subtype

test.subtype.markers <- FindAllMarkers(test, only.pos = T, logfc.threshold = 0.5, return.thresh = 0.01)
test.subtype.markers_selected <- test.subtype.markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)

subtype.med <- get_cl_medians(norm.dat[unique(test.subtype.markers_selected$gene),], 
                         droplevels(test$seurat_clusters))
subtype.dend.result <- build_dend(subtype.med, 
                          nboot = 1000, ncores = 20)

plot(subtype.dend.result$dend)

dend <- dend.result$dend
###attach cluster labels to the leaves of the tree 
labels(dend.labeled) <- consensus.cl.df[labels(dend), "cluster_label"]
plot(dend) 

test$seurat_clusters <- factor(test$seurat_clusters, levels = dend.result$pvclust.result$hclust$order)

test@meta.data$subclass <- plyr::mapvalues(test$seurat_clusters, from = curid, to = c('ASC','ODC','GABA',
                                                                                     'UKN','UKN',
                                                                                     'GABA_Gpr149','GABA_Gpr149',
                                                                                     'GABA_Meis2','GABA_Meis2','GABA_Meis2','GABA_Meis2',
                                                                                     'GABA_SLc5a7',
                                                                                     'GLUT.1',
                                                                                     'GABA_Sox6','GABA_Sox6','GABA_Sox6','GABA_Sox6','GABA_Sox6',
                                                                                     'GABA.Nfib','GABA.7','GABA.8','GABA.9',
                                                                                     'UKN.1',
                                                                                     'GABA.10',
                                                                                     'MGL',
                                                                                     'GABA.11','GABA.12','GABA.13',
                                                                                     'UKN.2',
                                                                                     'GABA.14',
                                                                                     'GLUT.2',
                                                                                     'GABA.15','GABA.16','GABA.17','GABA.18',
                                                                                     'GLUT.3',
                                                                                     'GABA.19',
                                                                                     'VEC',
                                                                                     'VLMC',
                                                                                     'COP',
                                                                                     'UKN.3',
                                                                                     'PER',
                                                                                     'GLUT.4',
                                                                                     'EPEN',
                                                                                     'OPC.2',
                                                                                     'ASC.3'))

test.clean <- subset(test, subtype != 'UKN')
Idents(test.clean) <- 'subtype'
test.clean$subtype<-factor(test.clean$subtype,levels = c(paste0('GABA.',1:19),
                                                         paste0('GLUT.',1:4),
                                                         paste0('ASC.',1:3),
                                                         paste0('ODC.',1:2),
                                                         paste0('OPC.',1:2),
                                                         'COP','MGL','EPEN','PER','VLMC','VEC'))
test.clean.markers <- FindAllMarkers(test.clean, only.pos = T, logfc.threshold = 0.5, return.thresh = 0.01)
test.clean.markers_selected <- test.clean.markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)
dot.markers<- test.clean.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
dot.markers<-dot.markers$gene

norm.dat <- test.clean@assays$integrated@data
norm.dat <- Matrix(norm.dat, sparse = T)

clean.med <- get_cl_medians(norm.dat[unique(test.clean.markers_selected$gene),], 
                            droplevels(test.clean$subtype))
clean.dend.result <- build_dend(clean.med, 
                                nboot = 1000, ncores = 20)
plot(clean.dend.result$dend)

Idents(test.clean) <- 'subclass'
test.clean.subclass.markers <- FindAllMarkers(test.clean, only.pos = T, logfc.threshold = 0.5, return.thresh = 0.01)
test.clean.subclass.markers_selected <- test.clean.subclass.markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)


norm.dat <- test.clean@assays$integrated@data
norm.dat <- Matrix(norm.dat, sparse = T)

clean.subclass.med <- get_cl_medians(norm.dat[unique(test.clean.subclass.markers_selected$gene),], 
                              droplevels(test.clean$subclass))
clean.subclass.dend.result <- build_dend(clean.subclass.med, 
                                  nboot = 1000, ncores = 20)
plot(clean.subclass.dend.result$dend)

test.clean.dend_marker = list()
test.clean.dend_marker$cl.28_19vs36_9 <- FindMarkers(test.clean, ident.1 = c(28,8,27,17,26,13,19), ident.2 = c(36,24,12,21,23,10,3,9), logfc.threshold = 0.5,return.thresh = 0.05)
test.clean.dend_marker$cl.12_21vs23_9 <- FindMarkers(test.clean, ident.1 = c(12,21), ident.2 = c(23,10,3,9), logfc.threshold = 0.5,return.thresh = 0.05)

cl.28_19vs36_9_left <- test.clean.dend_marker$cl.28_19vs36_9  %>%   
  filter(abs(pct.1-pct.2) > 0.3 & avg_log2FC >= 1) %>%
  top_n(50, avg_log2FC) %>% arrange(desc(avg_log2FC))
write.csv(cl.28_19vs36_9_left, file = '/mnt/snRNA_seq/snrna/my_fig/altas/cl.28_19vs36_9_left.csv')

cl.28_19vs36_9_right <- test.clean.dend_marker$cl.28_19vs36_9  %>%   
  filter(abs(pct.1-pct.2) > 0.3 & avg_log2FC <= -1) %>%
  top_n(50, avg_log2FC) %>% arrange(desc(-avg_log2FC))
write.csv(cl.28_19vs36_9_right, file = '/mnt/snRNA_seq/snrna/my_fig/altas/cl.28_19vs36_9_right.csv')


test.clean$seurat_clusters <- droplevels(test.clean$seurat_clusters)
test.clean@meta.data$subclass <- plyr::mapvalues(test.clean$seurat_clusters, 
                                                 from = levels(test.clean$seurat_clusters)[clean.dend.result$pvclust.result$hclust$order], 
                                                 to = c('Vascular','Vascular','Vascular','MGL','EPEN','ASC','ASC','ASC','ODC','ODC','OPC','OPC','OPC',
                                                        'GABA_Gpr149','GABA_Gpr149',
                                                        'GABA_Meis2','GABA_Meis2','GABA_Meis2','GABA_Meis2',
                                                        'GABA_Slc5a7',
                                                        'GLUT_Dach1',
                                                        'GABA_Sox6','GABA_Sox6','GABA_Sox6','GABA_Sox6','GABA_Sox6','GABA_Sox6',
                                                        'GABA_Il1rapl2',
                                                        'GLUT_Sv2b','GLUT_Ntng1',
                                                        'GABA_Zfhx3','GABA_Zfhx3','GABA_Zfhx3',
                                                        'GLUT_Zfhx3',
                                                        'GABA_Zfhx3','GABA_Zfhx3'))
Idents(test.clean) <- 'subclass'
test.subclass.markers <- FindAllMarkers(test.clean, only.pos = T, logfc.threshold = 0.5, return.thresh = 0.01)


#####DRAW############
DefaultAssay(test.clean)<-'integrated'
DefaultAssay(test.clean)<-'RNA'
dot.mk<-c('Cldn5','Cped1','Atp13a5','Cx3cr1','Tmem212','Inpp5d','Itih3','Slc1a3','Ptgds','Mog','Bcas1','Pdgfra',"Top2a",
          'Gpr149','Rarb','Meis2','Sphkap','Fstl4','Htr2c','Chat','Dach1','Sox6','Cntnap5c','Zfp804a',
          'Adarb2','Nell1','Gpc6','Reln','Trhde','Il1rapl2','Slc17a7','Sv2b','Slc17a6','Ntng1','Zfhx3','Rorb','Nwd2','Meis1','Vwc2l')
dot.mk2<-c('Cldn5','Cped1','Atp13a5','Cx3cr1','Tmem212','Inpp5d','Itih3','Slc1a3','Ptgds','Mog','Bcas1','Pdgfra',"Top2a",
           'Gpr149','Syt6','Foxp2','Meis2','Zeb2','Tshzl','Rarb','Uncl3c','Lhx8','Chat','Zfhx3','Trpc4',
           'Dpf3',"Nxph1",'Kcnh7','Fign','Ebf1','Sox6','Rorb','Nwd2','Meis1','Vwc2l','Reln','Trhde')
dot.mk3<-c('Cldn5','Cped1','Atp13a5','Cx3cr1','Tmem212','Inpp5d','Itih3','Slc1a3','Ptgds','Mog','Bcas1','Pdgfra',"Top2a",
          'Gpr149','Pbx3','Tshz1','Rarb','Fstl4','Meis2','Chat','Dach1','Cntnap5c','Casz1','Nwd2','Vwc2l',
          'Gpc6','Cdh4','Il1rapl2','Sv2b','Ntng1','Trpc4','Trhde','Cntnap5b','Rmst','Zfhx3','Kcnh7')

#fig1e
DotPlot(test.clean,features = c(unique(dot.mk3)),cols = c('white','steelblue'),col.min = 0)+
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_blank())
+
  NoLegend()+
  scale_color_gradient2(low='white',mid='seagreen',high='gold')

features<-c('Rbfox3','Gad2','Slc17a6',
            'Slc1a3','Mog','Bcas1','Pdgfra',
            'Cx3cr1','Cldn5','Atp13a5','Slc6a13','Tmem212')


VlnPlot(test.clean,features = c(unique(dot.mk3)),stack = T)

for(feature in features){
  FeaturePlot(test.clean,features = feature)+
    NoLegend()+
    theme(axis.title = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())+
    ggtitle(NULL)
  ggsave(filename = paste0('/media/zkmatrix/ESD-USB/',feature,'.umap.pdf'),height = 5,width = 5,dpi = 300)
}

#fig1d
plot(clean.subclass.dend.result$dend)

#fig1b,1c
DimPlot(test.clean,label = T,reduction = 'tsne')+NoLegend()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
my_cor<-c('#885748','#309D42','#356AAC','#EA7E26','#C72F22')
DimPlot(test.clean,reduction = 'tsne',group.by = 'orig.ident',cols = my_cor)+
  NoLegend()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())+
  ggtitle(NULL)
