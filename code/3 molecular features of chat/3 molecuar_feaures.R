rm(list = ls())
my_col <- c('#D52A29','#3D6CB4','#2FA148','#9467bd', '#a6d634','#efa4a5', '#c49c94','#ff7f0e', '#176CD8', '#8D574C', '#F2BF4A')

#fig3a
DimPlot(test.chat.ct,reduction = 'tsne',cols = my_col[1:3])+NoLegend()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())

#fig3b
DefaultAssay(test.chat.ct)<-'RNA'
markers<-FindAllMarkers(test.chat.ct, only.pos = T, logfc.threshold = 0.5, return.thresh = 0.01)
markers$pc<-markers$pct.1-markers$pct.2
markers.selected <- markers[markers$pc>=0.3,] %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DotPlot(test.chat.ct,features = c(unique(markers.selected$gene)))+
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_blank())+
  coord_flip()+
  scale_color_viridis_c()+
  scale_x_discrete(limits=rev(markers.selected$gene))

sgene<-c('Chat','Ngfr','Slc17a8')
DotPlot(test.chat.ct,features = sgene,cols = c('white','steelblue'))+
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_blank())
VlnPlot(test.chat.ct,features = sgene,stack = T)
VlnPlot(test.chat.ct,features='Chat')

VlnPlot(test.chat.ct, features = sgene, stack = TRUE, assay = 'RNA', group.by = 'seurat_clusters', flip = TRUE,
        pt.size = 1,fill.by = 'ident',cols = my_col[1:3],same.y.lims = T) +
  theme(axis.title = element_blank(),
    text = element_text(size = 20),
        axis.text.x = element_text(size = 20)) +
  NoLegend()

VlnPlot(test.chat.ct,'Chat')
VlnPlot(test.chat.ct,'Ngfr')
VlnPlot(test.chat.ct,'Slc17a8')
