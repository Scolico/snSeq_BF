glut<-subset(test.clean,subtype %in% c('GLUT.1','GLUT.2','GLUT.3','GLUT.4'))
glut$subtype<-factor(glut$subtype,levels = c('GLUT.1','GLUT.2','GLUT.3','GLUT.4'))
Idents(glut)<-'subtype'
VlnPlot(glut,features=c('Slc17a6','Slc17a7'),pt.size = 0,stack = T,
        fill.by = 'ident',same.y.lims = T,flip = T)+
  theme(axis.title.x = element_blank())+
  scale_fill_manual(values = c('#FFB61E','#FF4777','#FF461F','#1685A9'))+
  NoLegend()

DimPlot(test.clean,pt.size = .1,shuffle = T)+NoLegend()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
ggsave('umap2.pdf',height = 8,width = 8)

DimPlot(test.clean,pt.size = .1,group.by = 'orig.ident',shuffle = T)+
  #NoLegend()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())+
  scale_color_manual(values = c('#8D574C','#2FA148','#3D6CB4','#ff7f0e','#D52A29'))+
  ggtitle(NULL)
ggsave('umap.age3.pdf',height = 8,width = 8)


FeaturePlot(gaba,'Bmp4')

gaba<-subset(test.clean,subclass %in% c("GABA_Zfhx3","GABA_Meis2","GABA_Gpr149","GABA_Sox6","GABA_Il1rapl2"))
Idents(gaba)<-'subclass'
DefaultAssay(gaba)<-'integrated'
gaba.marker<-FindAllMarkers(gaba,only.pos = T)

DimPlot(test.clean,label = T)+NoLegend()
FeaturePlot(test.clean,'Tac1')
VlnPlot(test.clean,'Sst')
FeaturePlot(test.clean,'Npy')
FeaturePlot(test.clean,'Nos1')
FeaturePlot(test.clean,'Tacr3')
FeaturePlot(test.clean,'Kcnb2')
FeaturePlot(test.clean,'Sst')
FeaturePlot(test.clean,'Calb1')
FeaturePlot(test.clean,'Gpr149')

test.clean.markers_selected$diff<-test.clean.markers_selected$pct.1 - test.clean.markers_selected$pct.2

t<-data.frame(table(test.clean.markers_selected$gene))
uni<-t$Var1[t$Freq==1]
tt<-test.clean.markers_selected[test.clean.markers_selected$gene %in% uni,]

sub<-subset(test.clean,subtype %in% c('GABA.3','GABA.12'))
DimPlot(sub)
FeaturePlot(sub,'Pitx3')
VlnPlot(test.clean,'Drd1')

VlnPlot(test.clean,'Tac1',idents = 'GABA.12',group.by = 'orig.ident')
###fig.1f
ratio<-data.frame(table(test.clean$orig.ident,test.clean$spercific))
colnames(ratio)<-c('orig.ident','spercific','Freq')
ratio$orig.ident<-factor(ratio$orig.ident,levels =c('ctp4','ctp14','ctm3','ctm9','ctm15'))
ggplot(ratio,aes(x=spercific,y=Freq))+
  geom_bar(stat = 'identity',aes(fill=orig.ident),position = 'fill')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_blank())+
  scale_fill_manual(values = rev(c('#8D574C','#3D6CB4','#2FA148','#ff7f0e','#D52A29')))
