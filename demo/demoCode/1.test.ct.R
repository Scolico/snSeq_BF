library(dplyr)
set.seed(666)
##load data
test.clean <- load("G:/seq/test.clean.rda")
##randomly select cells to construct a subset for test
cells<-rownames(sample_frac(test.clean@meta.data,0.1))
data<-subset(test.clean,cells = cells)
##test_fig1b
DimPlot(data,pt.size = .1,shuffle = T)+NoLegend()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
#test_fig1c
DimPlot(data,pt.size = .1,group.by = 'Age',shuffle = T)+
  #NoLegend()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())+
  scale_color_manual(values = c('#8D574C','#2FA148','#3D6CB4','#ff7f0e','#D52A29'))+
  ggtitle(NULL)

###fig.1f
ratio<-data.frame(table(data$orig.ident,data$spercific))
colnames(ratio)<-c('orig.ident','spercific','Freq')
ratio$orig.ident<-factor(ratio$orig.ident,levels =c('ctp4','ctp14','ctm3','ctm9','ctm15'))
ggplot(ratio,aes(x=spercific,y=Freq))+
  geom_bar(stat = 'identity',aes(fill=orig.ident),position = 'fill')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_blank())+
  scale_fill_manual(values = rev(c('#8D574C','#3D6CB4','#2FA148','#ff7f0e','#D52A29','#a6d634')))
