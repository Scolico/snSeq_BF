setwd('G:\\paper1\\fig1.2\\gsva')
library(Seurat)
library(ggplot2)
data<-test.clean
data$orig.ident<-factor(data$orig.ident,
                        levels =c('ctp4','ctp14','ctm3','ctm9','ctm15') )
mycol<-c('#D52A29','#ff7f0e','#176CD8','#2FA148','#8D574C')

gene='Lncpint'
subtype='GABA.3'
VlnPlot(data,
        features=c(gene),
        idents = subtype,
        group.by = 'orig.ident',
        pt.size = 0,cols = mycol)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0))+
  ylab(gene)+
  ggtitle(subtype)+
  scale_x_discrete(labels=c('P4','P14','3M','9M','15M'))
ggsave(paste0(gene,'.',subtype,'.pdf'),width = 6,height = 4)

