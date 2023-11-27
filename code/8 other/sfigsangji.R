install.packages("treemapify")
library(ggplot2)
library(treemapify)
draw.data<-data.frame(table(draw$group,draw$subtype,draw$integrated_snn_res.0.8))
draw.data<-draw.data[draw.data$Freq!=0,]
ggplot(draw.data, aes(area = Freq, ,label=Var2,fill=Var3,
                      subgroup = Var1,subgroup2=Var2,subgroup3=Var3)) +
  geom_treemap()+
  geom_treemap_text(fontface = "italic", colour = "white", place = "centre",
                    grow = TRUE)

colnames(draw.data)<-c("group","subtype","integrated_snn_res.0.8",'Freq')
draw.data$group<-factor(draw.data$group,levels = rev(c(
  'GABA_Zfhx3','GABA_Gpr149','GABA_Sox6','GABA_Il1rapl2','GABA_Slc5a7','GABA_Meis2',
  'GLUT_Zfhx3','GLUT_Ntng1','GLUT_Dach1',
  'ODC','OPC','ASC','MGL')))

draw.data$subtype<-factor(draw.data$subtype,levels = rev(c(
  'GABA.1',
  'GABA.13',
  'GABA.14',
  'GABA.5',
  'GABA.6',
  'GABA.3',
  'GABA.12',
  'GABA.4',
  'GABA.7',
  'GABA.10',
  'GABA.16',
  'GABA.17',
  'GABA.18',
  'GABA.11',
  'GABA.15',
  'GABA.2',
  'GABA.19',
  'GABA.9',
  'GABA.8',
  'GLUT.1',
  'GLUT.2',
  'GLUT.3',
  'ODC.1',
  'ODC.2',
  'OPC.1',
  'OPC.2',
  'MGL',
  'COP',
  'ASC.1',
  'ASC.2',
  'ASC.3'
  
)))

draw.data$integrated_snn_res.0.8<-factor(draw.data$integrated_snn_res.0.8,levels = rev(c(
  29,
  10,
  6,
  25,
  26,
  17,
  2,
  7,
  9,
  16,
  21,
  13,
  14,
  15,
  19,
  18,
  23,
  4,
  12,
  20,
  30,
  24,
  27,
  1,
  5,
  3,
  11,
  0,
  22,
  8
)))

ggplot(data = draw.data,
       aes(axis1 = group, axis2 =subtype , axis3 = integrated_snn_res.0.8 ,
           y = Freq)) +
  #scale_x_discrete(limits = c("group","subtype","integrated_snn_res.0.8"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = subtype)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  NoLegend()+
  ggtitle("The affiliation of Norm cells after re_clustering integrated with AD cells")
ggsave('gg0.8.pdf',plot = gg0.8,width = 8,height = 10)