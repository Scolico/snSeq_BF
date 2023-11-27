terms<-c(
  #'GOMF_INTERLEUKIN_1_RECEPTOR_ACTIVITY',
         'GOBP_RESPONSE_TO_INTERLEUKIN_2',
         #'GOBP_RESPONSE_TO_INTERLEUKIN_15',
         'REACTOME_PI5P_REGULATES_TP53_ACETYLATION',
         'REACTOME_CLEC7A_INFLAMMASOME_PATHWAY',
         'GOBP_POSITIVE_REGULATION_OF_RESPONSE_TO_OXIDATIVE_STRESS'
         )
subclass<-c('GABA.3','GABA.12','GABA.15','GABA.18')
t<-gsva.aging.res[1:36]
t<-t[subclass]
t<-lapply(t,function(a){
  a[rownames(a) %in% terms,]
})

library(stringr)
subclass='GABA.3'
gaba<-data.frame(t(t[[subclass]]))
age<-str_extract(rownames(gaba),pattern = '.{4,5}_')
age<-str_remove(age,'_')
gaba$timepoint<-age
gaba$timepoint<-factor(gaba$timepoint,levels = c('ctp4','ctp14','ctm3','ctm9','ctm15'))
for(i in 1:4){
  new_data <- summarySE(gaba,measurevar = colnames(gaba)[i],groupvars = c('age'))
  new_data$value<-new_data[[3]]
  new_data$group<-colnames(gaba)[i]
  new_data<-new_data[,-3]
  assign(paste0('new_data_',i),new_data)
}

new_data<-rbind(new_data_1,new_data_2,new_data_3,new_data_4)

###平均值折线
#library(Rmisc)
ggplot(new_data,aes(x=age,y=value,group=group,color=group))+ 
  geom_line(size=.6,alpha=.8)+
  geom_point(size=2,alpha=.8,shape=20) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.15,size=.6,alpha=.8) + 
  scale_color_manual(
                    values=c('#BF242A','#FFB61E','#16A951','#1685A9','#50616D'))+
  labs(x=NULL,y='avarage(GSVA Score)')+
  theme_test()+
  theme(axis.text = element_text(color = 'black'),
        axis.title.x = element_blank())+
  ylim(-0.35,0.15)+
  ggtitle(subclass)+
  scale_x_discrete(limits=c('ctp4','ctp14','ctm3','ctm9','ctm15'),
                   labels=c('P4','P14','3M','9M','15M'))+
  NoLegend()



ggsave(paste0(subclass,'_gsva.mean4.pdf'),width =3,height = 3)


ggplot(new_data,aes(x=age,y=value))+
  geom_line()
