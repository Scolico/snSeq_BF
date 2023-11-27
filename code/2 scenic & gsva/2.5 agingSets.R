library(stringr)
aging.terms<-str_detect(names(my_genesets),pattern = 'INFLAM|INFLAM')
aging.terms.OXIDATIVE_STRESS<-names(my_genesets)[c(grep('OXIDATIVE_STRESS',
                                                        names(my_genesets)))]
aging.terms<-c(aging.terms.GROWTH_FACTOR,
               aging.terms.IMMUNE,
               aging.terms.INFLAM,
               aging.terms.INTERLEUKIN,
               aging.terms.LIPID_METABOLISM,
               aging.terms.OXIDATIVE_STRESS,
               aging.terms.OXIDATIVE_PHOSPHORYLATION,
               aging.terms.P53,
               aging.terms.RESPIRA,
               aging.terms.RIBOSOME)
x<-diff.terms[diff.terms$term %in% aging.terms,]
agingSets<-my_genesets[names(my_genesets) %in% aging.terms]

gsva.aging.matrix<-gsva.aging.res[73:108]
diff.top10<-lapply(gsva.aging.matrix,function(a){
  a<-a[order(a$diff,decreasing = T),]
  b<-a[1:10,]
})

library(pheatmap)
library(viridis)


for(i in names(diff.top10)){
  d<-diff.top10[[i]]
  ncol<-length(colnames(d))
  pdf(paste0('gsva.aging.top10.',i,'.pdf'),width =12 ,height =3 )
  dat<-d[,1:(ncol-1)]
  pheatmap(dat,
           cluster_cols = F,
           scale = 'row',fontsize = 8,
           color = inferno(100),
           cellwidth = 8,cellheight = 8,
           legend_breaks = c(min(dat),max(dat)),
           legend_labels = c('min','max'),
           main =paste0('Aging-related Terms in Developmental States for ',i))
  dev.off()
}
