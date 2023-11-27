dim(regulons_incidMat)
library(AUCell)
regulonAucThresholds[1:50]
dim(cellAUC)
rownames(cellAUC)[1:100]
regulonAucThresholds
BiocManager::install('SCENIC')
max(regulons_incidMat)

if (!requireNamespace("devtools", quietly=TRUE))install.packages("devtools")

devtools::install_github("aertslab/SCENIC")

packageVersion("SCENIC")

library(SCENIC)


a<-regulons[lengths(regulons)>=50]
rss.singleCell<-calcRSS(cellAUC,cellAnnotation = seu$subtype)
rss.time<-calcRSS(cellAUC,cellAnnotation = seu$Age)


head(as.data.frame(regulonAucThresholds))

regulonAucThresholds<-regulonAucThreshold

names(regulonAucThresholds)<-regulon_names
thresholds<-regulonAucThresholds[names(regulonAucThresholds) %in% rownames(cellAUC)]

binarizeAUC <- function(auc, thresholds)
{
  thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <-as.numeric(thresholds[x])
                                     names(which(auc[x,]>trh))
                                   }),names(thresholds))
  
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"  
  
  return(binaryRegulonActivity)
}



binaryRegulonActivity <- binarizeAUC(cellAUC, thresholds)
dim(binaryRegulonActivity)

BiocManager::install('pheatmap')
library(pheatmap)

cellClusters<-seu@meta.data

selectedResolution <- "subtype" # select resolution
# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellClusters), cellClusters[,selectedResolution]) 
#regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(cellAUC[,cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

regulonBinary_byCellType <- sapply(cellsPerCluster,
                                   function(cells) rowMeans(binaryRegulonActivity[,cells]))
# Scale expression:
regulonBinary_byCellType_Scaled <- t(scale(t(regulonBinary_byCellType), center = T, scale=T))


#pheatmap(regulonBinary_byCellType,scale = 'none')
pheatmap(regulonBinary_byCellType_Scaled[rownames(regulonBinary_byCellType_Scaled) %in% unique(regulons.top3$Regulon),],
         scale = 'none',cluster_cols = F,gaps_col = c(13),color = mako(100))

#pheatmap(regulonActivity_byCellType,scale = 'none')
pheatmap(regulonActivity_byCellType_Scaled[rownames(regulonActivity_byCellType_Scaled) %in% unique(regulons.top3$Regulon),],
         scale = 'none',cluster_cols = F,gaps_col = c(13),color = turbo(100))

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
dim(topRegulators)

regulons.top3<-topRegulators %>% group_by(CellType) %>% top_n(3, RelativeActivity)

library(ggplot2)
ggplot(regulons.top3,aes(y=Regulon,x=CellType,fill=RelativeActivity))+
  geom_point(aes(size=RelativeActivity,color=RelativeActivity))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  scale_x_discrete(limits=regulons.top3$CellType)+
  scale_y_discrete(limits=unique(regulons.top3$Regulon))

selectedResolution <- "Age" # select resolution
# Split the cells by cluster:
cellsPerCluster <- split(rownames(cellClusters), cellClusters[,selectedResolution]) 
#regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
# Calculate average expression:
regulonActivity_byAge <- sapply(cellsPerCluster,
                                function(cells) rowMeans(cellAUC[,cells]))
# Scale expression:
regulonActivity_byAge_Scaled <- t(scale(t(regulonActivity_byAge), center = T, scale=T))

regulonBinary_byAge <- sapply(cellsPerCluster,
                              function(cells) rowMeans(binaryRegulonActivity[,cells]))
# Scale expression:
regulonBinary_byAge_Scaled <- t(scale(t(regulonBinary_byAge), center = T, scale=T))


#pheatmap(regulonBinary_byCellType,scale = 'none')
pheatmap(regulonBinary_byAge,scale = 'row',cluster_cols = F,color = viridis(100))

#pheatmap(regulonActivity_byCellType_Scaled,scale = 'none')
pheatmap(regulonActivity_byAge,scale = 'row',cluster_cols = F,color = viridis(100))


rss.singleCell<-calcRSS(cellAUC,cellAnnotation = seu$subtype)
#rss.time<-calcRSS(cellAUC,cellAnnotation = seu$Age)

rssPlot.single <- plotRSS(rss.singleCell)
plotly::ggplotly(rssPlot.single$plot)
plotRSS_oneSet(rss.singleCell, setName = "GABA.13")

tfsToPlot <- c("Gli1") 
regulonsToPlot <- unlist(lapply(tfsToPlot, function(x) grep(paste0("^", x,"_"), rownames(cellAUC), value=TRUE)))

options(repr.plot.width=10, repr.plot.height=8) # To set the figure size in Jupyter
par(mfrow=c(2,3))
# Plot expression:
selectedEmbedding <- embeddings[["tsne"]]


AUCell_plotTSNE(selectedEmbedding, exprMat_log[tfsToPlot,], plots=c("Expression"), cex = .5)
# Plot regulon activity:
AUCell_plotTSNE(selectedEmbedding, cellAUC[regulonsToPlot,], plots=c("AUC"), cex = .5)


t.binary<-t(binaryRegulonActivity)
t.activity<-t(cellAUC)

regulons.top2<-regulons.top3 %>% group_by(CellType) %>% top_n(2, RelativeActivity)
my_col<-rev(rainbow(36))
my_col<-rep(my_col,each=2)
regulons.top2$col<-my_col
for(i in regulons.top2$Regulon){
  regulon<-i
  col<-regulons.top2$col[regulons.top2$Regulon==i]
  #seu.df$regulon  <- t.activity[,regulon]
  seu.df$regulon  <- t.binary[,regulon]
  data.x<-data.frame(seu.df, embeddings$umap@cell.embeddings)
  ggplot(data.x, aes(UMAP_1, UMAP_2,color=factor(regulon)))  +
    geom_point(size=.1,alpha=0.8) + 
    scale_color_manual(values = c('lightgrey',col))+
    theme(legend.position = "none") + theme_bw()+
    NoLegend()+
    ggtitle(regulon)
  ggsave(paste0(regulon,'.umap.pdf'),width = 8,height = 8)
}



####°´Ê±¼ä
tt<-data.frame(table(seu$Age,seu$spercific))
cellgroups1<-tapply(colnames(cellAUC),seu$subtype,list)

tt$index<-1:180
names(cellgroups)<-1:180

rss.time<-lapply(cellgroups1,function(x){
  print(paste0('start to analysis',names(x)))
  rss.time<-calcRSS(cellAUC[,c(colnames(cellAUC) %in% unlist(x))],cellAnnotation = seu$Age[rownames(seu@meta.data) %in% unlist(x)])
})
cell.exclude<-unique(names(cellgroups)[lengths(cellgroups)==1])
for(CellType in setdiff(names(rss.time),cell.exclude)){
  dd<-data.frame(rss.time[CellType])
  dd<-dd[,c(5,4,2,3,1)]
  dd<-t(scale(t(dd), center = T, scale=T))
  #pheatmap(dd,scale = 'none',cluster_cols = F,color = viridis(100))
  topReg.time <- reshape2::melt(as.matrix(dd))
  colnames(topReg.time) <- c('Regulon',"time", "RelativeActivity")
  topReg.time$time <- factor(as.character(topReg.time$time))
  topReg.time <- topReg.time[which(topReg.time$RelativeActivity>0),]
  dim(topReg.time)
  top.time<-topReg.time %>% group_by(time) %>% top_n(5,RelativeActivity)
  pdf(paste0(CellType,'.time.pdf'), width = 3.4, height = 4.46)
  ddd<-regulonActivity_byAge_Scaled[rownames(regulonActivity_byAge_Scaled) %in% unique(top.time$Regulon),]
  colnames(ddd)<-c('P4','P14','3M','9M','15M')
  pheatmap(ddd,
           scale = 'none',cluster_cols = F,color = inferno(100),angle_col = 90)
  dev.off()
}


