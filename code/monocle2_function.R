
#Load data
seu <- subset(test.clean,subtype=='GABA.15')
monocle2_data <- function(seu,color_by = "Age"){
  library(dplyr)
  library(ggplot2)
  library(Seurat)
  library(monocle)
  library(patchwork)
  library(reshape2)
  library(pheatmap)
  library(RColorBrewer)
  library(gplots)
  library(ggsci)
  library(ggpubr)
  library(viridis)
 #Import data via seurat results
  Expression_matrix <- as(as.matrix(seu@assays$RNA@counts), 'sparseMatrix')
  cell_metadata <- seu@meta.data
  pd <- new('AnnotatedDataFrame', data = cell_metadata) 
  gene_annotation <- data.frame(gene_short_name = row.names(Expression_matrix),
                                row.names = row.names(Expression_matrix))
  fd <- new('AnnotatedDataFrame', data = gene_annotation)   
  #Build the S4 object, CellDataSet
  cds <- newCellDataSet(Expression_matrix, 
                        phenoData = pd, 
                        featureData = fd, 
                        expressionFamily = negbinomial.size())
  #Estimated size factors and dispersions (required)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  #Selection of genes
  disp_table <- dispersionTable(cds)
  unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
  diff <- differentialGeneTest(cds[unsup_clustering_genes$gene_id,],fullModelFormulaStr="~age",cores=1) #fullModelFormulaStr="~age"根据需求自定义
  deg <- subset(diff, qval < 0.01)
  deg <- deg[order(deg$qval,decreasing=F),]
  ordergene <- top_n(deg,n=1000,desc(deg$qval))
  ordergene<- row.names(ordergene)
  cds <- setOrderingFilter(cds, ordergene) 
 #Dimension reduction & Sorting
 cds <- reduceDimension(
                         cds,
                         max_components = 2,
                         method = 'DDRTree')
  cds <- orderCells(cds)
  plot1 <- plot_cell_trajectory(cds,cell_size = 1)
  plot2 <- plot_cell_trajectory(cds,cell_size = 1,color_by = color_by)+scale_color_nejm()
  return(list(cds=cds,deg=deg,plot1=plot1,plot2=plot2,ordergene=ordergene))}
  
#Identification of genes differentially expressed in the pseudotime
run_Time_diff <- function(cores = 10){ 
  ordergene <- ordergene 
  cds <- orderCells(cds,root_state = 1,reverse = F)
  plot3 <- plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 1)
  
  Time_diff <- differentialGeneTest(cds[ordergene,], cores = cores, 
                                    fullModelFormulaStr = "~sm.ns(Pseudotime)")
  Time_diff <- Time_diff[,c(5,2,3,4,6)] 
  Time_genes <- Time_diff %>% filter(qval < 1e-4) %>% pull(gene_short_name) %>% as.character()
  p = plot_pseudotime_heatmap(cds[rownames(subset(Time_diff,qval < 1e-4)),], num_clusters=5, show_rownames=F, return_heatmap=T,
                              hmcols = colorRampPalette(c(rep("royalblue",3),"white",rep("firebrick3",3)))(56))
  return(list(Time_diff = Time_diff, Time_genes = Time_genes,p = p,plot3=plot3,cds=cds))}

#Heat map of differentially expressed genes of the pseudotime
heatmap_plot <- function(p,k=4){
  clusters <- cutree(p$tree_row,k=4)
  clustering<-data.frame(clusters)
  clustering[,1] <- as.character(clustering[,1])
  colnames(clustering) <- "gene_clusters"
  cluster1 <- subset(clustering,gene_clusters == 1)
  cluster2 <- subset(clustering,gene_clusters == 2)
  cluster3 <- subset(clustering,gene_clusters == 3)
  cluster4 <- subset(clustering,gene_clusters == 4)
  pseudotime_genes <- list(cluster1=row.names(cluster1),cluster2=row.names(cluster2),cluster3=row.names(cluster3),cluster4=row.names(cluster4))
  return(list(cluster1=cluster1,cluster2=cluster2,cluster3=cluster3,cluster4=cluster4,pseudotime_genes=pseudotime_genes))}
