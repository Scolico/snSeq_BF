##R version 4.1.2 and monocle 2.22.0 were used for analysis
library(dplyr)
set.seed(666)
##load data
test.clean <- load("G:/seq/test.clean.rda")
##randomly select cells to construct a subset for test
cells<-rownames(sample_frac(test.clean@meta.data,0.1))
data<-subset(test.clean,cells = cells)
#Load data
monocle2_data <- function(data,color_by = "Age"){
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
  Expression_matrix <- as(as.matrix(data@assays$RNA@counts), 'sparseMatrix')
  cell_metadata <- data@meta.data
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
  diff <- differentialGeneTest(cds[unsup_clustering_genes$gene_id,],fullModelFormulaStr="~Age",cores=1) #fullModelFormulaStr="~age"根据需求自定义
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
run_Time_diff <- function(ordergene,cds,cores = 10){ 
  ordergene <- ordergene 
  cds <- orderCells(cds,root_state = 1,reverse = F)
  plot3 <- plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 1)
  
  Time_diff <- differentialGeneTest(cds[ordergene,], cores = cores, 
                                    fullModelFormulaStr = "~sm.ns(Pseudotime)")
  Time_diff <- Time_diff[,c(5,2,3,4,6)] 
  Time_genes <- Time_diff %>% filter(qval < 1e-4) %>% pull(gene_short_name) %>% as.character()
  p = plot_pseudotime_heatmap(cds[rownames(subset(Time_diff,qval < 1e-4)),], num_clusters=5, show_rownames=F, return_heatmap=T)
  return(list(Time_diff = Time_diff, Time_genes = Time_genes,p = p,plot3=plot3,cds=cds))}
## Run GeneSwitches
library(GeneSwitches)
library(monocle)
library(Seurat)
library(SingleCellExperiment)
library(homologene)

gene <- homologene(gs_genelists$genenames, inTax = 9606, outTax = 10090)
my_genelists <- gs_genelists[gs_genelists$genenames %in% gene$`9606`,]
t <- 0
j = 1
for (i in my_genelists$genenames){
  t[j] <- gene[gene$`9606` == i, '10090']
  j <- j+1
}
my_genelists[['genenames']] <- t

#Input dataset

cds <- monocle2_data_test$cds
logexpr <- as.matrix(log2(exprs(cds) + 1))
#Convert from trajectory results
plot_monocle_State(cds)

sce_p1 <- convert_monocle2(cds, states = c(1,2), expdata = logexpr)

sce_p2 <- convert_monocle2(cds, states = c(1,3), expdata = logexpr)
#Binarize gene expression
sce_p1 <- binarize_exp(sce_p1, ncores = 10)
#Fit logistic regression & estimate switching time
sce_p1 <- find_switch_logistic_fastglm(sce_p1, downsample = T)
## filter top 15 best fitting switching genes among all the genes
p1.sg_allgenes <- filter_switchgenes(sce_p1, allgenes = T)
## filter top 10 best fitting switching genes among surface proteins and TFs only
p1.sg_surface <- filter_switchgenes(sce_p1, allgenes = F, genelists = my_genelists, genetype = c('Surface proteins'), topnum = 10)
p1.sg_tf <- filter_switchgenes(sce_p1, allgenes = F, genelists = my_genelists, genetype = c('TFs'), topnum = 10)
## combine switching genes and remove duplicated genes from sg allgenes
p1.sg_vis <- rbind(p1.sg_surface, p1.sg_allgenes[setdiff(rownames(p1.sg_allgenes), rownames(p1.sg_surface)),])
p1.sg_vis <- rbind(p1.sg_vis, p1.sg_tf[setdiff(rownames(p1.sg_tf), rownames(p1.sg_vis)),])
plot_timeline_ggplot(p1.sg_vis, timedata = sce_p1$Pseudotime, txtsize = 5)
ggsave('geneswitch.p1.pdf', width = 6, height = 6)

sce_p2 <- binarize_exp(sce_p2, ncores = 1)
sce_p2 <- find_switch_logistic_fastglm(sce_p2, downsample = T)
p2.sg_allgenes <- filter_switchgenes(sce_p2, allgenes = T, topnum = 10)
p2.sg_surface <- filter_switchgenes(sce_p2, allgenes = F, genelists = my_genelists, genetype = c('Surface proteins'), topnum = 10)
p2.sg_tf <- filter_switchgenes(sce_p2, allgenes = F, genelists = my_genelists, genetype = c('TFs'), topnum = 10)

p2.sg_vis <- rbind(p2.sg_surface, p2.sg_allgenes[setdiff(rownames(p2.sg_allgenes), rownames(p2.sg_surface)),])
p2.sg_vis <- rbind(p2.sg_vis, p2.sg_tf[setdiff(rownames(p2.sg_tf), rownames(p2.sg_vis)),])
plot_timeline_ggplot(p2.sg_vis, timedata = sce_p2$Pseudotime, txtsize = 5)



##### compare two trajectoy
sg.p1_p2 <- common_genes(p1.sg_vis, p2.sg_vis, r2cutoff = 0.1, path1name = 'branch1', path2name = 'branch2')
common_genes_plot(sg.p1_p2, timedata = sce_p1$Pseudotime)

sg_disgs.p1_p2 <- distinct_genes(p1.sg_vis, p2.sg_vis, r2cutoff = 0.1, path1name = 'branch1', path2name = 'branch2',
                                 path1time = sce_p1$Pseudotime, path2time = sce_p2$Pseudotime)
plot_timeline_ggplot(sg_disgs.p1_p2, timedata = sce_p1$Pseudotime, color_by = 'Paths', iffulltml = F, txtsize = 5)


