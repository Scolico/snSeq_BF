library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(igraph)
#load demo.ad as data0
data0<-demo.ad
DefaultAssay(data0)<-'RNA'
data<-data0
table(data$maintype,data$group)
maintype<-c('GABA',"Glu",'ODC','OPC','ASC','MGL')
data<-subset(data,maintype %in% c('GABA',"Glu",'ODC','OPC','ASC','MGL'))

#choose subset to analyze
data<-subset(data0,group=='GABA_Sox6')
# for maintype
#data<-subset(data0,maintype=='GABA')

table(data$orig.ident)
seurat_obj <- SetupForWGCNA(
  data,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)
# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("orig.ident","group"), # specify the columns in seurat_obj@meta.data to group by
  k = 15, # nearest-neighbors parameter
  ident.group = 'group' # set the Idents of the metacell seurat object
)
table(seurat_obj@misc[["tutorial"]][["wgcna_metacell_obj"]]@meta.data$orig.ident,seurat_obj@misc[["tutorial"]][["wgcna_metacell_obj"]]@meta.data$group)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "GABA_Sox6", # the name of the group of interest in the group.by column
  group.by='group' # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  setDatExpr = FALSE, # set this to FALSE since we did this above
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=7,
  setDatExpr=FALSE,
  minModuleSize = 30
)
PlotDendrogram(seurat_obj, main='GABA_Sox6 Dendrogram')

# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj
  #group.by.vars ='orig.ident' 
)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

seurat_obj$AD<-1
seurat_obj$AD[seurat_obj$orig.ident %in% c('ctm9','ctm15')]<-0
seurat_obj$time<-9
seurat_obj$time[seurat_obj$orig.ident %in% c('adm15','ctm15')]<-15

# list of traits to correlate
cur_traits <- c('AD', 'time')

seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits = cur_traits,
  group.by='group'
)

PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 2,
  text_digits = 2,
  text_color = 'white',
  high_color = 'yellow',
  mid_color = 'black',
  low_color = 'purple',
  #plot_max = 0.2,
  combine=F
)

mt.cor<-t(seurat_obj@misc[["tutorial"]][["mt_cor"]][["cor"]][["GABA_Sox6"]])

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'group', group_name = 'GABA_Sox6'
)

