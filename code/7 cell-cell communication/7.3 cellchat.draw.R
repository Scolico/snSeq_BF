library(CellChat)
library(patchwork)

#1.Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#2.Compare the number of interactions and interaction strength 
#among different cell populations
#Differential number of interactions or interaction strength 
#among different cell populations

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T,measure = "weight")

#3.heatmap
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

#4.
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, top = 0.5,
                   weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", 
                                       names(object.list)[i]))
}

#5.Differential number of interactions or interaction strength 
#among different cell types

group.cellType <- c(rep("GABA_Slc5a7", 4), rep('GABA_Sox6', 4), rep("GABA_Zfhx3", 4))
group.cellType <- factor(group.cellType, levels = c("GABA_Slc5a7", "GABA_Sox6", "GABA_Zfhx3"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)

#6. Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

#7.
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "GABA_Zfhx3", signaling.exclude = NULL)
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "GABA_Slc5a7", signaling.exclude = NULL)
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "GABA_Sox6", signaling.exclude = NULL)
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "GABA_Meis2", signaling.exclude = NULL)

gg5 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "GABA_Il1rapl2", signaling.exclude = NULL)
gg6 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "GABA_Gpr149", signaling.exclude = NULL)

gg7 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "GLUT_Zfhx3", signaling.exclude = NULL)
gg8 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "GLUT_Ntng1", signaling.exclude = NULL)
gg9 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "GLUT_Dach1", signaling.exclude = NULL)

gg10 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "MGL", signaling.exclude = NULL)
gg11 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "OPC", signaling.exclude = NULL)
gg12 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "ODC", signaling.exclude = NULL)
gg13 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "ASC", signaling.exclude = NULL)

patchwork::wrap_plots(plots = list(gg1,gg2,gg3,gg4))
patchwork::wrap_plots(plots = list(gg5,gg6,gg7,gg8))
patchwork::wrap_plots(plots = list(gg9))
patchwork::wrap_plots(plots = list(gg10,gg11,gg12,gg13))

#8.Identify signaling networks with larger (or less) difference as well as signaling groups 
#based on their functional/structure similarity
#Identify signaling groups based on their functional similarity
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5,dot.alpha = 0.8,
                            color.use = c('#F3533A','#FA9F42','#8AD879','#5ACFC9'))
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)

#Identify signaling groups based on structure similarity
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5,
                            )
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

#9.Compute and visualize the pathway distance in the learned joint manifold
rankSimilarity(cellchat, type = "functional",x.rotation=1)

#10. Identify and visualize the conserved and context-specific signaling pathways
#Compare the overall information flow of each signaling pathway
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = T,axis.gap=F)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = T)
gg1 + gg2

#11. Compare outgoing (or incoming) signaling associated with each cell population
library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(
  object.list[[i]], pattern = "incoming",
  signaling = pathway.union, title = names(object.list)[i], 
  width = 5, height = 14)
ht2 = netAnalysis_signalingRole_heatmap(
  object.list[[i+1]], pattern = "incoming", signaling = pathway.union, 
  title = names(object.list)[i+1], width = 5, height = 14)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(
  object.list[[i]], pattern = "incoming", signaling = pathway.union, 
  title = names(object.list)[i], width = 5, height = 14, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(
  object.list[[i+1]], pattern = "incoming", signaling = pathway.union, 
  title = names(object.list)[i+1], width = 5, height = 14, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(
  object.list[[i]], pattern = "all", signaling = pathway.union, 
  title = names(object.list)[i], width = 5, height = 14, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(
  object.list[[i+1]], pattern = "all", signaling = pathway.union, 
  title = names(object.list)[i+1], width = 5, height =14, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#12.Identify dysfunctional signaling by comparing the communication probabities
levels(cellchat@idents[["NL"]])
netVisual_bubble(cellchat, sources.use = c(1:13), targets.use = c(5), 
                 comparison = c(1, 2), angle.x = 45,
                 signaling =c('NEGR') )
netVisual_bubble(cellchat, sources.use = c(5), targets.use = c(2,3,5,8), 
                 comparison = c(2, 1), angle.x = 45,n.colors = 100,
                 signaling =c('NRXN') )
netVisual_bubble(cellchat, sources.use = c(2:10,13), targets.use = c(5), 
                 comparison = c(2, 1), angle.x = 45,n.colors = 100,
                 signaling =c('NRXN') )

#> Comparing communications on a merged object

gg1 <- netVisual_bubble(
  cellchat, sources.use = c(1:13), targets.use = c(5),  comparison = c(1, 2), 
  max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T,
  signaling =c('NRXN','NEGR'))
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(
  cellchat, sources.use = c(1:13), targets.use = c(5),  comparison = c(1, 2), 
  max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T,
  signaling =c('NRXN','NEGR'))
#> Comparing communications on a merged object
gg1 + gg2

#13.Identify dysfunctional signaling by using differential expression analysis
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "LS"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "LS",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "NL",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(
  cellchat, pairLR.use = pairLR.use.up, 
  sources.use = 5, targets.use = c(1:13), 
  comparison = c(1, 2),  angle.x = 90, 
  remove.isolate = T,
  title.name = paste0("Up-regulated signaling in ", names(object.list)[2])
  )
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(
  cellchat, pairLR.use = pairLR.use.down, 
  sources.use = 5, targets.use = c(1:13), comparison = c(1, 2),  
  angle.x = 90, remove.isolate = T,
  title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2


# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

computeEnrichmentScore(net.down, species = 'mouse')

#14. plot pathways
pathways.show <- c("NRXN") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

# show all the significant signaling pathways from fibroblast to immune cells
pathways.show <- c("CCL") 
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(
    object.list[[i]], sources.use = c(5), targets.use = c(2:10),slot.name = "netP", 
    title.name = paste0("Signaling pathways sending from Chat - ", names(object.list)[i]), 
    legend.pos.x = 10,signaling=pathways.show)
}


#15
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) # set factor level
plotGeneExpression(cellchat, signaling = 'BMP', split.by = "datasets", colors.ggplot = T)
#> The default behaviour of split.by has changed.
#> Separate violin plots are now plotted side-by-side.
#> To restore the old behaviour of a single split violin,
#> set split.plot = TRUE.
#>       
#> This message will be shown once per session.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.