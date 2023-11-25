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

#Input datasets
load('/mnt/snRNA_seq/pseudo.rda')

cds <- pseudo.result$cds
logexpr <- as.matrix(log2(exprs(cds) + 1))
#Convert from trajectory results
plot_monocle_State(pseudo.result$cds)

sce_p1 <- convert_monocle2(cds, states = c(1,2,3), expdata = logexpr)

sce_p2 <- convert_monocle2(cds, states = c(1,2,4), expdata = logexpr)
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

sce_p2 <- binarize_exp(sce_p2, ncores = 10)
sce_p2 <- find_switch_logistic_fastglm(sce_p2, downsample = T)
p2.sg_allgenes <- filter_switchgenes(sce_p2, allgenes = T, topnum = 10)
p2.sg_surface <- filter_switchgenes(sce_p2, allgenes = F, genelists = my_genelists, genetype = c('Surface proteins'), topnum = 10)
p2.sg_tf <- filter_switchgenes(sce_p2, allgenes = F, genelists = my_genelists, genetype = c('TFs'), topnum = 10)

p2.sg_vis <- rbind(p2.sg_surface, p2.sg_allgenes[setdiff(rownames(p2.sg_allgenes), rownames(p2.sg_surface)),])
p2.sg_vis <- rbind(p2.sg_vis, p2.sg_tf[setdiff(rownames(p2.sg_tf), rownames(p2.sg_vis)),])
plot_timeline_ggplot(p2.sg_vis, timedata = sce_p2$Pseudotime, txtsize = 5)
ggsave('geneswitch.p2.pdf', width = 6, height = 6)


##### compare two trajectoy
sg.p1_p2 <- common_genes(p1.sg_vis, p2.sg_vis, r2cutoff = 0.1, path1name = 'branch3', path2name = 'branch1')
common_genes_plot(sg.p1_p2, timedata = sce_p1$Pseudotime)

sg_disgs.p1_p2 <- distinct_genes(p1.sg_vis, p2.sg_vis, r2cutoff = 0.1, path1name = 'branch3', path2name = 'branch1',
                                 path1time = sce_p1$Pseudotime, path2time = sce_p2$Pseudotime)
plot_timeline_ggplot(sg_disgs.p1_p2, timedata = sce_p1$Pseudotime, color_by = 'Paths', iffulltml = F, txtsize = 5)
ggsave('geneswitch.p1_p2.pdf', width = 6, height = 6)


