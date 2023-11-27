#### read scanpy data####
# @param confg

source('/mnt/snRNA_seq/snrna/code/required_func.R')
source('/mnt/snRNA_seq/snrna/code/config.R')
source('/mnt/snRNA_seq/snrna/code/paper_load_sc_data.R')

file_path <- '/mnt/snRNA_seq/snrna/result'
sample_name <- c('ctm15', 'ctm3', 'ctm9', 'ctp14', 'ctp4')

filepath = lapply(sample_name, function(sample) {
  paste0('/mnt/data/filtered/',sample, '_paper.h5ad')
})

paper_list <- load_h5(filepath)
names(paper_list) <- sample_name

paper.combine.norm_2k <- integrate_norm(paper_list)
# paper.combine.norm_5k <- integrate_norm(paper_list, nfeatures = 5000)

saveRDS(paper.combine.norm_2k, file = '/mnt/snRNA_seq/rds/22-1-5/paper.combine.norm_2k.rds')
# saveRDS(paper.combine.norm_5k, file = '/mnt/snRNA_seq/rds/22-1-5/paper.combine.norm_5k.rds')

for (i in c('umap', 'tsne')){
  write.table(Embeddings(paper.combine.norm_2k, reduction = i), file.path(file_path, paste0('paper.combine_', i,'.txt')), row.names = colnames(paper.combine.norm_2k))
}

paper.combine.norm_2k$orig.ident <- factor(paper.combine.norm_2k$orig.ident, levels = c('ctp4', 'ctp14', 'ctm3', 'ctm9', 'ctm15', 'adm9', 'adm15'))
plot_dim(paper.combine, group.by = 'seurat_clusters', pt.size = 0.01)
ggsave('paper_cluster_tsne.pdf', width = 10, height = 10, units = 'cm', path = '/mnt/snRNA_seq/snrna/result')

plot_dim(paper.combine, group.by = 'sd_type', pt.size = 0.01)
ggsave('paper_maintype_tsne.pdf', width = 10, height = 10, units = 'cm', path = '/mnt/snRNA_seq/snrna/result')

paper.combine.ct <- integrate_norm(paper_list)

