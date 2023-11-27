#### read scanpy data####

# @param config

seu_to_loom <- function(seu, filename){
  require(reticulate)
  sc <- import('scanpy')
  lp <- import('loompy')
  np <- import('numpy')
  row_attrs <- data.frame(Gene = np$array(rownames(seu)))
  col_attrs <- data.frame(CellID = np$array(colnames(seu)))
  lp$create(paste0('/mnt/data/loom/', filename, '.loom'), as.matrix(seu@assays$RNA@counts), row_attrs, col_attrs)
}

load_h5 <- function(filepath){
  require(Seurat)
  require(reticulate)
  require(future)
  require(future.apply)
  print('load scanpy data')
  seurat_list = list()
  plan('multiprocess', workers = 10)
  seurat_list <- future_lapply(filepath, function(file){
    sc <- import('scanpy')
    print(paste0(file, ' is doning'))
    sample <- sc$read_h5ad(file)
    count <- t(as.matrix(sample$X))
    rownames(count) <- unlist(rownames(sample$var))
    colnames(count) <- unlist(rownames(sample$obs))
    seurat <- CreateSeuratObject(counts = count, project = 'MS/VDB', min.cells = 3, min.features = 200)
    seurat <- AddMetaData(seurat, data.frame(sample$obs$doublet_info, row.names = colnames(seurat)), 'doublet_info')
    
    seurat <- AddMetaData(seurat, sample$obs$iqr_outlier, 'iqr_outliers')
    seurat <- AddMetaData(seurat, sample$obs$std_outlier, 'std_outliers')
    seurat <- AddMetaData(seurat, data.frame(sample$obs$pct_counts_hb, row.names = colnames(seurat)), 'pct_counts_hb')
    seurat <- AddMetaData(seurat, data.frame(sample$obs$pct_counts_ribo, row.names = colnames(seurat)), 'pct_counts_ribo')
    seurat <- AddMetaData(seurat, data.frame(sample$obs$pct_counts_mt, row.names = colnames(seurat)), 'pct_counts_mt')
    
    seurat[['sex']] <- colnames(seurat) %in% WhichCells(seurat, expression = Xist > 0)
    
    filter_feature <- rownames(seurat)[-which(rownames(seurat) %in% c(config@mt.genes, config@rb.genes, config@hb.genes, config@sex.genes))]
    # filter_cell = Cells(subset(seurat, std_outliers == FALSE & doublet_info =='False'))
    seurat <- subset(seurat, features = filter_feature)
    print(paste0(file, ' finished'))
    return(seurat)
  })
  return(seurat_list)
}

integrate_norm <- function(seurat_list, nfeatures = 2000){
  require(future)
  require(future.apply)
  options(future.globals.maxSize = 20000 * 1024^2)
  plan('multiprocess', workers = 10)
  seurat_list <-future_lapply(seurat_list, function(x){
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, nfeatures = nfeatures)
  })
  #plan('sequential')
  features = SelectIntegrationFeatures(object.list = seurat_list, nfeatures = nfeatures)
  #plan('multiprocess', workers = 10)
  seu.anchor <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
  #plan('sequential')
  seu.combine <- IntegrateData(anchorset = seu.anchor, dims = 1:30)
  DefaultAssay(seu.combine) <- 'integrated'
  seu.combine <- subset(seu.combine, doublet_info == 'False')
  #plan('multiprocess', workers = 10)
  seu.combine <- ScaleData(seu.combine, vars.to.regress = 'pct_counts_mt', verbose = FALSE)
  #plan('sequential')
  seu.combine <- RunPCA(seu.combine, npcs = 30, verbose = FALSE)
  seu.combine <-RunUMAP(seu.combine, reduction = 'pca', dims = 1:30)
  seu.combine <- FindNeighbors(seu.combine, reduction = 'pca', dims = 1:30)
  #plan('multiprocess', workers = 10)
  seu.combine <- FindClusters(seu.combine, resolution = seq(0.8, 1.2, 0.1), method = 'igraph')
  #plan('sequential')
  seu.combine <- RunTSNE(seu.combine, reduction = 'pca', dims = 1 :30, seed.use = 2021)
  return(seu.combine)
}

integrate_SCT <- function(seurat_list, nfeatures = 2000){
  require(future)
  require(future.apply)
  require(glmGamPoi)
  options(future.globals.maxSize = 20000 * 1024^2)
  plan('multiprocess', workers = 10)
  seurat_list <- future_lapply(seurat_list, FUN = SCTransform, method = 'glmGamPoi')
  features <- SelectIntegrationFeatures(seurat_list, nfeatures = 3000)
  seurat_list <- PrepSCTIntegration(seurat_list, anchor.features = features)
  seurat.anchors <- FindIntegrationAnchors(seurat_list, normalization.method = 'SCT', anchor.features = features, dims = 1:npc)
  seurat.combine.sct <- IntegrateData(seurat.anchors, normalization.method = 'SCT', dims = 1:npc) 
  
  seurat.combine.sct <- RunPCA(seurat.combine.sct, npcs = 50, verbose = FALSE) 
  seurat.combine.sct <- RunUMAP(seurat.combine.sct, reduction = 'pca', dims = 1:npc)
  seu.combine.sct <- FindNeighbors(seu.combine.sct, reduction = 'pca', dims = 1:30)
  #plan('multiprocess', workers = 10)
  seu.combine.sct <- FindClusters(seu.combine.sct, resolution = seq(0.8, 1.2, 0.1), method = 'igraph')
  #plan('sequential')
  seu.combine.sct <- RunTSNE(seu.combine.sct, seed.use = 2021)
  return(seurat.combine.sct)
}



# 
# saveRDS(seurat_list, file = '/mnt/snRNA_seq/rds/22-1-5/paper.rds')
# 
# library(Hmisc)
# cc_genes <- lapply(cc.genes, function(x){
#   x <- unlist(lapply(x, function(y){
#     return(capitalize(tolower(y)))
#   }))
#   return(x)
# })
# saveRDS(cc_genes, file = '/mnt/snRNA_seq/rds/ccgenes.rds')
# 
# lp <- import('loompy')
# row_attrs <- as.data.frame(rownames(seurat@assays$SCT@counts))
# names(row_attrs) <- 'Genes'
# col_attrs <- as.data.frame(rownames(seurat@meta.data))
# names(col_attrs) <- ' CellID'
# lp$create(paste0("/mnt/data/loom/",'chat',".loom"), Matrix(seurat@assays$SCT@counts, sparse = TRUE), row_attrs, col_attrs)
