data0<-data
data0<- NormalizeData(data0, normalization.method = "LogNormalize", scale.factor = 10000)
data0<- FindVariableFeatures(data0, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(data0)
data0 <- ScaleData(data0, features = all.genes)
data0 <- RunPCA(data0)

data0 <- FindNeighbors(data0, dims = 1:10)
data0 <- FindClusters(data0, resolution = 0.2)
data0 <- RunUMAP(data0, dims = 1:10)
DimPlot(data0)
markers<-FindAllMarkers(data0)

