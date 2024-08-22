library(Seurat)
library(Matrix)

so <- readRDS('data/integrated_preped_so05.rds')

# so[['RNA']] <- NULL
# DefaultAssay(so) <- 'RNA'
# SaveH5Seurat(so, filename = "data/integrated_so08.h5Seurat", overwrite = T)
# Convert("data/integrated_so08.h5Seurat", dest = "h5ad")

so$barcode <- colnames(so)
so$UMAP_1 <- so@reductions$umap_rpca@cell.embeddings[,1]
so$UMAP_2 <- so@reductions$umap_rpca@cell.embeddings[,2]
write.csv(so@meta.data, file='scvelo/data/metadata.csv', quote=F, row.names=F)

message('Saved mata data')

write.csv(so@reductions$pca@cell.embeddings, file='scvelo/data/pca.csv', quote=F, row.names=F)
message('Saved PCA')

write.csv(so@reductions$umap_rpca@cell.embeddings, file='scvelo/data/umap.csv', quote=F, row.names=F)
message('Saved UMAP')

norm_counts_matrix <- Matrix(GetAssayData(so, assay='SCT', layer='scale.data'), sparse = T)
writeMM(norm_counts_matrix, file='scvelo/data/norm_counts.mtx')

message('Saved pearson residuals')

write.table(
  data.frame('gene'=rownames(norm_counts_matrix)),file='scvelo/data/features.csv',
  quote=F,row.names=F,col.names=F
)

message('Saved gene names')

write.table(
  data.frame('barcodes'=colnames(norm_counts_matrix)),file='scvelo/data/barcodes.csv',
  quote=F,row.names=F,col.names=F
)

message('Saved barcodes')

expression_matrix <- ReadMtx(
  mtx = "scvelo/data/norm_counts.mtx", features = "scvelo/data/features.csv", feature.column = 1,
  cells = "scvelo/data/barcodes.csv"
)

message('Successfully loaded expression_matrix')

