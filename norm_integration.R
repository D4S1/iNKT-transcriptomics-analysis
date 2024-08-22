# NORMALIZATION
library(Seurat)
library(dplyr)
library(tidyr)
library(SeuratWrappers)


so_list_sct_norm <- readRDS('data/so_list_afqc.rds')

so_list_sct_norm <- mapply(function(so, pref){
  so <- RenameCells(so, add.cell.id=pref)
}, so_list_sct_norm, c('P1', 'P2', 'P3', 'P4'), SIMPLIFY = FALSE)

merged_sct_norm <- merge(so_list_sct_norm[[1]], y = so_list_sct_norm[-1],
                         project = "iNKT",
                         )

merged_sct_norm <- NormalizeData(merged_sct_norm, assay = "ADT", normalization.method = "CLR")
merged_sct_norm <- SCTransform(merged_sct_norm, assay = 'RNA', verbose = FALSE)
merged_sct_norm <- RunPCA(merged_sct_norm, npcs = 30, verbose = F)

# before integration do this:
# options(future.globals.maxSize = 8000 * 1024^2)
# otherwise you will get this error:
# Error in getGlobalsAndPackages(expr, envir = envir, globals = globals) : 
# The total size of the 10 globals exported for future expression (‘FUN()’) is 2.78 GiB.. 
# This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize').
# The three largest globals are ‘object.list’ (2.78 GiB of class ‘list’), ‘NNHelper’
# (93.29 KiB of class ‘function’) and ‘FUN’ (21.25 KiB of class ‘function’)

options(future.globals.maxSize = 8000 * 1024^2)

merged_sct_norm<- IntegrateLayers(
  object = merged_sct_norm,
  method = RPCAIntegration,
  new.reduction = "integrated_rpca",
  normalization.method = "SCT",
  verbose = F
)

# for scvi to work you need to install: reticulate (via CRAN), scvi-tools (via conda env)
# https://github.com/pybind/pybind11/discussions/3453
# merged_sct_norm<- IntegrateLayers(
#   object = merged_sct_norm,
#   method = scVIIntegration,
#   new.reduction = "integrated_scvi",
#   normalization.method = "SCT",
#   conda_env = "/home/jd438446/miniconda3/envs/scvi-env",
#   verbose = F
# )


merged_sct_norm <- FindNeighbors(merged_sct_norm, dims = 1:30, reduction = "integrated_rpca")
merged_sct_norm <- FindClusters(merged_sct_norm, resolution = 0.5, cluster.name = "clusters_rpca")
merged_sct_norm <- RunUMAP(merged_sct_norm, dims = 1:30, reduction = "integrated_rpca", reduction.name = "umap_rpca")

saveRDS(merged_sct_norm, 'data/integrated_so05.rds')
