# NORMALIZATION
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# setwd(<project dir>)
so_list_sct_norm <- readRDS('data/so_list_afqc.rds')


so_list_sct_norm <- mapply(function(so, pref){
  # future: dsb method for antibodies normalization
  so <- NormalizeData(so, assay = "ADT", normalization.method = "CLR")
  
  # Transformed data will be available in the SCT assay,
  # which is set as the default after running sctransform
  so <- SCTransform(so, assay = 'RNA', verbose = FALSE)
  so <- RenameCells(so, add.cell.id=pref)
  return (so)
}, so_list_sct_norm, c('P1', 'P2', 'P3', 'P4'), SIMPLIFY = FALSE)

# future me: collapse option is not in Seurat yet, 
# despite being mentioned in documentation.
# Seurat v5 documentation recommend using JoinLayers 
# to merge RNA/ADT assays but this function doesn't support SCTAssay.

# update: merge work just fine without collapse, we still get
# one merged matrix, despite having list of pools

merged_sct_norm <- merge(so_list_sct_norm[[1]], y = so_list_sct_norm[-1],
                   project = "iNKT",
                   merge.data=TRUE,
                   )

VariableFeatures(merged_sct_norm$SCT) <- rownames(merged_sct_norm$SCT@scale.data)
saveRDS(merged_so, 'data/merged_sct_norm.rds')
