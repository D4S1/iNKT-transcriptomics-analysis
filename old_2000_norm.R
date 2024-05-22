library(Seurat)
library(dplyr)
library(tidyr)


setwd("/home/jd438446/iNKT_project/")
so_list_old_2000_norm <- readRDS('data/so_list_afqc.rds')

so_list_old_2000_norm <- mapply(function(so, pref){
  so <- NormalizeData(so, assay = "ADT", normalization.method = "CLR")
  
  # Normalization of raw counts
  so <- NormalizeData(so, assay = 'RNA', verbose = FALSE)
  so <- FindVariableFeatures(so)
  so <- ScaleData(so, assay='RNA', verbose=FALSE)
  
  # Prep for merge
  so <- RenameCells(so, add.cell.id=pref)
  return (so)
}, so_list_old_2000_norm, c('P1', 'P2', 'P3', 'P4'), SIMPLIFY = FALSE)

merged_old_2000_norm <- merge(so_list_old_2000_norm[[1]], y = so_list_old_2000_norm[-1],
                         project = "iNKT",
                         merge.data = TRUE,
)
merged_old_2000_norm<- JoinLayers(merged_old_2000_norm, assay = 'RNA')
saveRDS(merged_old_2000_norm, 'data/merged_old_2000_norm.rds')
