library(Seurat)
library(dplyr)
library(tidyr)
library(SingleR)
library(presto)

# OPTIONS

file <- 'data/merged_old_2000_norm.rds'
method_name <- 'old_2000_norm_' # sct_norm
assay <-'RNA' # SCT
regressed_file <- NULL # if path provided, script read the file instead of calculating new object
resolution <- 0.3


# Original and regressed effect objects creation

original_so <- readRDS(file)
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
original_so <- CellCycleScoring(original_so, s.features = s_genes, g2m.features = g2m_genes)

if(is.null(regressed_file)){
  regressed_so <- original_so
  regressed_so <- ScaleData(regressed_so, vars.to.regress = c('S.Score', 'G2M.Score'), features = rownames(regressed_so))
  saveRDS(regressed_so, paste0('data/', method_name, 'regressed.rds'))
}else{
  regressed_so <- readRDS(regressed_file)
}

so_labels <- c('original', 'regressed_effect')
so_list <- c(original_so, regressed_so)

message('Generation of objects finished')

# Donor treatment field in meta data

so_list <- lapply(so_list, function(so){
  so[['donor_treatment']] <- paste(so@meta.data[['treatment']], so@meta.data[['donor']])
  return (so)
})

# PCA

so_list <- lapply(so_list, RunPCA)

message('PCA finished')

# Clustering

so_list <- lapply(so_list, function(so){
  so <- FindNeighbors(object = so, dims = 1:10)
  so <- FindClusters(object = so, resolution=resolution, verbose=FALSE) # Louvain algorithm
  # renumber the clusters starting from 1
  so$seurat_clusters <- as.factor(as.numeric(so$seurat_clusters))
  Idents(so) <- "seurat_clusters"
  so <- RunUMAP(object = so, dims = 1:10, verbose=FALSE)
  return (so)
})

message('Clustering finished')

# NK/T annotation

hpca_ref_old <- HumanPrimaryCellAtlasData()
# hpca_ref_new <- celldex::HumanPrimaryCellAtlasData()

so_list <- lapply(so_list, function(so){
  singleR_old <- SingleR(test = so[[assay]]$scale.data,
                             ref = hpca_ref_old,
                             labels = hpca_ref_old$label.main,
                             assay.type.test=1
                             )
  so[["NK_T"]]<-singleR_old$pruned.labels
  return (so)
  })

message('NK/T annotation finished')

saveRDS(so_list, paste0('data/', method_name, 'preped.rds'))

# Markers
if(assay == 'SCT'){
  so_list<-lapply(so_list, PrepSCTFindMarkers)
}

markers_df_list <-lapply(so_list, FindAllMarkers)
saveRDS(markers_df_list, paste0('data/', method_name, 'markers.rds'))

message('Searching for markers finished')
