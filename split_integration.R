# NORMALIZATION
library(Seurat)
library(dplyr)
library(tidyr)
library(SeuratWrappers)


so_list <- readRDS('new_data/so_list_afqc.rds')

so_list <- mapply(function(so, pref){
  so <- RenameCells(so, add.cell.id=pref)
}, so_list, c('P1', 'P1', 'P2', 'P2','P3', 'P3', 'P4', 'P4'), SIMPLIFY = FALSE)


split_so_by_treatment <- function(so_list) {
  # Initialize lists for IL2 and IL15 objects
  il2 <- c()
  il15 <- c()
  
  # Iterate over the list of Seurat objects
  for (obj in so_list) {
    # Check the unique treatment value in the metadata
    treatment <- unique(obj@meta.data$treatment)
    
    # Ensure only one treatment per object
    if (length(treatment) != 1) {
      stop(paste("Seurat object", obj, "has multiple or no treatments specified."))
    }
    
    # Split the object into appropriate list based on treatment
    if (treatment == "IL2") {
      il2 <- c(il2, list(obj))
    } else if (treatment == "IL15") {
      il15 <- c(il15, list(obj))
    } else {
      stop(paste("Unknown treatment type '", treatment, "' found in object", obj, "."))
    }
  }
  
  # Return a list containing both lists
  return(list(IL2 = il2, IL15 = il15))
}
merge_treatment <- function(so_list, treatment){
  # merge all object from the same treatment into one
  merged_so <- merge(so_list[[1]], y = so_list[-1], project = treatment)
  
  # normalize the data
  merged_so <- NormalizeData(merged_so, assay = "ADT", normalization.method = "CLR")
  merged_so <- SCTransform(merged_so, assay = 'RNA', verbose = FALSE)
  merged_so <- RunPCA(merged_so, npcs = 30, verbose = F)
  
  options(future.globals.maxSize = 8000 * 1024^2)
  
  merged_so <- IntegrateLayers(
    object = merged_so,
    method = RPCAIntegration,
    new.reduction = "integrated_rpca",
    normalization.method = "SCT",
    verbose = F
  )
  return(merged_so)
}

merged_so <- merge(so_list[[1]], y = so_list[-1], project = "iNKT")

merged_so <- NormalizeData(merged_so, assay = "ADT", normalization.method = "CLR")
merged_so <- SCTransform(merged_so, assay = 'RNA', verbose = FALSE)
merged_so <- RunPCA(merged_so, npcs = 30, verbose = F)

# before integration do this:
# options(future.globals.maxSize = 8000 * 1024^2)
# otherwise you will get this error:
# Error in getGlobalsAndPackages(expr, envir = envir, globals = globals) : 
# The total size of the 10 globals exported for future expression (‘FUN()’) is 2.78 GiB.. 
# This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize').
# The three largest globals are ‘object.list’ (2.78 GiB of class ‘list’), ‘NNHelper’
# (93.29 KiB of class ‘function’) and ‘FUN’ (21.25 KiB of class ‘function’)

options(future.globals.maxSize = 8000 * 1024^2)

integrated<- IntegrateLayers(
  object = merged_so,
  method = RPCAIntegration,
  new.reduction = "integrated_rpca",
  normalization.method = "SCT",
  verbose = F
)

integrated <- FindNeighbors(integrated, dims = 1:30, reduction = "integrated_rpca")
integrated <- FindClusters(integrated, resolution = 0.5, reduction = "integrated_rpca", cluster.name = "clusters")
integrated <- RunUMAP(integrated, dims = 1:30, reduction = "integrated_rpca", reduction.name = "umap_rpca")
saveRDS(integrated, 'new_data/integrated_rpca_so05.rds')



integrated2<- IntegrateLayers(
  object = merged_so,
  method = HarmonyIntegration,
  new.reduction = "integrated_harmony",
  normalization.method = "SCT",
  verbose = F
)

integrated2 <- FindNeighbors(integrated2, dims = 1:30, reduction = "integrated_harmony")
integrated2 <- FindClusters(integrated2, resolution = 0.5, reduction = "integrated_harmony", cluster.name = "clusters")
integrated2 <- RunUMAP(integrated2, dims = 1:30, reduction = "integrated_harmony", reduction.name = "umap_harmony")

saveRDS(integrated2, 'new_data/integrated_harmony_so05.rds')
