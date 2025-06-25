# NORMALIZATION
library(Seurat)
library(dplyr)
library(tidyr)
library(SeuratWrappers)


so_list <- readRDS('neo/so_list_afqc.rds')

so_list <- mapply(function(so, pref){
  so <- RenameCells(so, add.cell.id=pref)
}, so_list, c('P1', 'P1', 'P2', 'P2','P3', 'P3', 'P4', 'P4'), SIMPLIFY = FALSE)


merged_so <- merge(so_list[[1]], y = so_list[-1], project = "iNKT")

merged_so <- NormalizeData(merged_so, assay = "ADT", normalization.method = "CLR")
merged_so <- SCTransform(merged_so, assay = 'RNA', verbose = F, variable.features.n=3000)
merged_so <- RunPCA(merged_so, npcs = 30, verbose = F)

# before integration do this:
# options(future.globals.maxSize = 8000 * 1024^2)
# otherwise you will get this error:
# Error in getGlobalsAndPackages(expr, envir = envir, globals = globals) : 
# The total size of the 10 globals exported for future expression (‘FUN()’) is 2.78 GiB.. 
# This exceeds the maximum allowed size of 500.00 MiB (option 'future.globals.maxSize').
# The three largest globals are ‘object.list’ (2.78 GiB of class ‘list’), ‘NNHelper’
# (93.29 KiB of class ‘function’) and ‘FUN’ (21.25 KiB of class ‘function’)

saveRDS(merged_so,"neo/merged_so.rds")


# HARMONY INTEGRATION
integrated <- IntegrateLayers(
  object = merged_so,
  method = HarmonyIntegration,
  new.reduction = "harmony",
  normalization.method = "SCT",
  verbose = F
)

integrated <- FindNeighbors(integrated,
                            dims = 1:30,
                            reduction = "harmony",
                            graph.name = "harmony_graph")

for (i in 2:15){
  integrated <- FindClusters(integrated,
                             resolution = i/10,
                             graph.name = "harmony_graph",
                             cluster.name = paste0("harmony_clusters_tr_0", i),
                             algorithm = 4
                             )
}


# SCVI INTEGRATION - look at scvi integration python scripts
# seurat object
integrated <- readRDS("neo/merged_so.rds")

# scvi
latent <- read.csv("neo/latent_embedding_nlayers12_nlatent10.csv", row.names = 1)
latent <- latent[colnames(integrated), ]
colnames(latent) <- paste0("scvi_", seq_len(ncol(latent)))

scvi_reduction <- CreateDimReducObject(
  embeddings = as.matrix(latent),
  key = "scvi_",
  assay = DefaultAssay(integrated)
)
integrated[["scvi"]] <- scvi_reduction

# Umap
umap <- read.csv('neo/scvi_umap_coords(2).csv', row.names = 1)
umap_reduction <- CreateDimReducObject(
  embeddings = as.matrix(umap),
  key = "umap_",
  assay = DefaultAssay(integrated)
)
integrated[["umap"]] <- umap_reduction


# Clustering
clusters <- read.csv("neo/clusters(2).csv", row.names = 1)
integrated@meta.data['clusters'] <- clusters


# ==============================================================================
old_so <-  readRDS("new_data/integrated_harmony_so05.rds")
embeddings <- Embeddings(old_so, reduction = "integrated_harmony")[, 1:30]
cluster_labels <- as.integer(as.factor(old_so@meta.data[["clusters"]]))

# Silhouette score (on Euclidean distance matrix)
sil <- silhouette(cluster_labels, dist(embeddings))
avg_sil <- summary(sil)$avg.width

# Calinski-Harabasz and Davies-Bouldin indices
crit <- intCriteria(traj = embeddings, part = cluster_labels,
                    crit = c("Calinski_Harabasz", "Davies_Bouldin"))

# avg_sil, crit$calinski_harabasz, crit$davies_bouldin
# 0.085, 2657.654, 2.084 

# ==============================================================================
embeddings <- Embeddings(integrated, reduction = "scvi")[,1:10]
cluster_labels <- as.integer(as.factor(integrated@meta.data[['clusters']]))

sil <- silhouette(cluster_labels, dist(embeddings))
avg_sil <- summary(sil)$avg.width

# Calinski-Harabasz and Davies-Bouldin indices
crit <- intCriteria(traj = embeddings, part = cluster_labels,
                    crit = c("Calinski_Harabasz", "Davies_Bouldin"))

avg_sil
crit$calinski_harabasz
crit$davies_bouldin

# ==============================================================================
  
library(cluster)
library(clusterCrit)

evaluate_clustering <- function(so,
                                embedding_reduction = "harmony",
                                embedding_dims = 1:30,
                                i_range = 2:15,
                                cluster_prefix = "harmony_clusters_tr_0") {
  
  # Get the embedding matrix
  embeddings <- Embeddings(so, reduction = embedding_reduction)[, embedding_dims]
  
  # Placeholder for results
  results <- data.frame(
    resolution = numeric(),
    silhouette = numeric(),
    calinski_harabasz = numeric(),
    davies_bouldin = numeric(),
    n_clusters = integer()
  )
  
  for (i in i_range) {
    cluster_col <- paste0(cluster_prefix, i)
    
    # Extract cluster labels and ensure they're integers
    cluster_labels <- as.integer(as.factor(so@meta.data[[cluster_col]]))
    
    # Silhouette score (on Euclidean distance matrix)
    sil <- silhouette(cluster_labels, dist(embeddings))
    avg_sil <- summary(sil)$avg.width
    
    # Calinski-Harabasz and Davies-Bouldin indices
    crit <- intCriteria(traj = embeddings, part = cluster_labels,
                        crit = c("Calinski_Harabasz", "Davies_Bouldin"))
    
    # Store results
    results <- rbind(results, data.frame(
      n_clusters = length(unique(cluster_labels)),
      resolution = i / 10,
      silhouette = avg_sil,
      calinski_harabasz = crit$calinski_harabasz,
      # davies_bouldin = crit$davies_bouldin
    ))
  }
  
  return(results)
}

df <- evaluate_clustering(integrated)
# n_clusters resolution silhouette calinski_harabasz davies_bouldin
# 1           4        0.2 0.13743769          5161.093       1.926670
# 2           6        0.3 0.09852999          3691.320       2.228078
# 3           7        0.4 0.09949547          3426.422       2.202743
# 4           9        0.5 0.10367612          2987.194       1.952908
write.csv(df, "neo/clustering_harmony_metrics.csv", row.names = FALSE)

