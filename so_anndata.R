library(Seurat)
library(reticulate)
library(sceasy)
library(Matrix)

# SETUP PYTHON
use_virtualenv("/home/hoshi/ME/Coding/inkt/env", T)
sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)
anndata <- reticulate::import("anndata", convert = FALSE)
scipy_sparse <- import("scipy.sparse", convert = FALSE)

# SELECT HVG
hvg <- merged_so$SCT@var.features
merged_so <- merged_so[hvg]
saveRDS(merged_so, "../data/hvg_so.rds")

# LOAD OBJECT
merged_so <- readRDS("data/integrated_hvg.rds")

# get scaled_data
data <- as.matrix(GetAssayData(merged_so, assay = 'SCT', layer = 'data'))

# save to CSV
write.csv(data, file = "scaled_data.csv")


# CONVERT TO ANNDATA
counts <- t(Matrix(GetAssayData(merged_so, assay='SCT', layer='counts'), sparse = T))
counts <- scipy_sparse$csr_matrix(r_to_py(counts))
obs <- merged_so@meta.data[c('pool', 'donor', 'treatment', 'donor_treatment')]
var <- data.frame(gene = rownames(merged_so))
rownames(var) <- rownames(merged_so)

adata <- anndata$AnnData(
  X = counts,
  obs = obs,
  var = var,
)
adata$write("data/hvg.h5ad")

# LOAD EMBEDDINGS
# latent <- read.csv("data/latent_embedding_treatment.csv", row.names = 1, header = 1)
latent <- as.matrix(latent[colnames(merged_so), ])

scvi_reduction <- CreateDimReducObject(
  embeddings = latent,
  key = "scvi_",
  assay = DefaultAssay(merged_so)
)

# Add to Seurat object
merged_so[["scvi"]] <- scvi_reduction
merged_so <- RunUMAP(merged_so, reduction = "scvi", dims = 1:30)
DimPlot(merged_so, reduction = "umap", group.by = "treatment")

merged_so <- FindNeighbors(merged_so, reduction = "scvi", dims = 1:30)
merged_so <- FindClusters(merged_so, resolution = 0.5, reduction = "scvi", cluster.name = "clusters")
DimPlot(merged_so, reduction = "umap", group.by = "seurat_clusters")

cluster_name <- "clusters"

#################################################################
#                   EVALUATION METRICS
#################################################################

library(cluster)
library(clusterCrit)

# Harmony
harmony_embeddings <- Embeddings(merged_so, reduction = "harmony")

harmony_clusters <- merged_so@meta.data$...
harmony_clusters <- as.integer(as.factor(harmony_clusters))

sil_harmony <- silhouette(clusters, dist(harmony_embeddings))
avg_sil_harmony <- summary(sil_harmony)$avg.width
crit_harmony <- intCriteria(traj = harmony_embeddings, part = harmony_clusters, crit = c("Calinski_Harabasz", "Davies_Bouldin"))

# scVI
scvi_embeddings <- Embeddings(merged_so, reduction = "scvi")

scvi_clusters <- merged_so@meta.data$seurat_clusters
scvi_clusters <- as.integer(as.factor(scvi_clusters))

sil_scvi <- silhouette(scvi_clusters, dist(scvi_embeddings))
avg_sil_scvi <- summary(sil_scvi)$avg.width
crit_scvi <- intCriteria(traj = scvi_embeddings, part = scvi_clusters, crit = c("Calinski_Harabasz", "Davies_Bouldin"))

cat("Harmony: Sil =", avg_sil_harmony, ", CH =", crit_harmony$calinski_harabasz, ", DB =", crit_harmony$davies_bouldin, "\n")
cat("scVI: Sil =", avg_sil_scvi, ", CH =", crit_scvi$calinski_harabasz, ", DB =", crit_scvi$davies_bouldin, "\n")



#################################################################
#                       VISUALIZATION
#################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(ggmosaic)
library(pals)
library(grid)

general_palette <- c(
  "#E41A1C", "#008080", "#FF7F00", "#6A33C2", "#FFFF33",
  "#09005e", "#4DAF4A", "#F643D8", "#A65628", "#0000FF",
  "#F781BF", "#33a02c", "#984EA3", "#EEE685", "#377EB8",
  "#ff0000", "#565656", "#FF83FA", "#778B00", "#F6BB43",
  "#43F661", "#A52A3C", "#15EACA", "#da7964", "#999999",
  "#36648B", "#FB8072", "#9c65a1", "#454B1B", "#FF7F00"
)

phase_palette <- c('#377EB8', "#4DAF4A", "#F6BB43")
treatment_palette <- c('#9784bc', '#522f96')
donor_palette <- brewer.paired(8)[c(2,4,6,8)]
donor_treatment_palette <- brewer.paired(8)

subset_metadata <- c('donor_treatment', 'clusters', 'donor', 'treatment')
tmp <- merged_so@meta.data[subset_metadata]
tmp$clusters <- factor(tmp$seurat_clusters, levels = 0:(length(unique(tmp$clusters))-1))


dn_order <- c('donor1', 'donor2', 'donor3', 'donor4')
tmp$donor <- factor(tmp$donor, levels = rev(dn_order))
dn_comp <- ggplot(data = tmp) + 
  geom_mosaic(aes(x = product(cluster_name), fill = donor)) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 24),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 24)
  ) +
  scale_fill_manual(values = donor_palette, limits = dn_order)

tr_comp <- ggplot(data = tmp) + 
  geom_mosaic(aes(x = product(seurat_clusters), fill = treatment)) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 24),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 24)
  ) +
  scale_fill_manual(values = treatment_palette, limits = c('IL2', 'IL15'))

dntr_order <- c("IL2 donor1", "IL15 donor1", "IL2 donor2", "IL15 donor2", 
                "IL2 donor3", "IL15 donor3", "IL2 donor4", "IL15 donor4")

tmp$donor_treatment <- factor(tmp$donor_treatment, levels = rev(dntr_order))

dn_tr_comp <- ggplot(data = tmp) + 
  geom_mosaic(aes(x = product(seurat_clusters), fill = donor_treatment)) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 24),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 24)
  ) +
  scale_fill_manual(values = donor_treatment_palette, limits = dntr_order)


comp_plot <- grid.arrange(dn_comp, tr_comp, dn_tr_comp, ncol = 2)


