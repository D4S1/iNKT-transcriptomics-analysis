library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(SingleR)
library(ggmosaic)

# OPTIONS
so_labels <- c('original', 'regressed_effect')
so_list <- readRDS('data/old_2000_norm_preped.rds')
method_name <- 'old_2000_norm_'
dir_path <- 'graphics2/'

# Hot fix, needed only with my files.
# so_list <- lapply(so_list, function(so){
#   names(so@meta.data)[names(so@meta.data) == "NK_T_old"] <- "NK_T"
#   return (so)
# })

# STATS PLOTS

mapply(function(so, name){
  p_cl <- DimPlot(object = so, reduction = "umap") + ggtitle(paste0('Clustering in ', name, ' data'))
  
  p_tr <- DimPlot(object = so, reduction = "umap", split.by='treatment', ncol=2) + ggtitle('A')
  p_dn <- DimPlot(object = so, reduction = "umap", split.by = 'donor', ncol=2) + ggtitle('B')
  p_dn_t <- DimPlot(object = so, reduction = "umap", split.by = 'donor_treatment', ncol=4) + ggtitle('C')
  p_p <- DimPlot(object = so, reduction = "umap", split.by='pool', ncol=2) + ggtitle('D')
  
  stat_plot <- grid.arrange(p_tr, p_dn, p_dn_t, p_p, ncol = 2)
  ggsave(paste(dir_path, method_name,  name, "_stat.png", sep="") , stat_plot, height = 10, width = 16)
  ggsave(paste(dir_path, method_name,  name, "_clustering.png", sep="") , p_cl, height = 8, width = 12)
}, so_list, so_labels)


# PHASE PLOTS

mapply(function(so, name){
  p_phase <- DimPlot(object = so, reduction = "umap", split.by='Phase')  + ggtitle('A')
  p_dn_tr <- DimPlot(object = so, reduction = "umap", group.by='Phase', split.by='donor_treatment', ncol=4) + ggtitle('B')
  p_p <- DimPlot(object = so, reduction = "umap", group.by='Phase', split.by='pool', ncol=2) + ggtitle('C')
  p_dn <- DimPlot(object = so, reduction = "umap", group.by='Phase', split.by='donor', ncol=2) + ggtitle('D')
  stat_plot <- grid.arrange(p_phase, p_dn, p_dn_tr, p_p, ncol = 2)
  
  ggsave(paste(dir_path, method_name, name, "_phase_stat.png", sep="") , stat_plot, height = 10, width = 16)
}, so_list, so_labels)

mapply(function(so, name){
  so <- RunPCA(so, features = c(cc.genes$s.genes, cc.genes$g2m.genes))
  cc_plot <- DimPlot(so, reduction='pca', group.by = 'Phase') + ggtitle(paste(name, 'cell cycle phase pca', sep=" "))
  ggsave(paste(dir_path,method_name, name, '_cc_pca.png', sep='' ), cc_plot) 
}, so_list, so_labels)


# NK/T ANNOTATION

mapply(function(so, name){
  ann_plot <- DimPlot(object = so, reduction = "umap", group.by='NK_T') +
    ggtitle(paste(name, 'NK and T annotation', sep=" "))
  
  ggsave(paste(dir_path,method_name, name, "_nk_t_ann.png", sep="") , ann_plot, height = 8, width = 10)
}, so_list, so_labels)


# CLUSTER COMPOSITION PLOTS

mapply(function(so, name){

  tr_comp <- ggplot(data = so@meta.data) +
    geom_mosaic(aes(x = product(seurat_clusters), fill=treatment)) +
    ggtitle('A')
  
  dn_tr_comp <- ggplot(data = so@meta.data) +
    geom_mosaic(aes(x = product(seurat_clusters), fill=donor_treatment)) +
    ggtitle('B')
  
  phase_comp <- ggplot(data = so@meta.data) +
    geom_mosaic(aes(x = product(seurat_clusters), fill=Phase)) +
    ggtitle('C')
  
  nk_comp <- ggplot(data = so@meta.data) +
    geom_mosaic(aes(x = product(seurat_clusters), fill=NK_T)) +
    ggtitle('D') +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

  comp_plot <- grid.arrange(tr_comp, dn_tr_comp, phase_comp, nk_comp, ncol = 2)
  ggsave(paste0(dir_path,method_name, name, '_comp.png'), comp_plot, height = 10, width = 14)
}, so_list, so_labels)


# EXPRESSION PLOTS

# Literature markers
# ZBTB16->PLZF    TBX21->T-bet    RORC->ROR gamma t
literature_markers <- c('ZBTB16', 'TBX21', 'GATA3', 'RORC', 'RORA', 'IL4', 'IFNG')

mapply(function(so, name){
  exp_plot <- FeaturePlot(so, features=literature_markers, ncol=4)
  vln_plot <- VlnPlot(so, features = literature_markers, split.by = 'seurat_clusters', ncol=4)

  ggsave(paste(dir_path,method_name, name, '_lit_markers_vln.png', sep=""), vln_plot, height = 10, width = 15)
  ggsave(paste(dir_path,method_name, name, '_lit_markers_exp.png', sep=""), exp_plot, height = 8, width = 16)
}, so_list, so_labels)


# Antibodies
adt_labels <- paste0("adt_", c("CD122", "CD62L","CD314", "CD223", "CD279", "CD150", "CD319", "TCR"))

mapply(function(so, name){
  exp_plot <- FeaturePlot(so, features = adt_labels, max.cutoff="q95", ncol=4) & scale_color_viridis_c()
  vln_plot <- VlnPlot(so, features = adt_labels, split.by = 'seurat_clusters', ncol=4)
  
  ggsave(paste(dir_path,method_name, name, '_adt_vln.png', sep=""), vln_plot, height = 10, width = 15)
  ggsave(paste(dir_path,method_name, name, '_adt_exp.png', sep=""), exp_plot, height = 12, width = 18)
}, so_list, so_labels)


# Interleukin
mapply(function(so, name){
  exp_plot <- FeaturePlot(so, features=c('IL2', 'IL15'))
  vln_plot <- VlnPlot(so, features = c('IL2', 'IL15'), split.by = 'seurat_clusters')
  
  ggsave(paste(dir_path,method_name, name, '_IL_vln.png', sep=""), vln_plot, height = 10, width = 15)
  ggsave(paste(dir_path,method_name, name, '_IL_exp.png', sep=""), exp_plot, height = 8, width = 16)
}, so_list, so_labels)

# CLUSTER MARKERS

cluster_markers <- readRDS('data/sct_norm_markers.rds')

marker_list <- lapply(cluster_markers, function(m){
  genes <- m %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_max(avg_log2FC, n = 3) %>%
    ungroup()
  
  gene_list <- c(genes$gene)
  unique_gene_list <- unique(gene_list)

  return(gene_list)
})

mapply(function(so, name, markers){
  exp_plot <- FeaturePlot(so, features=markers, ncol=4)
  vln_plot <- VlnPlot(so, features = markers, split.by = 'seurat_clusters', ncol=4)
  
  ggsave(paste(dir_path,method_name, name, '_markers_vln.png', sep=""), vln_plot, height = 26, width = 20)
  ggsave(paste(dir_path,method_name, name, '_markers_exp.png', sep=""), exp_plot, height = 26, width = 20)
}, so_list, so_labels, marker_list)

saveRDS(cluster_markers, 'data/old_norm_markers.rds')