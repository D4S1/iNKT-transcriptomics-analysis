library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(SingleR)
library(ggmosaic)
library(pals)
library(grid)

# OPTIONS

so <- readRDS('data/integrated_preped_so05.rds')
method_name <- 'integrated05'
dir_path <- 'final_graphics/'
file_extn <- 'pdf' # or 'png'

if (!dir.exists(dir_path)){
  dir.create(dir_path)
}

mysave <- function(extension, filename, plot, height, width, grid=F){
  if(extension=='pdf'){
    save_func <- pdf
  }else if(extension=='png'){
    save_func <-png
  }else{
    message('Wrong file extension. Allowed extensions: pdf, png')
  }
  save_func(filename, height=height, width=width)
  if (grid){
    grid.draw(plot)
  }else{
  print(plot)
  }
  dev.off()
}


# STATS PLOTS

cp <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33","#A65628", "#F781BF", "#999999", "#F643D8", "#F6BB43", "#43F661", "#15EACA")

stat_cl <- DimPlot(object = so, reduction = "umap_rpca") + ggtitle(paste0('Clustering in ', ' data')) + scale_colour_manual(values=cp)
 
stat_tr <- DimPlot(object = so, reduction = "umap_rpca", split.by='treatment', ncol=2) + ggtitle('A') + scale_colour_manual(values=cp)
stat_dn <- DimPlot(object = so, reduction = "umap_rpca", split.by = 'donor', ncol=2) + ggtitle('B') + scale_colour_manual(values=cp)
stat_dn_tr <- DimPlot(object = so, reduction = "umap_rpca", split.by = 'donor_treatment', ncol=4) + ggtitle('C') + scale_colour_manual(values=cp)
stat_p <- DimPlot(object = so, reduction = "umap_rpca", split.by='pool', ncol=2) + ggtitle('D') + scale_colour_manual(values=cp)

stat_plot <- grid.arrange(stat_tr, stat_dn, stat_dn_tr, stat_p, ncol = 2)

mysave(file_extn, paste(dir_path, method_name, "_stat.", file_extn, sep="") , stat_plot, height = 10, width = 16, grid=T)
mysave(file_extn, paste(dir_path, method_name, "_clustering.", file_extn, sep="") , stat_cl, height = 8, width = 12)


# PHASE PLOTS
cp_phase <- c('#377EB8', "#4DAF4A", "#F6BB43")
  
phase_phase <- DimPlot(object = so, reduction = "umap_rpca", split.by='Phase') + scale_colour_manual(values=cp)
phase_tr <- DimPlot(object = so, reduction = "umap_rpca", group.by='Phase', split.by='treatment', ncol=4) + ggtitle('A') + scale_colour_manual(values=cp_phase)
phase_dn_tr <- DimPlot(object = so, reduction = "umap_rpca", group.by='Phase', split.by='donor_treatment', ncol=4) + ggtitle('B') + scale_colour_manual(values=cp_phase)
phase_p <- DimPlot(object = so, reduction = "umap_rpca", group.by='Phase', split.by='pool', ncol=2) + ggtitle('C') + scale_colour_manual(values=cp_phase)
phase_dn <- DimPlot(object = so, reduction = "umap_rpca", group.by='Phase', split.by='donor', ncol=2) + ggtitle('D') + scale_colour_manual(values=cp_phase)

phase_plot <- grid.arrange(phase_tr, phase_dn, phase_dn_tr, phase_p, ncol = 2)
 
mysave(file_extn, paste0(dir_path, method_name, "_phase_stat.", file_extn) , phase_plot, height = 10, width = 16, grid=T)
mysave(file_extn, paste0(dir_path, method_name, "_phase.", file_extn) , phase_phase, height = 8, width = 12)

# NK/T ANNOTATION
ann_cp <- c("#008080", "#454B1B", "#33a02c", "#43F661", "#15EACA", "#565656", "#6A33C2", "#F643D8", "#ff0000", "#778B00", "#ff7f00", "#0000FF", "#A52A3C", "#da7964", "#09005e", "#FF83FA", "#EEE685", "#9c65a1", "#36648B")
ann_plot_main <- DimPlot(object = so, reduction = "umap_rpca", group.by='NK_T_main', alpha = 0.6) + ggtitle('NK and T annotation') + scale_colour_manual(values=cp[c(3,2,5,11,7)])
ann_plot_fine <- DimPlot(object = so, reduction = "umap_rpca", group.by='NK_T_fine', alpha=0.5) + ggtitle('NK and T annotation') + scale_colour_manual(values=ann_cp)

mysave(file_extn, paste0(dir_path,method_name, "_nk_t_main_ann.", file_extn) , ann_plot_main, height = 8, width = 12)
mysave(file_extn, paste0(dir_path,method_name, "_nk_t_fine_ann.", file_extn) , ann_plot_fine, height = 8, width = 12)


# CLUSTER COMPOSITION PLOTS
tmp <- so@meta.data[c('donor_treatment', 'clusters_rpca', 'Phase', 'donor', 'treatment', 'NK_T_main', 'NK_T_fine', 'sctype_cell_ann')]
tmp$clusters_rpca <- factor(tmp$clusters_rpca, levels = 0:12)
tmp$sctype_cell_ann <- factor(tmp$sctype_cell_ann)

tmp_frac <- tmp %>%
  group_by(clusters_rpca, sctype_cell_ann) %>% 
  tally() %>% # or mutate(n = n())
  group_by(clusters_rpca) %>% 
  mutate(frac = n / sum(n)) %>% 
  ungroup()

tmp_frac_colored <- tmp_frac %>% 
  group_by(sctype_cell_ann) %>% 
  summarize(frac = max(frac, na.rm = TRUE)) %>% 
  arrange(desc(frac))

col1 <- c("#fff78a", "#F9D949", "#508d4e", "#3aa6b9", '#C8ACD6', "#921A40", "#FDE767", "#96c9f4", '#ecb390', "#0f67b1", "#d24545", "#FF8551", "#7469b6",'#ad88c6', '#c75b7a')
col0 <- "#C7C8CC"
tmp_frac_colored$color <- col0
tmp_frac_colored$color[1:length(col1)] <- col1

cols <- tmp_frac_colored$color
names(cols) <- tmp_frac_colored$sctype_cell_ann

dn_order <- c('donor1', 'donor2', 'donor3', 'donor4')
tmp$donor <- factor(tmp$donor, levels = rev(dn_order))
dn_comp <- ggplot(data = tmp) + geom_mosaic(aes(x = product(clusters_rpca), fill=donor)) +
  ggtitle('A') +
  scale_fill_manual(values=brewer.paired(8)[c(2,4,6,8)], limits=dn_order)

tr_comp <- ggplot(data = tmp) + geom_mosaic(aes(x = product(clusters_rpca), fill=treatment)) +
  ggtitle('B') +
  scale_fill_manual(values=c('#9784bc', '#522f96'), limits=c('IL2', 'IL15'))

dntr_order <- c("IL2 donor1", "IL15 donor1", "IL2 donor2", "IL15 donor2", "IL2 donor3", "IL15 donor3", "IL2 donor4","IL15 donor4")
tmp$donor_treatment <- factor(tmp$donor_treatment, levels = rev(dntr_order))
dn_tr_comp <- ggplot(data = tmp) + geom_mosaic(aes(x = product(clusters_rpca), fill=donor_treatment)) +
  ggtitle('C') +
  scale_fill_manual(values=brewer.paired(n=8), limits=dntr_order)

phase_comp <- ggplot(data = tmp) + geom_mosaic(aes(x = product(clusters_rpca), fill=Phase)) +
  ggtitle('D') +
  scale_fill_manual(values=cp_phase)

singler_m_comp <- ggplot(data = tmp) + geom_mosaic(aes(x = product(clusters_rpca), fill=NK_T_main)) +
  ggtitle('SingleR annotation main labels') +
  scale_fill_manual(values=cp[c(3,2,5,11,7)])

singler_f_comp <- ggplot(data = tmp) + geom_mosaic(aes(x = product(clusters_rpca), fill=NK_T_fine)) +
  ggtitle('SingleR annotation fine labels') +
  scale_fill_manual(values=ann_cp) +
  theme(axis.text.y = element_blank())

tmp$sctype_cell_ann_rev <- factor(tmp$sctype_cell_ann, levels = rev(levels(tmp$sctype_cell_ann)))
sctype_f_comp <- ggplot(data = tmp) + geom_mosaic(aes(x = product(clusters_rpca), fill=sctype_cell_ann_rev)) +
  scale_fill_manual(NULL, values=cols, limits = levels(tmp$sctype_cell_ann)) +
  ggtitle('ScType cell type annotation composition') +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x="Cluster", y='Cell type') + 
  theme(axis.text.y = element_blank())

sctype_m_comp <- DimPlot(so, reduction = "umap_rpca", label = TRUE, repel = TRUE, group.by = 'sctype_cluster_ann') +ggtitle('ScType annotation per cluster')


comp_plot <- grid.arrange(dn_comp, tr_comp, dn_tr_comp, phase_comp, ncol = 2)
mysave(file_extn, paste0(dir_path,method_name, '_comp.', file_extn), comp_plot, height = 10, width = 16, grid=T)
mysave(file_extn, paste0(dir_path,method_name, '_singleR_m_comp.', file_extn), singler_m_comp, height = 8, width = 12)
mysave(file_extn, paste0(dir_path,method_name, '_singleR_f_comp.', file_extn), singler_f_comp, height = 8, width = 12)
mysave(file_extn, paste0(dir_path,method_name, '_sctype_m_comp.', file_extn), sctype_m_comp, height = 8, width = 12)
mysave(file_extn, paste0(dir_path,method_name, '_sctype_f_comp.', file_extn), sctype_f_comp, height = 8, width = 12)


# EXPRESSION PLOTS

# Literature markers
# ZBTB16->PLZF  TBX21->T-bet  RORC->ROR gamma t
literature_markers <- c('ZBTB16', 'TBX21', 'GATA3', 'RORC', 'IL4', 'IFNG')

exp_plot <- FeaturePlot(so, features=literature_markers, ncol=3)
vln_plot <- VlnPlot(so, features = literature_markers, split.by = 'clusters_rpca', ncol=3)

mysave(file_extn, paste0(dir_path,method_name, '_lit_markers_vln.', file_extn), vln_plot, height = 10, width = 15)
mysave(file_extn, paste0(dir_path,method_name, '_lit_markers_exp.', file_extn), exp_plot, height = 10, width = 15)



# Antibodies
adt_labels <- paste0("adt_", c("CD122", "CD62L","CD314", "CD223", "CD279", "CD150", "CD319", "TCR"))

adt_exp_plot <- FeaturePlot(so, features = adt_labels, max.cutoff="q95", ncol=4) & scale_color_viridis_c()
adt_vln_plot <- VlnPlot(so, features = adt_labels, split.by = 'clusters_rpca', ncol=4)

mysave(file_extn, paste0(dir_path,method_name, '_adt_vln.', file_extn), adt_vln_plot, height = 10, width = 20)
mysave(file_extn, paste0(dir_path,method_name, '_adt_exp.', file_extn), adt_exp_plot, height = 10, width = 20)



# Interleukin

il_exp_plot <- FeaturePlot(so, features=c('IL2', 'IL15'))
il_vln_plot <- VlnPlot(so, features = c('IL2', 'IL15'), split.by = 'clusters_rpca')

mysave(file_extn, paste0(dir_path,method_name, '_IL_vln.', file_extn), il_vln_plot, height = 8, width = 16)
mysave(file_extn, paste0(dir_path,method_name, '_IL_exp.', file_extn), il_exp_plot, height = 8, width = 16)


# CLUSTER MARKERS
cluster_markers <- read.csv('data/integrated_markers05_min075.csv')

genes <- cluster_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_max(avg_log2FC, n = 3) %>%
  ungroup()

unique_genes <- genes$gene
# write.csv(genes, 'data/i05_markers_min075_top20.csv')

for (i in seq(1, 39, 3)){
  cluster = (i-1)/3
  exp_plot <- FeaturePlot(so, features=unique_genes[i:(i+2)], ncol=3)
  vln_plot <- VlnPlot(so, features=unique_genes[i:(i+2)], split.by = 'seurat_clusters', ncol=3) 
  mysave(file_extn, paste0(dir_path,method_name, '_markers_vln_c', cluster, '.', file_extn), vln_plot, height = 5, width = 15)
  mysave(file_extn, paste0(dir_path,method_name, '_markers_exp_c', cluster, '.', file_extn), exp_plot, height = 5, width = 15)
}



nkg2d_exp <- FeaturePlot(so, features = 'KLRK1', alpha=0.6) + ggtitle('NKG2D expression')
nkg2d_vln <- VlnPlot(so, features = 'KLRK1', pt.size=0) + ggtitle('NKG2D expression density in clusters')

mysave(file_extn, paste0(dir_path,method_name, '_NKG2D_exp.', file_extn), nkg2d_exp, height = 8, width = 12)
mysave(file_extn, paste0(dir_path,method_name, '_NKG2D_vln.', file_extn), nkg2d_vln, height = 8, width = 12)
