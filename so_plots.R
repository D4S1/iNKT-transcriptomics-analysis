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

so <- readRDS('new_data/integrated_harmony_preped05.rds')
method_name <- 'harmony05'
dir_path <- 'new_figures/'
file_extn <- 'pdf' # or 'png'
reduction <- "umap_harmony"

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

# COLOR PALLETS:

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


# STATS PLOTS

stat_cl <- DimPlot(object = so, reduction = reduction, label = T) + ggtitle('Clustering in data') + theme_minimal(base_size = 20) + scale_colour_manual(values=general_palette)
 
stat_tr <- DimPlot(object = so, reduction = reduction, split.by='treatment', ncol=2) +
  scale_colour_manual(values=general_palette) +
  theme_minimal(base_size = 20)

stat_dn <- DimPlot(object = so, reduction = reduction, split.by = 'donor', ncol=2)+
  theme_minimal(base_size = 20) +
  scale_colour_manual(values=general_palette)
  
stat_dn_tr <- DimPlot(object = so, reduction = reduction, split.by = 'donor_treatment', ncol=4)+
  theme_minimal(base_size = 20) +
  scale_colour_manual(values=general_palette)

stat_p <- DimPlot(object = so, reduction = reduction, split.by='pool', ncol=2)+
  theme_minimal(base_size = 20) +
  scale_colour_manual(values=general_palette)

stat_plot <- grid.arrange(stat_tr, stat_dn, stat_dn_tr, stat_p, ncol = 2)

mysave(file_extn, paste(dir_path, method_name, "_stat.", file_extn, sep="") , stat_plot, height = 10, width = 20, grid=T)
mysave(file_extn, paste(dir_path, method_name, "_clustering.", file_extn, sep="") , stat_cl, height = 8, width = 12)


# PHASE PLOTS

phase_phase <- DimPlot(object = so, reduction = reduction, split.by='Phase') +
  theme_minimal(base_size = 20) +
  scale_colour_manual(values=general_palette)

phase_tr <- DimPlot(object = so, reduction = reduction, group.by='Phase', split.by='treatment', ncol=4) +
  ggtitle('A') +
  scale_colour_manual(values=phase_palette)

phase_dn_tr <- DimPlot(object = so, reduction = reduction, group.by='Phase', split.by='donor_treatment', ncol=4) +
  ggtitle('B') +
  scale_colour_manual(values=phase_palette)

phase_p <- DimPlot(object = so, reduction = reduction, group.by='Phase', split.by='pool', ncol=2) +
  ggtitle('C') +
  scale_colour_manual(values=phase_palette)

phase_dn <- DimPlot(object = so, reduction = reduction, group.by='Phase', split.by='donor', ncol=2) +
  ggtitle('D') +
  scale_colour_manual(values=phase_palette)

phase_plot <- grid.arrange(phase_tr, phase_dn, phase_dn_tr, phase_p, ncol = 2)
 
mysave(file_extn, paste0(dir_path, method_name, "_phase_stat.", file_extn) , phase_plot, height = 10, width = 16, grid=T)
mysave(file_extn, paste0(dir_path, method_name, "_phase.", file_extn) , phase_phase, height = 8, width = 12)


# NK/T ANNOTATION
# ann_general_palette <- c("#008080", "#454B1B", "#33a02c", "#43F661", "#15EACA", "#565656", "#6A33C2", "#F643D8", "#ff0000", "#778B00", "#ff7f00", "#0000FF", "#A52A3C", "#da7964", "#09005e", "#FF83FA", "#EEE685", "#9c65a1", "#36648B", "#FFD700")
ann_plot_main <- DimPlot(object = so, reduction = reduction, group.by='NK_T_main', alpha = 0.6) +
  ggtitle('NK and T annotation') +
  scale_colour_manual(values=general_palette)

ann_plot_fine <- DimPlot(object = so, reduction = reduction, group.by='NK_T_fine', alpha=0.5) +
  ggtitle('NK and T annotation') +
  scale_colour_manual(values=general_palette)

mysave(file_extn, paste0(dir_path,method_name, "_nk_t_main_ann.", file_extn) , ann_plot_main, height = 8, width = 12)
mysave(file_extn, paste0(dir_path,method_name, "_nk_t_fine_ann.", file_extn) , ann_plot_fine, height = 8, width = 12)


# CLUSTER COMPOSITION PLOTS
subset_metadata <- c('donor_treatment', 'clusters', 'Phase', 'donor', 'treatment', 'NK_T_main', 'NK_T_fine', 'sctype_cluster_ann', 'sctype_cell_ann')
tmp <- so@meta.data[subset_metadata]
tmp$clusters <- factor(tmp$clusters, levels = 0:(length(unique(tmp$clusters))-1))
tmp$sctype_cell_ann <- factor(tmp$sctype_cell_ann)

tmp_frac <- tmp %>%
  group_by(clusters, sctype_cell_ann) %>% 
  tally() %>% # or mutate(n = n())
  group_by(clusters) %>% 
  mutate(frac = n / sum(n)) %>% 
  ungroup()

tmp_frac_colored <- tmp_frac %>% 
  group_by(sctype_cell_ann) %>% 
  summarize(frac = max(frac, na.rm = TRUE)) %>% 
  arrange(desc(frac))

# col1 <- c("#fff78a", "#F9D949", "#508d4e", "#3aa6b9", '#C8ACD6', "#921A40", "#FDE767", "#96c9f4", '#ecb390', "#0f67b1", "#d24545", "#FF8551", "#7469b6",'#ad88c6', '#c75b7a')
# col0 <- "#C7C8CC"
tmp_frac_colored$color <- "#C7C8CC"
tmp_frac_colored$color[1:15] <- general_palette[1:15]

cols <- tmp_frac_colored$color
names(cols) <- tmp_frac_colored$sctype_cell_ann

dn_order <- c('donor1', 'donor2', 'donor3', 'donor4')
tmp$donor <- factor(tmp$donor, levels = rev(dn_order))
dn_comp <- ggplot(data = tmp) + 
  geom_mosaic(aes(x = product(clusters), fill = donor)) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 24),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 24)
  ) +
  scale_fill_manual(values = donor_palette, limits = dn_order)

tr_comp <- ggplot(data = tmp) + 
  geom_mosaic(aes(x = product(clusters), fill = treatment)) +
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
  geom_mosaic(aes(x = product(clusters), fill = donor_treatment)) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 24),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 24)
  ) +
  scale_fill_manual(values = donor_treatment_palette, limits = dntr_order)

phase_comp <- ggplot(data = tmp) + 
  geom_mosaic(aes(x = product(clusters), fill = Phase)) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 24),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 24)
  ) +
  scale_fill_manual(values = phase_palette)


singler_m_comp <- ggplot(data = tmp) + geom_mosaic(aes(x = product(clusters), fill=NK_T_main)) +
  ggtitle('SingleR annotation main labels') +
  scale_fill_manual(values=general_palette)

singler_f_comp <- ggplot(data = tmp) + geom_mosaic(aes(x = product(clusters), fill=NK_T_fine)) +
  ggtitle('SingleR annotation fine labels') +
  scale_fill_manual(values=general_palette) +
  theme(axis.text.y = element_blank(), text=element_text(size=20))

tmp$sctype_cell_ann_rev <- factor(tmp$sctype_cell_ann, levels = rev(levels(tmp$sctype_cell_ann)))
sctype_f_comp <- ggplot(data = tmp) + geom_mosaic(aes(x = product(clusters), fill=sctype_cell_ann_rev)) +
  scale_fill_manual(NULL, values=cols, limits = levels(tmp$sctype_cell_ann)) +
  ggtitle('ScType cell type annotation composition') +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x="Cluster", y='Cell type') + 
  theme(axis.text.y = element_blank(), text=element_text(size=20))

sctype_m_comp <- DimPlot(so, reduction = reduction, label = TRUE, label.size = 6, repel = TRUE, group.by = 'sctype_cluster_ann') +
  ggtitle('ScType annotation per cluster') +
  theme_minimal(base_size = 20) +
  theme(
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 24),
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 24)
  ) +
  scale_colour_manual(values = general_palette, limits = unique(tmp$sctype_cluster_ann))


comp_plot <- grid.arrange(dn_comp, tr_comp, dn_tr_comp, phase_comp, ncol = 2)
mysave(file_extn, paste0(dir_path,method_name, '_comp.', file_extn), comp_plot, height = 15, width = 25, grid=T)
mysave(file_extn, paste0(dir_path,method_name, '_singleR_m_comp.', file_extn), singler_m_comp, height = 8, width = 16)
mysave(file_extn, paste0(dir_path,method_name, '_singleR_f_comp.', file_extn), singler_f_comp, height = 8, width = 16)
mysave(file_extn, paste0(dir_path,method_name, '_sctype_m_comp.', file_extn), sctype_m_comp, height = 8, width = 16)
mysave(file_extn, paste0(dir_path,method_name, '_sctype_f_comp.', file_extn), sctype_f_comp, height = 8, width = 16)


# EXPRESSION PLOTS

# Literature markers
# ZBTB16->PLZF  TBX21->T-bet  RORC->ROR gamma t
literature_markers <- c('ZBTB16', 'TBX21', 'GATA3', 'RORC', 'IL4', 'IFNG')

exp_plot <- FeaturePlot(so, features=literature_markers, ncol=3) + theme_minimal(base_size = 20) 
vln_plot <- VlnPlot(so, features = literature_markers, split.by = 'clusters', ncol=3) & scale_fill_manual(values=general_palette)
mysave(file_extn, paste0(dir_path,method_name, '_lit_markers_vln.', file_extn), vln_plot, height = 10, width = 15)
mysave(file_extn, paste0(dir_path,method_name, '_lit_markers_exp.', file_extn), exp_plot, height = 10, width = 15)



# Antibodies
adt_labels <- paste0("adt_", c("CD122", "CD62L","CD314", "CD223", "CD279", "CD150", "CD319", "TCR"))

adt_exp_plot <- FeaturePlot(so, features = adt_labels, max.cutoff="q95", ncol=4) & scale_color_viridis_c()
adt_vln_plot <- VlnPlot(so, features = adt_labels, split.by = clusters, ncol=4) + theme_minimal(base_size = 20)

mysave(file_extn, paste0(dir_path,method_name, '_adt_vln.', file_extn), adt_vln_plot, height = 10, width = 20)
mysave(file_extn, paste0(dir_path,method_name, '_adt_exp.', file_extn), adt_exp_plot, height = 10, width = 20)



# Interleukin

il_exp_plot <- FeaturePlot(so, features=c('IL2', 'IL15')) + theme_minimal(base_size = 20)
il_vln_plot <- VlnPlot(so, features = c('IL2', 'IL15'), split.by = clusters) + theme_minimal(base_size = 20)

mysave(file_extn, paste0(dir_path,method_name, '_IL_vln.', file_extn), il_vln_plot, height = 8, width = 16)
mysave(file_extn, paste0(dir_path,method_name, '_IL_exp.', file_extn), il_exp_plot, height = 8, width = 16)


# CLUSTER MARKERS
cluster_markers <- read.csv('new_data/integrated_rpca_markers05_min075.csv')

genes <- cluster_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_max(avg_log2FC, n = 3) %>%
  ungroup()

unique_genes <- genes$gene
# write.csv(genes, 'data/i05_markers_min075_top20.csv')

for (i in seq(1, length(unique_genes), by = 3)){
  cluster = (i-1)/3
  selected_features <- unique_genes[i:(i+2)]

  vln_plot <- VlnPlot(so, features=selected_features, ncol=3, pt.size=0) & scale_fill_manual(values=general_palette)

  mysave(file_extn, paste0(dir_path,method_name, '_markers_vln_c', cluster, '.', file_extn), vln_plot, height = 5, width = 15)
}



nkg2d_exp <- FeaturePlot(so, features = 'KLRK1', alpha=0.6) + ggtitle('NKG2D expression') + theme_minimal(base_size = 20)
nkg2d_vln <- VlnPlot(so, features = 'KLRK1', pt.size=0) + ggtitle('NKG2D expression density in clusters')

mysave(file_extn, paste0(dir_path,method_name, '_NKG2D_exp.', file_extn), nkg2d_exp, height = 8, width = 12)
mysave(file_extn, paste0(dir_path,method_name, '_NKG2D_vln.', file_extn), nkg2d_vln, height = 8, width = 12)


p <- ggplot(tmp, aes(x = treatment, fill = treatment)) +
  geom_bar(position = "dodge", stat = "count") +
  scale_fill_manual(values=treatment_palette, limits=c('IL2', 'IL15')) +
  labs(x = "treatment", y = "Count") +
  theme_minimal(base_size = 20)


p <- ggplot(tmp, aes(x = donor, fill=donor)) +
  geom_bar(position = "dodge", stat = "count") +
  scale_fill_manual(values=donor_palette) +
  labs(x = "Donors", y = "Count") +
  theme_minimal(base_size = 20)



# inkt1:
p <- VlnPlot(so, features=c("KIR3DL1", "LAG3", "TYROBP"), ncol=3, pt.size=0) & scale_fill_manual(values=general_palette)
mysave('pdf', 'new_figures/inkt1_p1.pdf', p, height=5, width=15)
p <- VlnPlot(so, features=c('CCL3','CCL4','CCL5'), ncol=3, pt.size=0) & scale_fill_manual(values=general_palette)
mysave('pdf', 'new_figures/inkt1_p2.pdf', p, height=5, width=15)

# inkt2
p <- VlnPlot(so, features=c('TNFRSF4',  'TNFRSF18'), ncol=2, pt.size=0)& scale_fill_manual(values=general_palette)
mysave('pdf', 'new_figures/inkt2_p1.pdf', p, height=5, width=10)
p <- VlnPlot(so, features=c('IL2RA', 'IL17RB'), ncol=2, pt.size=0)& scale_fill_manual(values=general_palette)
mysave('pdf', 'new_figures/inkt2_p2.pdf', p, height=5, width=10)

#inkt17
p <- VlnPlot(so, features=c('RUNX2', 'CAMK4', 'EXOC6B'), ncol=3, pt.size=0)& scale_fill_manual(values=general_palette)
mysave('pdf', 'new_figures/inkt17_p1.pdf', p, height=5, width=15)
p <- VlnPlot(so, features=c('CRY1', 'ZFP36L1', 'MYO1D'), ncol=3, pt.size=0)& scale_fill_manual(values=general_palette)
mysave('pdf', 'new_figures/inkt17_p2.pdf', p, height=5, width=15)
  
# steady-state
p <- VlnPlot(so, features=c('IL7R', 'FOXP1', 'BTG1'), ncol=3, pt.size=0)& scale_fill_manual(values=general_palette)
mysave('pdf', 'new_figures/inkt_ss_p1.pdf', p, height=5, width=15)
p <-  VlnPlot(so, features=c('RIPOR2', 'PDCD4'), ncol=2, pt.size=0)& scale_fill_manual(values=general_palette)
mysave('pdf', 'new_figures/inkt_ss_p2.pdf', p, height=5, width=10)
p <- VlnPlot(so, features=c('KLF2', 'STK38', 'MAF'), ncol=3, pt.size=0)& scale_fill_manual(values=general_palette)
mysave('pdf', 'new_figures/inkt_ss_p3.pdf', p, height=5, width=15)
p <- VlnPlot(so, features=c('IER2','LIME1'), ncol=2, pt.size=0)& scale_fill_manual(values=general_palette)
mysave('pdf', 'new_figures/inkt_ss_p4.pdf', p, height=5, width=10)
p <- VlnPlot(so, features=c('CD69','adt_CD62L'), ncol=2, pt.size=0)& scale_fill_manual(values=general_palette)
mysave('pdf', 'new_figures/inkt_ss_p5.pdf', p, height=5, width=10)

