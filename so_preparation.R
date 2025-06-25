library(Seurat)
library(dplyr)
library(tidyr)
library(SingleR)
library(presto)
library(HGNChelper)
library(openxlsx)
library(ramify)

# OPTIONS
so <- readRDS('path/integrated.rds')
assay <-'SCT' # SCT
so <- integrated

# Cell cycle scoring
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
so <- CellCycleScoring(so, s.features = s_genes, g2m.features = g2m_genes)

# SingleR annotation
hpca_ref <- celldex::HumanPrimaryCellAtlasData()

singleR_fine <- SingleR(test = so[['SCT']]$data, # Ważne! Do annotacji singleR NIE UŻYWAĆ scalowanych wartości
                       ref = hpca_ref,
                       labels = hpca_ref$label.fine,
                       assay.type.test=1
)
so[["SingleR_main"]]<- replace(singleR_main$pruned.labels, is.na(singleR_main$pruned.labels), 'Unassigned')
# so[["SingleR_fine"]]<- replace(singleR_fine$pruned.labels, is.na(singleR_fine$pruned.labels), 'Unassigned')
message('NK/T annotation finished')

cell_annotation <- data.frame(
  row.names = colnames(so),
  clusters = so[["clusters"]][,1],
  celltype = so[["SingleR_fine"]][,1],
  stringsAsFactors = F
)

cluster_annotation <- cell_annotation %>%
  group_by(clusters, celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(clusters) %>%
  slice_max(order_by = count, n = 1, with_ties = FALSE) %>%
  ungroup()

res_df <- merge(cell_annotation, cluster_annotation, by = "clusters")
rownames(res_df) <- rownames(cell_annotation)

so[['sR_cluster_annotation']] <- res_df['celltype.y']
colnames(so@meta.data)[colnames(so@meta.data) == "cell_annotation"] <- "scT_cell_annotation"

# ScType annotation

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
gs_list <- gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

scRNAseqData_scaled <- GetAssayData(so, layer = "scale.data")

# run ScType
es_max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = T, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)


cL_resutls <- do.call("rbind", lapply(unique(so@meta.data$clusters), function(cl){
  es_max.cl = sort(rowSums(es_max[ ,rownames(so@meta.data[so@meta.data$clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es_max.cl), scores = es_max.cl, ncells = sum(so@meta.data$clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unassigned"

so@meta.data$cluster_annotation <- ""

for(j in unique(sctype_scores$cluster)){
  cl_type <- sctype_scores[sctype_scores$cluster==j,]; 
  so@meta.data$cluster_annotation[so@meta.data$clusters == j] = as.character(cl_type$type[1])
}
so[['cell_annotation']] <- rownames(es_max)[argmax(es_max, rows = F)]

saveRDS(so, "path/annotated.rds")

# Markers

so <- PrepSCTFindMarkers(so)
markers_df <-FindAllMarkers(so, only.pos = T, min.pct = 0.75)
write.csv(markers_df, 'new_data/integrated_harmony_markers05_min075.csv')

saveRDS(so, 'new_data/integrated_harmony_preped05.rds')
