library(Seurat)
library(dplyr)
library(tidyr)
library(SingleR)
library(presto)
library(HGNChelper)
library(openxlsx)
library(ramify)

# OPTIONS
so <- readRDS('data/integrated_so05.rds')
assay <-'SCT' # SCT


# Original and regressed effect objects creation

so <- readRDS(file)
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
so <- CellCycleScoring(so, s.features = s_genes, g2m.features = g2m_genes)

# SingleR annotation
hpca_ref <- celldex::HumanPrimaryCellAtlasData()

singleR_main <- SingleR(test = so[[assay]]$data, # Ważne! Do annotacji singleR NIE UŻYWAĆ scalowanych wartości
                         ref = hpca_ref,
                         labels = hpca_ref$label.main,
                         assay.type.test=1
  )
singleR_fine <- SingleR(test = so[[assay]]$data, # Ważne! Do annotacji singleR NIE UŻYWAĆ scalowanych wartości
                       ref = hpca_ref,
                       labels = hpca_ref$label.fine,
                       assay.type.test=1
)
so[["NK_T_main"]]<- replace(singleR_main$pruned.labels, is.na(singleR_main$pruned.labels), 'unassigned')
so[["NK_T_fine"]]<- replace(singleR_fine$pruned.labels, is.na(singleR_fine$pruned.labels), 'unassigned')
message('NK/T annotation finished')

# ScType annotation

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
gs_list <- gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

scRNAseqData_scaled <- as.matrix(so[["SCT"]]@scale.data)

# run ScType
es_max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive)


cL_resutls <- do.call("rbind", lapply(unique(so@meta.data$clusters_rpca), function(cl){
  es_max.cl = sort(rowSums(es_max[ ,rownames(so@meta.data[so@meta.data$clusters_rpca==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es_max.cl), scores = es_max.cl, ncells = sum(so@meta.data$clusters_rpca==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])

so@meta.data$sctype_cluster_ann <- ""

for(j in unique(sctype_scores$cluster)){
  cl_type <- sctype_scores[sctype_scores$cluster==j,]; 
  so@meta.data$sctype_cluster_ann[so@meta.data$clusters_rpca == j] = as.character(cl_type$type[1])
}
so[['sctype_cell_ann']] <- rownames(es_max)[argmax(es_max, rows = F)]

# Markers

so <- PrepSCTFindMarkers(so)
markers_df <-FindAllMarkers(so, only.pos = T, min.pct = 0.75)
write.csv(markers_df, 'data/integrated_markers05_min075.csv')

saveRDS(so, 'data/integrated_preped_so05.rds')
