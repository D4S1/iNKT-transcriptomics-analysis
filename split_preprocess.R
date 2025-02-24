library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Matrix)
library(gridExtra)

# OPTIONS
path <- "/home/ajank/iNKT_differentiation/data/CITEseq/" # to dir with data (pools directories)
project_dir <- "/home/jd438446/iNKT_project/"


# Directory setting
setwd(path)


# Loadning H5 files
filtered_h5 <- list.files(path, pattern="filtered.*\\.h5", recursive = TRUE)
filtered_h5 <- lapply(filtered_h5, Read10X_h5)

# filtered_h5[[1]][[1]]
# first [[i]] indicate pool number
# second [[j]] j in {1,2} indicates gene expression (1) or Antibodies (2)
# another way is to use filtered_h5$'Gene Expression' and filtered_h5$'Custom'

# get the donors_id filenames (1 per pool). In those file there are mapping barcode -> donor 0 / donor 1
pool_donors <- list.files(path, pattern='ids.tsv', recursive = TRUE)
pool_donors <- lapply(pool_donors, function(path){
  # Function takes path and returns data frame
  # were row are cell barcodes and column keeps 
  # donors (pool notation) for each cell
  
  donors_ids <- read.csv(path, sep="\t", header=TRUE)[1:2]
  # donors_ids has columns: cell, donor_id
  
  donors_ids <- filter(donors_ids, donor_id == 'donor0' | donor_id == 'donor1')
  donors_ids <- data.frame(pool_donor=donors_ids$donor_id, row.names=donors_ids$cell)
  return (donors_ids)
})

setwd(project_dir)

# ======CREATING ANNOTATION======

# in annotation we want: donor, treatment, pool, donor_treatment

# columns: pool, pool_donor, donor, treatment
donors_translator_df <- read.csv('donor_translator.csv', sep=';', header=TRUE)

annotations <- lapply(1:4, function(pool_num){
  # Function takes pool number(int)
  # Returns merged data frame row: cell, columns: pool, donor, treatment
  pool_annotation <- filter(donors_translator_df, pool==paste0('pool', pool_num))
  
  an_frame <-data.frame(cell=rownames(pool_donors[[pool_num]]),pool_donor=pool_donors[[pool_num]]$pool_donor)
  annotation <- merge(an_frame, pool_annotation, by='pool_donor', all.x = TRUE)
  annotation <- annotation[order(annotation$cell),]
  rownames(annotation) <- an_frame$cell
  annotation[['donor_treatment']] <- paste(annotation[['treatment']], annotation[['donor']])
  return (annotation[c(-1, -2)])
})


# Match matrices to annotation
filtered_h5 <- mapply(function(pool_data, annotation){
  # check that all cells listed in rows of annotation are represented in pool date
  stopifnot(length(setdiff(rownames(annotation), colnames(pool_data$'Gene Expression'))) == 0)
  
  # omits rows from pool data that are not present in annotation
  names(pool_data) <- c('Gene Expression', 'ADT')
  pool_data$'Gene Expression' <- pool_data$'Gene Expression'[,rownames(annotation)]
  pool_data$'ADT' <- pool_data$'ADT'[,rownames(annotation)]
  return(pool_data)
}, filtered_h5, annotations, SIMPLIFY = FALSE)


# language R
# given a list of filtered_h5 and annotation (df with columns: donor, treatment, pool, donor_treatment)
# we want to split annotation 2 two df based on treatment values
# and then we want to split h5 object same way so barcodes from each treatment
# are in one h5 object
# and return 2 vectors one of them keeping h5 objects splited by treatments
# and 2nd of them keeping the annotations

split_pools <- function(filtered_h5s, annotations){
  # Initialize lists to hold the split h5 objects and annotations
  split_h5 <- c()
  split_annotations <- c()
  
  # Iterate through each pair of h5 and annotation
  for (i in seq_along(annotations)) {
    
    # Split the annotation by treatment
    annotation <- annotations[[i]]
    h5 <- filtered_h5s[[i]]
    treatments <- unique(annotation$treatment)
    
    anns <- lapply(treatments, function(treatment){
      return(annotation[annotation$treatment == treatment,])
    })
    
    h5s <- lapply(anns, function(ann){
      h5$'Gene Expression' <- h5$'Gene Expression'[,rownames(ann)]
      h5$'ADT' <- h5$'ADT'[,rownames(ann)]
      return(h5)
    })
    
    split_h5 <- c(split_h5, h5s)
    split_annotations <- c(split_annotations, anns)
  }
  
  # Return the split h5 objects and annotations
  return(list(h5_objects = split_h5, annotations = split_annotations))
}

res <- split_pools(filtered_h5, annotations)
filtered_h5 <- res$h5_objects
annotations <- res$annotations

# ======CREATING SEURAT OBJECTS======

so_list <- mapply(function(pool_data, annotation){
  so <- CreateSeuratObject(
    pool_data$'Gene Expression',
    meta.data = annotation,
    project = annotation$donor_treatment[[1]]
  )
  so[["ADT"]] <- CreateAssayObject(counts=pool_data$'ADT')
  return (so)
}, filtered_h5, annotations)

saveRDS(so_list, file="new_data/so_list_b4qc.rds")

# so_list <- readRDS("data/so_list_b4qc.rds")


# ======QUALITY CONTROL======

# plot before filtering
b4_nrows <- lapply(so_list, ncol)

so_list <- lapply(so_list, function(pool_so){
  pool_so <- PercentageFeatureSet(pool_so, "^MT-", col.name = "percent_mito")
  pool_so <- PercentageFeatureSet(pool_so, "^RP[SL]", col.name = "percent_ribo")
  return (pool_so)
})

create_plot <- function(so_list, filename){
  plot_list <- lapply(so_list, function(so){
    vln_mito <- VlnPlot(
      so, 
      features = "percent_mito", 
      split.by = 'donor',
      pt.size=0,
      split.plot = TRUE,
    )
    vln_ribo <- VlnPlot(
      so, 
      features = "percent_ribo", 
      split.by = 'donor',
      pt.size=0,
      split.plot = TRUE,
    )
    scatter_plot <- FeatureScatter(
      so,
      "nCount_RNA",
      "nFeature_RNA",
      group.by = "donor",
      plot.cor = FALSE,
      pt.size = 0.5,
    )
    return (grid.arrange(vln_mito, vln_ribo, scatter_plot, ncol = 3))
    
  })
  final_plot <- do.call(grid.arrange, c(plot_list, ncol = 1))
  ggsave(paste('graphics/', filename, '.png', sep='' ), final_plot, width = 15, height = 15)
  return (final_plot)
}


# b4qc_plot <- create_plot(so_list, 'b4qc_plot')

so_list <- lapply(so_list, function(so){
  max_mito <- quantile(so@meta.data$percent_mito, 0.95)
  max_ribo <- quantile(so@meta.data$percent_ribo, 0.95)
  min_nCount <- quantile(so@meta.data$nCount_RNA, 0.05)
  max_nCount <- quantile(so@meta.data$nCount_RNA, 0.95)
  min_nFeature <- quantile(so@meta.data$nFeature_RNA, 0.05)
  
  so <- subset(so, percent_mito <= max_mito & 
                 percent_ribo <= max_ribo &
                 nCount_RNA >= min_nCount &
                 nCount_RNA <= max_nCount &
                 nFeature_RNA >= min_nFeature
  )
  return (so)
})

# plot after filtering
# afqc_plot <- create_plot(so_list, 'afqc_plot')
af_nrows <- lapply(so_list, ncol)

cat("statistics after filtering",file="qc_stats.txt",sep="\n")
for (i in 1:length(so_list)){
  cat(paste("Pool: ", i, sep=' '),file="qc_stats.txt",append=TRUE)
  cat(paste("Number of Cells before Filtering: ", b4_nrows[[i]], sep=' '),file="qc_stats.txt",append=TRUE)
  cat(paste("Number of Cells after Filtering: ", af_nrows[[i]], sep=' '),file="qc_stats.txt",append=TRUE)
  percent = round((b4_nrows[[i]] - af_nrows[[i]]) / b4_nrows[[i]] * 100, 1)
  cat(paste("% of filtered out: ", percent, "%", sep=' '),file="qc_stats.txt",append=TRUE)
}

saveRDS(so_list, file='new_data/so_list_afqc.rds')
