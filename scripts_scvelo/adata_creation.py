import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

counts_path = '/home/jd438446/iNKT_project/scvelo/data/norm_counts.mtx'
meta_path = '/home/jd438446/iNKT_project/scvelo/data/metadata.csv'
pca_path = '/home/jd438446/iNKT_project/scvelo/data/pca.csv'
umap_path = '/home/jd438446/iNKT_project/scvelo/data/umap.csv'
features_path = '/home/jd438446/iNKT_project/scvelo/data/features.csv'

# load sparse matrix:
X = io.mmread(counts_path)

# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv(meta_path)

# load gene names:
with open(features_path, 'r') as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.obs['clusters_rpca'] = adata.obs['clusters_rpca'].astype('category')
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv(pca_path)
umap = pd.read_csv(umap_path)

pca.index = adata.obs.index
umap.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = umap.to_numpy()


# save dataset as anndata format
adata.write('data/integrated08_adata.h5ad')
