import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad


scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2


adata_path = '/home/jd438446/iNKT_project/scvelo/data/i05_adata_velocity.h5ad'
adata = sc.read_h5ad(adata_path)


dn_tr_labels = list(adata.obs['donor_treatment'].unique())
dn_tr_adatas = [adata[adata.obs['donor_treatment'] == label] for label in dn_tr_labels]

a_il2 = adata[adata.obs['treatment'] == 'IL2']
a_il15 = adata[adata.obs['treatment'] == 'IL15']



scv.pl.velocity_embedding_grid(adata, basis='umap', color='clusters_rpca', save='cluster_grid.png', title='RNA velocity grid', scale=0.25)
scv.pl.velocity_embedding_stream(adata, basis='umap', color='clusters_rpca', save='cluster_stream.png', title='RNA velocity stream')

scv.pl.velocity_embedding_stream(a_il2, basis='umap', color='clusters_rpca', save='il2_velocity.png', title='RNA velocity stream for IL2 treatment')
scv.pl.velocity_embedding_stream(a_il15, basis='umap', color='clusters_rpca', save='il15_velocity.png', title='RNA velocity stream for IL15 treatment')

for label, adata_obj in zip(dn_tr_labels, dn_tr_adatas):
  scv.pl.velocity_embedding_stream(adata_obj, basis='umap', color='clusters_rpca', save=f'{label}_velocity.png', title=f'RNA velocity stream for {label}')
