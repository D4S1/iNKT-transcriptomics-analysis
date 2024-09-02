import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

loom_dir = '/home/ajank/iNKT_differentiation/data/CITEseq/'
adata_path = '/home/jd438446/iNKT_project/scvelo/data/i05_adata.h5ad'

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2

adata = sc.read_h5ad(adata_path)

# load loom files for spliced/unspliced matrices for each sample:
pools_labels = ['P1_', 'P2_', 'P3_', 'P4_']
looms = [
          scv.read(f'{loom_dir}pool1/velocyto/pool1.loom', cache=True),
          scv.read(f'{loom_dir}pool2/velocyto/pool2.loom', cache=True),
          scv.read(f'{loom_dir}pool3/velocyto/pool3.loom', cache=True),
          scv.read(f'{loom_dir}pool4/velocyto/pool4.loom', cache=True),
]

for label, loom in zip(pools_labels, looms):
  # rename barcodes in order to merge:
  barcodes = [f"{label}{bc.split(':')[1][:-1]}-1" for bc in loom.obs.index.tolist()]
  loom.obs.index = barcodes
  # make variable names unique
  loom.var_names_make_unique()


# concatenate the three loom
loom = looms[0].concatenate(looms[1:])

# merge matrices into the original adata object
adata = scv.utils.merge(adata, loom)
print('Looms merge complited')

# pre-process
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
print('Moments calculation complited')

# compute velocity
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)
print('Velocity calculation complited')

adata.write('data/i80_adata_velocity.h5ad')
