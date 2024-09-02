#!/bin/bash

ann_path='/home/jd438446/iNKT_project/scvelo/data/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
ranger_out='/home/ajank/iNKT_differentiation/data/CITEseq/'

# velocyto run10x $ranger_out $ann_path'pool1'
velocyto run10x $ranger_out'pool2' $ann_path
velocyto run10x $ranger_out'pool3' $ann_path
velocyto run10x $ranger_out'pool4' $ann_path
