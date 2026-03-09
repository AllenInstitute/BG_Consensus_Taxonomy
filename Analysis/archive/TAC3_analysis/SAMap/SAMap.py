from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import scanpy as sc
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import os
import glob

##
cxg_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/xspecies"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/samap"
os.chdir(analysis_dir)

##
h5ad_files = glob.glob(f"{analysis_dir}/*for_samap_pr.h5ad")

## --------------------------------------------
## Preprocess data for SAMap
## Move layers["UMIs"] to X
# for filename in h5ad_files:
#     adata = ad.read_h5ad(filename)
#     if "UMIs" in adata.layers.keys():
#         print("Processing", filename)
#         adata.X = adata.layers["UMIs"].copy()
#         del adata.layers["UMIs"]
#         adata.write_h5ad(filename)

## Read data with SAMap
# sams = dict(zip(["hu", "mq", "ms", "mm"], h5ad_files))  ## Human, Macaque, Marmoset, Mouse

## --------------------------------------------
## Load in SAM objects
sams = {}
for file, species in zip(h5ad_files, ["hu", "mm", "ms", "mq"]):  ## Human, Macaque, Marmoset, Mouse
    print(file)
    sams[species] = SAM()
    sams[species].load_data(file)

## Annotations
keys = {"hu": "Cluster", "mm": "AIT21.cluster", "ms": "Cluster", "mq": "Cluster"}

## Run SAMap
sm = SAMAP(sams, keys=keys, f_maps = 'maps/')

## Run alignment
sm.run(pairwise=True, NUMITERS=3)
samap = sm.samap

##
sm.query_gene_pairs('hu_LHX8')

sm.scatter()
plt.savefig("/home/nelson.johansen/myplot.png")


D,MappingTable = get_mapping_scores(sm, keys, n_top = 0)

## Remove log1p from uns to fix errors
del sm.sams["hu"].adata.uns["log1p"]
del sm.sams["mq"].adata.uns["log1p"]
del sm.sams["ms"].adata.uns["log1p"]
gpf = GenePairFinder(sm, keys=keys)
gene_pairs = gpf.find_all(align_thr=0.2)

gene_pairs = gpf.find_all(n='hu_STR-BF TAC3-PLPP4-LHX8 GABA', align_thr=0.1)

##
adatas = []
for species in sm.sams.keys():
    adatas.append(sm.sams[species].adata)

combined = ad.concat(adatas, label="species", keys=list(sm.sams.keys()), index_unique=None)
combined.X = combined.layers["X_disp"].copy()

## --------------------------------------------
## Preprocess data for SAMap
h5ad_files = glob.glob(f"{analysis_dir}/*for_samap.h5ad")
adatas = []
for filename in h5ad_files:
    adatas.append(ad.read_h5ad(filename))

##
combined_orig = ad.concat(adatas, label="species", keys=["hu", "mq", "ms", "mm"], index_unique=None, join="outer", merge="unique")

## 
combined_orig.obsm["X_samap"] = combined[combined_orig.obs_names].obsm["X_umap_samap"].copy()
sc.pp.normalize_total(combined_orig, target_sum=1e6)
sc.pp.log1p(combined_orig)

##
combined_orig.write_h5ad(os.path.join(cxg_dir, "combined_for_samap_pr.h5ad"))