import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")  ## non-GUI backend that supports tostring_rgb()
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import to_hex

from tqdm import tqdm
from pathlib import Path
import keras
import crested

import tensorflow as tf
import anndata as ad
import os
import pysam
import scanpy as sc

from scipy.stats import spearmanr, pearsonr
import seaborn as sns

## -------------------------
variant_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalGanglia/scripts/EvoRegionhDEGsEnrichment/Locations"
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"
model_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested"
cxg_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/"

## -------------------------
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)
cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

## Create color map
color_map = dict(zip(cluster_meta.Group, cluster_meta.color_hex_group))

## -------------------------
models = {
    "human": tf.keras.models.load_model(os.path.join(model_dir, "human/finetune/12.keras"), compile=False),
    "macaque": tf.keras.models.load_model(os.path.join(model_dir, "macaque/finetune/11.keras"), compile=False)
}

genomes = {
    "human": pysam.FastaFile("/allen/programs/celltypes/workgroups/rnaseqanalysis/references/human/10x/grch38.p2/genome/fasta/genome.fa"),
    "macaque": pysam.FastaFile("/allen/programs/celltypes/workgroups/rnaseqanalysis/references/macaque/ncbi/mmul10/genome/fasta/genome.fa")
}


adatas = {
    "human": ad.read_h5ad(os.path.join(data_dir, "human/crested_adata/human_basalganglia_hmba_pre-print_finetune_crested.h5ad")),
    "macaque": ad.read_h5ad(os.path.join(data_dir, "macaque/crested_adata/macaque_basalganglia_hmba_pre-print_finetune_crested.h5ad"))
}

## -------------------------
## Gather top 100 peaks per species
for species in adatas.keys():
    print(f"Processing {species}...")
    crested.pp.sort_and_filter_regions_on_specificity(
        adatas[species], top_k=100, method="proportion"
    )
    adatas[species].var["sequence"] = [genomes[species].fetch(chrom, start, end).upper() for chrom, start, end in tqdm(zip(adatas[species].var['chr'], adatas[species].var['start'], adatas[species].var['end']))]

## Process cross-species embeddings
model_embedding_dict = {"human": {}, "macaque": {}} 
for species in models.keys():
    model = models[species]
    model_embedding_dict[species] = {}
    for target_species in adatas.keys():
        print(f"Predicting {target_species} using model trained on {species}...")
        adata = adatas[target_species]
        ## -------------------------
        model_embedding = crested.tl.extract_layer_embeddings(
            input = adata.var["sequence"].tolist(),
            model = model,
            layer_name = "global_average_pooling1d",
        )
        ## Store embeddings in adata
        model_embedding_dict[species][target_species] = model_embedding


species = "human"
target_species = "macaque"

print(f"Storing embeddings: model from {species}, data from {target_species}...")

adatas[species].var["species"] = species
adatas[target_species].var["species"] = target_species

##
embedding_adata = ad.AnnData(
    X=None,
    obs = pd.concat([adatas[species].var] + [adatas[target_species].var]),
    var = pd.concat([adatas[species].obs] + [adatas[target_species].obs]),
)
## Merge embeddings from species against both target species
embedding_adata.obsm[f"X_embedding_from_{species}_model"] = np.vstack([
    model_embedding_dict[species][target_species] 
    for target_species in adatas.keys()
])

embedding_adata.obs["celltype"] = embedding_adata.obs["Class name"]
for tax_level in ["Neighborhood", "Class", "Subclass"]:
    embedding_adata.obs[tax_level.lower()] = embedding_adata.obs["celltype"].map(dict(zip(cluster_meta.Group, cluster_meta[tax_level])))

## -------------------------
## Plotting
plot_adata = embedding_adata[embedding_adata.obs["neighborhood"] != "Nonneuronal", :].copy()

sc.pp.neighbors(plot_adata, use_rep=f"X_embedding_from_{species}_model")
sc.tl.umap(plot_adata, min_dist=1.0, random_state=42)

## Plot UMAP
colors = plot_adata.obs["Class name"].map(color_map).values.astype(str)
colors[pd.isna(colors)] = "lightgrey"  # Handle NaN values

# Define marker shapes for each species
marker_map = {
    "human": "o",       # circle
    "macaque": "s",     # square
    "mouse": "^",       # triangle
    "marmoset": "D",    # diamond
}

## Plot UMAP colored by cell type
# Loop over species and plot each with different marker
for sp in np.unique(plot_adata.obs.species):
    plt.figure(figsize=(12, 12))
    idx = plot_adata.obs.species == sp
    plt.scatter(
        plot_adata.obsm["X_umap"][idx, 0],
        plot_adata.obsm["X_umap"][idx, 1],
        c=np.array(colors)[idx],          # colors already mapped to cell types
        s=15,
        alpha=0.8,
        marker=marker_map.get(sp, "o"),   # default to 'o' if species not in map
        label=sp
    )
    ## Add color legend
    plt.title("UMAP of CRESted embedding of top 100 Human peaks")
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.savefig(os.path.join(analysis_dir, "crested", "figures", f"{species}_xspecies_{sp}_crested_umap.png"), dpi=900)







# ## Add to cellxgene as anndata
# embedding_adata = ad.AnnData(
#     X=predictions,
#     obs = adata_top.var.copy(),
#     var = pd.DataFrame(index=adata_top.obs_names),
#     obsm={"X_embedding": model_embedding},
#     layers={"accessibility": adata_top.X.T},
# )
# embedding_adata.obs_names_make_unique()
# embedding_adata.obs["celltype_prediction_gini"] = calc_gini(embedding_adata.X).max(axis=1)
# embedding_adata.obs["celltype_prediction"] = embedding_adata.var_names[embedding_adata.X.argmax(axis=1)]

# for tax_level in ["Neighborhood", "Class", "Subclass"]:
#     embedding_adata.obs[tax_level.lower()] = embedding_adata.obs["celltype"].map(dict(zip(cluster_meta.Group, cluster_meta[tax_level])))
#     embedding_adata.obs[f"{tax_level.lower()}_pred"] = embedding_adata.obs["celltype_prediction"].map(dict(zip(cluster_meta.Group, cluster_meta[tax_level])))