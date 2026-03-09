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
models = [
    tf.keras.models.load_model(os.path.join(model_dir, "human/finetune/12.keras"), compile=False),
    tf.keras.models.load_model(os.path.join(model_dir, "macaque/finetune/11.keras"), compile=False)
]

adatas = [
    ad.read_h5ad(os.path.join(data_dir, "human/crested_adata/human_basalganglia_hmba_pre-print_finetune_crested.h5ad")),
    ad.read_h5ad(os.path.join(data_dir, "macaque/crested_adata/macaque_basalganglia_hmba_pre-print_finetune_crested.h5ad"))
]

genomes = {
    "human": pysam.FastaFile("/allen/programs/celltypes/workgroups/rnaseqanalysis/references/human/10x/grch38.p2/genome/fasta/genome.fa"),
    "macaque": pysam.FastaFile("/allen/programs/celltypes/workgroups/rnaseqanalysis/references/macaque/ncbi/mmul10/genome/fasta/genome.fa")
}

## -------------------------
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)
cluster_meta["Group"].replace(" ", "_", regex=True, inplace=True)

## Create color map
color_map = dict(zip(cluster_meta.Group, cluster_meta.color_hex_group))

## -------------------------
cross_preds_human = {}
cross_preds_macaque = {}
cross_preds = [cross_preds_human, cross_preds_macaque]
peak_values = {}

species = ["human", "macaque"]

for adata_, sp in zip(adatas, species):
    adata_filtered = adata_.copy()
    crested.pp.sort_and_filter_regions_on_specificity(
        adata_filtered, top_k=500, method="proportion"
    )
    genome = genomes[sp]
    ## Fetch sequences
    seqs = [genome.fetch(chrom, start, end).upper() for chrom, start, end in tqdm(zip(adata_filtered.var['chr'], adata_filtered.var['start'], adata_filtered.var['end']))]
    ## Predict with all models
    for model, cr_preds in zip(models, cross_preds):
        cr_preds[sp] = {}
        peak_values[sp] = {}
        preds = crested.tl.predict(seqs, model)
        multiplier = 0
        for cls in list(adata_filtered.obs_names):
            cr_preds[sp][cls] = preds[multiplier*top_k:(multiplier+1)*top_k]
            peak_values[sp] = adata_filtered.X.T[multiplier*top_k:(multiplier+1)*top_k]
            multiplier+=1


## -------------------------
species_data = {
    'human':   {'cross_preds': cross_preds_human,   'adata': adatas[0]},
    'macaque': {'cross_preds': cross_preds_macaque, 'adata': adatas[1]},
}

_,_, corr_combined_moma, mouse_labels, macaque_labels = compare_species_correlations(
    species_pair=['human', 'macaque'],
    reordered_classes=cluster_meta["Group"].tolist(),
    species_data=species_data,
    method='spearman'
)

# --- Build DataFrames from numpy arrays ---
corr_mouse_macaque = pd.DataFrame(corr_combined_moma, index=mouse_labels, columns=macaque_labels)

# --- Melt each matrix to long form ---
def melt_corr(df, species1, species2):
    df = df.copy()
    df = df.stack().reset_index()
    df.columns = ['CellType1', 'CellType2', 'Correlation']
    df['SpeciesPair'] = f'{species1}–{species2}'
    return df

df_mm = melt_corr(corr_mouse_macaque, 'Human', 'Macaque')

# --- Combine and filter to only matching cell types (diagonal) ---
all_corrs = df_mm
all_corrs = all_corrs[all_corrs['CellType1'] == all_corrs['CellType2']].copy()
all_corrs['CellType'] = all_corrs['CellType1']

# --- Set plotting order ---
species_order = ['Human–Macaque']
celltype_order = [ct for ct in cluster_meta["Group"].tolist() if ct in all_corrs['CellType'].unique()]
# all_corrs['SpeciesPair'] = pd.Categorical(all_corrs['SpeciesPair'], categories=species_order, ordered=True)
all_corrs['CellType'] = pd.Categorical(all_corrs['CellType'], categories=celltype_order, ordered=True)

# --- Prepare color palette and mapping ---
palette = sns.color_palette('Set2', n_colors=3)
species_to_color = dict(zip(species_order, palette))

# --- Begin plot ---
plt.figure(figsize=(max(6, 0.2 * len(celltype_order)), 5))

# Map cell types to numeric positions on x-axis
x_pos = np.arange(len(celltype_order))
celltype_to_x = dict(zip(celltype_order, x_pos))

# Plot line + dot per species pair
for species_pair in species_order:
    df = all_corrs[all_corrs['SpeciesPair'] == species_pair].dropna()
    df = df.sort_values('CellType')
    x = [celltype_to_x[ct] for ct in df['CellType']]
    y = df['Correlation'].values
    color = species_to_color[species_pair]
    # Line
    plt.plot(x, y, label=species_pair, color=color, linewidth=1.5, alpha=0.7)
    # Dots
    plt.scatter(x, y, color=color, s=40, zorder=3)

# --- Final plot aesthetics ---
plt.xticks(ticks=x_pos, labels=celltype_order, rotation=90)
plt.ylim(0, 1.05)
plt.xlabel('Cell Type')
plt.ylabel('Sequence Prediction Spearman Correlation')
plt.title('Per Cell Type Correlation Across Species Pairs')
plt.legend(title='Species Pair', bbox_to_anchor=(1.02, 1), loc='upper left')
plt.grid(True, axis='y', linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(analysis_dir, "crested", "figures",'xspecies_model_corr.pdf'), bbox_inches='tight')








## -------------------------
## -------------------------
def compute_multi_correlation_matrices(cross_preds_a, cross_preds_b, 
                                       adata_a, adata_b, 
                                       species_key1='Human', 
                                       species_key2='Macaque',
                                       method='spearman', 
                                       log_transform=True,
                                       axis=1):
    """
    Computes three correlation matrices:
    - One for species_key1
    - One for species_key2
    - One for all classes combined across both species
    Parameters
    ----------
    cross_preds_a, cross_preds_b : dict
        Dictionaries of predictions for dataset A and B.
    adata_a, adata_b : AnnData
        AnnData objects with .obs_names to infer shape.
    species_key1, species_key2 : str
        The species keys to extract predictions.
    method : str
        'spearman' or 'pearson'.
    log_transform : bool
        Apply log1p transform to prediction values.
    axis : int
        Whether to extract correlations across columns (1) or rows (0).
    Returns
    -------
    dict of np.ndarray
        Dictionary with keys 'species1', 'species2', and 'combined' for the 3 correlation matrices.
    """
    def get_matrix(preds_a, preds_b, labels_a, labels_b, classes_a, classes_b):
        n_a, n_b = len(classes_a), len(classes_b)
        matrix = np.zeros((n_a, n_b))
        ##
        all_vals = np.zeros((len(classes_a), len(classes_b),2, 100)) ## FIX HARDCODING
        ##
        for i, cls_a in enumerate(preds_a):
            vals_a = preds_a[cls_a]
            vals_a = vals_a[:, i]
            ##
            if log_transform:
                vals_a = np.log1p(vals_a)
            ##
            for j, _ in enumerate(classes_b):
                vals_b = preds_b[cls_a]
                vals_b = vals_b[:, j]
                ##
                if log_transform:
                    vals_b = np.log1p(vals_b)
                ##
                if method == 'spearman':
                    corr, _ = spearmanr(vals_a, vals_b)
                else:
                    corr, _ = pearsonr(vals_a, vals_b)
                matrix[i, j] = corr
                all_vals[i,j,0] = vals_a
                all_vals[i,j,1] = vals_b
        ##
        return matrix, all_vals
    ##
    preds1_a = cross_preds_a[species_key1]
    preds1_b = cross_preds_b[species_key1]
    corr1, all_vals_a = get_matrix(preds1_a, preds1_b, preds1_a.keys(), preds1_b.keys(), list(adata_a.obs_names), list(adata_b.obs_names))
    ##
    preds2_a = cross_preds_a[species_key2]
    preds2_b = cross_preds_b[species_key2]
    corr2, all_vals_b = get_matrix(preds2_b, preds2_a, preds2_b.keys(), preds2_a.keys(), list(adata_b.obs_names), list(adata_a.obs_names))
    ##
    n_a, n_b = len(list(adata_a.obs_names)), len(list(adata_b.obs_names))
    matrix = np.zeros((n_a, n_b))
    ##
    for i in range(n_a):
        for j in range(n_b):
            vals_a = np.concatenate((all_vals_a[i,j,0], all_vals_b[j,i,1]))
            vals_b = np.concatenate((all_vals_a[i,j,1], all_vals_b[j,i,0]))
            if method == 'spearman':
                corr, _ = spearmanr(vals_a, vals_b)
            else:
                corr, _ = pearsonr(vals_a, vals_b)
            matrix[i, j] = corr
    ##
    return {
        'species1': corr1,
        'species2': corr2,
        'combined': matrix
    }


def plot_correlation_heatmap(corr_matrix, y_labels, x_labels, species_pair, figsize=(12, 12), cmap='coolwarm'):
    """
    Plots a heatmap of the correlation matrix.
    Parameters
    ----------
    corr_matrix : np.ndarray
        The correlation matrix to visualize.
    species_key : str
        Key to access the correct species in both prediction dicts.
    figsize : tuple
        Size of the figure.
    cmap : str
        Colormap to use for the heatmap.
    """
    plt.figure(figsize=figsize)
    sns.heatmap(corr_matrix, xticklabels=x_labels, yticklabels=y_labels, cmap=cmap, annot=False, fmt=".2f")
    plt.xlabel(f'Predictions from {species_pair[1]}')
    plt.ylabel(f'Predictions from {species_pair[0]}')
    plt.title('Cross-prediction Correlation')
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_dir, "crested", "figures",f'cross_species_corr_{"_".join(species_pair)}.pdf'), bbox_inches='tight')
    plt.close()

def compare_species_correlations(species_pair, reordered_classes, species_data, method='pearson'):
    """
    Compute and plot reordered correlation heatmap between two species' predictions.
    Parameters
    ----------
    species_pair : list of str
        Two species names, e.g. ['Mouse', 'Macaque'].
    reordered_classes : list of str
        Desired order of classes to plot.
    species_data : dict
        Dictionary with keys 'Mouse', 'Macaque', 'Human', each mapping to:
            {
                'cross_preds': ...,
                'adata': ...
            }
    method : str
        'pearson' or 'spearman'
    """
    cp_a = species_data[species_pair[0]]['cross_preds']
    cp_b = species_data[species_pair[1]]['cross_preds']
    adata_a = species_data[species_pair[0]]['adata']
    adata_b = species_data[species_pair[1]]['adata']
    ##
    corrs = compute_multi_correlation_matrices(
        cp_a, cp_b,
        adata_a, adata_b,
        species_key1=species_pair[0],
        species_key2=species_pair[1],
        method=method
    )
    #
    corr_a = corrs['species1']
    corr_b = corrs['species2']
    corr_combined = corrs['combined']  # 2D array
    ##
    obs_names_a = list(adata_a.obs_names)
    obs_names_b = list(adata_b.obs_names)
    ##
    index_order_a = [obs_names_a.index(cls) for cls in reordered_classes if cls in obs_names_a]
    index_order_b = [obs_names_b.index(cls) for cls in reordered_classes if cls in obs_names_b]
    ##
    reordered_obs_names_a = [obs_names_a[i] for i in index_order_a]
    reordered_obs_names_b = [obs_names_b[i] for i in index_order_b]
    # Slice all correlation matrices
    corr_a = corr_a[np.ix_(index_order_a, index_order_b)]
    corr_b = corr_b[np.ix_(index_order_b, index_order_a)]
    corr_combined = corr_combined[np.ix_(index_order_a, index_order_b)]
    # Plot
    plot_correlation_heatmap(corr_combined, reordered_obs_names_a, reordered_obs_names_b, species_pair)
    return corr_a, corr_b, corr_combined, reordered_obs_names_a, reordered_obs_names_b
