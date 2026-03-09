import os
import anndata as ad
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib
import pertpy as pt
import scanpy as sc

from tqdm import tqdm

## -------------------------
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/xspecies/BasalGanglia/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis"
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/anno_tables/"

## -------------------------
cluster_meta = pd.read_excel(os.path.join(anno_dir, "HMBA_BG_consensus_annotation.xlsx"), sheet_name="consensus_anno_pre-print")
cluster_meta = cluster_meta.drop_duplicates(subset=["Group"]).reset_index(drop=True)

## Create color map
color_map = dict(zip(cluster_meta.Group, cluster_meta.color_hex_group))

## -------------------------
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/xspecies/BasalGanglia"
adata = ad.read_h5ad(os.path.join(data_dir, "HMBA_Human_Macaque_Marmoset_Siletti_BG_alignment_Neuron.h5ad"))

## -------------------------
## propagate consensus annotations to all cells

from sklearn.neighbors import KNeighborsClassifier
from collections import Counter
from sklearn.neighbors import KNeighborsClassifier
def distribute_evenly(a, k):
    n = len(a)  # Number of elements
    r = [0] * n  # Initialize result vector with zeros
    # First pass: Try to distribute evenly based on the initial total
    avg = k // n
    for i in range(n):
        r[i] = min(a[i], avg)
    remaining = k - sum(r)  # Calculate remaining after initial distribution
    # Distribute the remaining as evenly as possible
    while remaining > 0:
        # Calculate the updated average for the remaining amount
        updates = 0
        for i in range(n):
            if r[i] < a[i]:
                potential_add = min(a[i] - r[i], remaining // (n - updates))
                if potential_add > 0:
                    r[i] += potential_add
                    remaining -= potential_add
                    updates += 1
                    if remaining == 0:
                        break
        # In case there are no updates, but still remaining (to avoid infinite loop)
        if updates == 0:
            break
    return r

def select_k_cells(adata, obs_column, k):
    a = adata.obs[obs_column].value_counts()
    categories = a.index
    r = distribute_evenly(a,k)
    inds=[]
    for i,cat in enumerate(categories):
        cat_cells = adata.obs.index[adata.obs[obs_column] == cat].tolist()
        inds.append(np.random.choice(cat_cells, size=r[i], replace=False))
    selected_indices=[item for sublist in inds for item in sublist]
    return adata[selected_indices,:]

def consensus_knn_transfer(adata, label_col, train_mask, test_mask, n_neighbors=20, n_iter=10):
    preds_list = []
    adata.obs[label_col + '_propagated'] = list(adata.obs[label_col])
    train_data = adata[train_mask]
    # Determine number of cells to sample from each class
    min_count = train_data.obs[label_col].value_counts().min()
    k_val = len(train_data.obs[label_col].unique()) * min_count * 2
    for _ in tqdm(range(n_iter)):
        # Downsample training data evenly
        equalized_train = select_k_cells(train_data, label_col, k_val)
        neigh = KNeighborsClassifier(n_neighbors=n_neighbors)
        X_train = equalized_train.obsm['X_scVI']
        y_train = equalized_train.obs[label_col]
        neigh.fit(X_train, y_train)
        test_data = adata[test_mask]
        X_test = test_data.obsm['X_scVI']
        preds = neigh.predict(X_test)
        preds_list.append(preds)
    # Stack predictions from all iterations and take majority vote for each test cell.
    preds_array = np.vstack(preds_list)  # shape: (n_iter, n_test)
    consensus_preds = [
        Counter(preds_array[:, j]).most_common(1)[0][0]
        for j in range(preds_array.shape[1])
    ]
    adata.obs.loc[test_mask, label_col + '_propagated'] = consensus_preds
    return adata


adata.obs['Group'] = adata.obs['Group'].astype(str)

train_mask = adata.obs['Group']!='nan'
test_mask  = adata.obs['Group']=='nan'
adata = consensus_knn_transfer(adata, 'Group', train_mask, test_mask,
                               n_neighbors=30, n_iter=20)

## -------------------------
adata.obs["study_id"] = adata.obs["organism"].map({
    "Human": "HMBA",
    "Macaque": "HMBA",
    "Marmoset": "HMBA",
    "Homo sapiens": "Siletti"
})

## -------------------------
## Plot UMAP with consensus annotation colors
colors = adata.obs["Group_propagated"].map(color_map).values.astype(str)
colors[pd.isna(colors)] = "lightgrey"  # Handle NaN values
colors[adata.obs["study_id"] != "Siletti"] = "#add8e6"  # Highlight mouse cells only

## Plot UMAP colored by cell type
plt.figure(figsize=(12,12))
adata_plot = adata[adata.obs["study_id"] != "Siletti", :]
plt.scatter(adata_plot.obsm["X_umap"][:,0], 
            adata_plot.obsm["X_umap"][:,1], 
            c="#d3d3d3",
            s=0.1, 
            alpha=0.2)
adata_plot = adata[adata.obs["study_id"] == "Siletti", :]
plt.scatter(adata_plot.obsm["X_umap"][:,0], 
            adata_plot.obsm["X_umap"][:,1], 
            c=adata_plot.obs["Group_propagated"].map(color_map).values.astype(str),
            s=0.1, 
            alpha=0.8)
## Add color legend
plt.title("UMAP of HMBA BG Neurons - Consensus Annotation", fontsize=16)
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.savefig(os.path.join(analysis_dir, "homology", "figures", f"xpsecies_siletti_neuron_v2.png"), dpi=900)

## -------------------------
## Milo analysis

## Initialize object for Milo analysis
milo = pt.tl.Milo()
mdata = milo.load(adata)

## Milo setup
sc.pp.neighbors(mdata["rna"], use_rep="X_scVI", n_neighbors=150)
milo.make_nhoods(mdata["rna"], prop=0.1)

# mdata["rna"].obsm["nhoods"]
# mdata["rna"][mdata["rna"].obs["nhood_ixs_refined"] != 0].obs[["nhood_ixs_refined", "nhood_kth_distance"]]

## Neighborhood size distribution
nhood_size = np.array(mdata["rna"].obsm["nhoods"].sum(0)).ravel()

plt.figure(figsize=(7,5))
plt.hist(nhood_size, bins=100)
plt.xlabel("# cells in nhood")
plt.ylabel("# nhoods")
plt.title("Neighborhood size distribution")
plt.tight_layout()
plt.savefig(os.path.join(analysis_dir, "homology", "figures", f"xpsecies_mouse_neuron_milo_nhood_size.png"), dpi=900)

adata.obs["load_milo"] = adata.obs["load_id"]
adata.obs["load_milo"] = adata.obs["load_milo"].astype(str)
adata.obs["load_milo"][adata.obs["load_milo"] == 'nan'] = adata.obs["sample_id"][adata.obs["load_milo"] == 'nan']

## Number of cells from each sample counted in a neighbourhood
mdata = milo.count_nhoods(mdata, sample_col="load_id")
mdata["milo"]

## Reorder categories (by default, the last category is taken as the condition of interest)
mdata["rna"].obs["study_id"] = mdata["rna"].obs["study_id"].astype('category').cat.reorder_categories(["HMBA", "Siletti"])
milo.da_nhoods(mdata, design="~study_id", solver="pydeseq2", model_contrasts="study_idHMBA-study_idSiletti")

mdata["milo"].obs
mdata["milo"].var

## Basic plots of DA results
old_figsize = plt.rcParams["figure.figsize"]
plt.rcParams["figure.figsize"] = [10, 5]
plt.subplot(1, 2, 1)
plt.hist(mdata["milo"].var.PValue, bins=50)
plt.xlabel("P-Vals")
plt.subplot(1, 2, 2)
plt.plot(mdata["milo"].var.logFC, -np.log10(mdata["milo"].var.SpatialFDR), ".")
plt.xlabel("log-Fold Change")
plt.ylabel("- log10(Spatial FDR)")
plt.tight_layout()
plt.rcParams["figure.figsize"] = old_figsize
plt.savefig(os.path.join(analysis_dir, "homology", "figures", f"xpsecies_human_neuron_milo_da_results.png"), dpi=900)

## -------------------------
## Visualize on embedding
## --- 

## Build neighborhood graph structed on UMAP
milo.build_nhood_graph(mdata)

plt.rcParams["figure.figsize"] = [12, 12]
milo.plot_nhood_graph(
    mdata,
    alpha=0.1,  # SpatialFDR level (1%)
    min_size=0.1,  # Size of smallest dot
)
plt.savefig(os.path.join(analysis_dir, "homology", "figures", f"xpsecies_mouse_neuron_milo_nhood_graph_contrast.png"), dpi=900)

## Annotate neighborhoods based on Group_propagated
milo.annotate_nhoods(mdata, anno_col="Group_propagated")

## Label mixed neighborhoods as "Mixed"
mdata["milo"].var["nhood_annotation"] = mdata["milo"].var["nhood_annotation"].cat.add_categories("Mixed")
mdata["milo"].var.loc[mdata["milo"].var["nhood_annotation_frac"] < 0.6, "nhood_annotation"] = "Mixed"

milo.plot_da_beeswarm(mdata, 
                      alpha=0.1)
plt.savefig(os.path.join(analysis_dir, "homology", "figures", f"xpsecies_mouse_neuron_milo_da_beeswarm.png"), dpi=900)

sc.pl.violin(mdata["milo"].T, "logFC", groupby="nhood_annotation", rotation=90, orient="h", show=False)

## Define ordered categories for violin plot
sorted_annos = cluster_meta.Group.tolist() + ["Mixed"]
sorted_annos = [anno for anno in sorted_annos if anno in mdata["milo"].T.obs["nhood_annotation"].cat.categories]

mdata["milo"].T.obs["nhood_annotation"] = mdata["milo"].T.obs["nhood_annotation"].cat.reorder_categories(sorted_annos, ordered=True)

## Color per observation based on annotation
colors = [color_map[anno] if anno in color_map else "lightgrey" for anno in sorted_annos]

# --- Build color map per annotation ---
palette = {anno: color_map.get(anno, "lightgrey") for anno in sorted_annos}

# --- Plot ---
plt.figure(figsize=(8, 12))
sns.violinplot(
    data=mdata["milo"].T.obs,
    x="logFC",
    y="nhood_annotation",
    hue="nhood_annotation",
    palette=palette,
    orient="h",
    order=sorted_annos,
    inner="quart",
    density_norm="width"
)
plt.axvline(x=0, color="black", linestyle="--")
plt.legend([], [], frameon=False)  # optional: hide redundant hue legend
plt.tight_layout()
plt.savefig(
    os.path.join(analysis_dir, "homology", "figures", "xpsecies_mouse_neuron_milo_da_violin.pdf"),
    dpi=900
)

## -------------------------
anno_col = "nhood_annotation"
nhood_adata = mdata["milo"].T.copy()
alpha = 0.1

##
sorted_annos = (
            nhood_adata.obs[[anno_col, "logFC"]].groupby(anno_col).median().sort_values("logFC", ascending=True).index
        )

##
anno_df = nhood_adata.obs[[anno_col, "logFC", "SpatialFDR"]].copy()
anno_df["is_signif"] = anno_df["SpatialFDR"] < alpha
anno_df = anno_df[anno_df[anno_col] != "nan"]

##
obs_col = nhood_adata.uns["annotation_obs"]
if palette is None:
    palette = dict(
        zip(
            mdata[feature_key].obs[obs_col].cat.categories,
            mdata[feature_key].uns[f"{obs_col}_colors"],
            strict=False,
        )
    )
sns.violinplot(
    data=anno_df,
    y=anno_col,
    x="logFC",
    order=sorted_annos,
    inner=None,
    orient="h",
    palette=palette,
    linewidth=0,
    scale="width",
)
sns.stripplot(
    data=anno_df,
    y=anno_col,
    x="logFC",
    order=sorted_annos,
    size=2,
    hue="is_signif",
    palette=["grey", "black"],
    orient="h",
    alpha=0.5,
)
plt.legend(loc="upper left", title=f"< {int(alpha * 100)}% SpatialFDR", bbox_to_anchor=(1, 1), frameon=False)
plt.axvline(x=0, ymin=0, ymax=1, color="black", linestyle="--")
plt.savefig(os.path.join(analysis_dir, "homology", "figures", f"xpsecies_mouse_neuron_milo_da_violin.png"), dpi=900)
