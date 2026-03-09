import os
from collections import defaultdict

import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

from Bio import Phylo

##
cactus_path = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/liftover/"
analysis_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/cactus/"

## Load species and ancestral .h5ad
adata_files = {
    "species": ad.read_h5ad(os.path.join(analysis_dir, "analysis", "all_species_zoonomia_overlap_HALPER.h5ad")),
    "ancestral": ad.read_h5ad(os.path.join(analysis_dir, "analysis", "all_species_zoonomia_overlap_HALPER_ancestral.h5ad"))
}

## Merge all adatas
merged_adata = ad.concat(
    adata_files.values(),
    axis=1,           
    join="outer", ## include all features (vars)
    merge="unique",   
)
del merged_adata.obs['source']
merged_adata.obs_names_make_unique()

## Compute alignment 90%
peak_length = 501  # Fixed length from HALPER
merged_adata.layers[">=90%_aligned"] = (merged_adata.X/peak_length >= 0.9).astype(int)

## --------------------------------------------
## Setup tree for evo distance walking
## --------------------------------------------
tree = Phylo.read(os.path.join(cactus_path, "zoonomia_447_no_length.nh"), "newick")  # or .newick

def tree_to_dataframe(tree):
    records = []
    ##
    def walk(clade, parent_name=None):
        node_name = clade.name if clade.name else f"internal_{len(records)}"
        is_leaf = len(clade.clades) == 0
        records.append({
            "node_name": node_name,
            "parent_name": parent_name,
            "branch_length": clade.branch_length,
            "is_leaf": is_leaf
        })
        for child in clade.clades:
            walk(child, node_name)
    ##
    walk(tree.root)
    return pd.DataFrame(records)

phylo_df = tree_to_dataframe(tree)
phylo_df["branch_length"] = phylo_df["branch_length"].fillna(0)

# ---------------------------------------------------------------------
# Visualization utilities
# ---------------------------------------------------------------------

# --------------------------------------------------------------
# Compute coordinates for each clade using depth (X) and order (Y)
# --------------------------------------------------------------
def get_clade_positions(tree):
    """
    Compute (x, y) coordinates for each clade in a rectangular layout.
    Returns a dict {clade: (x, y)}.
    """
    depths = tree.depths()
    terminals = tree.get_terminals()
    y_positions = {term: i for i, term in enumerate(terminals)}
    coords = {}
    ## Recursive function to calculate coordinates
    def calc_coords(clade):
        x = depths[clade]
        if clade.is_terminal():
            y = y_positions[clade]
        else:
            child_y = [calc_coords(c)[1] for c in clade.clades]
            y = np.mean(child_y)
        coords[clade] = (x, y)
        return x, y
    ##
    calc_coords(tree.root)
    return coords

from collections import deque
def get_polar_coords(tree):
    """
    Compute polar coordinates (angle θ, radius r) for each clade.
    θ is evenly spaced among leaves, r is normalized cumulative branch length.
    Returns a dict {clade: (θ, r)}.
    """
    depths = tree.depths()
    max_depth = max(depths.values()) or 1.0
    ##
    terminals = tree.get_terminals()
    n_leaves = len(terminals)
    leaf_angles = {term: 2 * np.pi * i / n_leaves for i, term in enumerate(terminals)}
    ##
    coords = {}
    ## Recursive function to assign polar coordinates
    def calc_coords(clade):
        r = depths[clade] / max_depth
        if clade.is_terminal():
            theta = leaf_angles[clade]
        else:
            child_thetas = [calc_coords(c)[0] for c in clade.clades]
            theta = np.mean(child_thetas)
        coords[clade] = (theta, r)
        return theta, r
    ##
    calc_coords(tree.root)
    return coords

def create_gradient_cmap(base_color, n_colors=100):
    """Create a light-to-base color gradient colormap."""
    colors = ["white", base_color]
    return LinearSegmentedColormap.from_list(f"{base_color}_gradient", colors, N=n_colors)

# ---------------------------------------------------------------------
# Visualization setup
# ---------------------------------------------------------------------

species_colors = {
    "human": "blue",
    "macaque": "green",
    "marmoset": "pink",
}

# ---------------------------------------------------------------------
# Plotting loop
# ---------------------------------------------------------------------

## Arguments
cmap_name = "plasma"
cmap_obj = plt.get_cmap(cmap_name)
norm = mpl.colors.Normalize(vmin=0, vmax=1)
width_scale = 1.0
color_by_species = False  # If True, use species-specific colormaps; else use shared cmap

##
species_alignment_df = {}

##
coord_dict = get_polar_coords(tree)


for species in ["human", "macaque", "marmoset"]:
    print(f"--- Rendering {species} ---")
    ## Extract alignment layer and compute per-node counts
    alignment = merged_adata[merged_adata.obs.species == species, :].layers[">=90%_aligned"].toarray()
    nodes = merged_adata.var_names
    n_aligned = alignment.sum(axis=0)
    ##
    node_counts = pd.DataFrame({
        "node_name": nodes,
        "n_peaks_aligned": n_aligned,
        "fraction_peaks": n_aligned / alignment.shape[0]
    })
    ##
    node_to_count = dict(zip(node_counts["node_name"], node_counts["n_peaks_aligned"]))
    ## Attach peak info to clades
    for clade in tree.find_clades():
        clade.n_peaks = node_to_count.get(clade.name, 0)
        clade.n_peaks_frac = clade.n_peaks / alignment.shape[0]
    ## Store for reference
    species_alignment_df[species] = pd.DataFrame({
        "node_name": [clade.name for clade in tree.find_clades()],
        "n_peaks": [clade.n_peaks for clade in tree.find_clades()],
        "n_peaks_frac": [clade.n_peaks_frac for clade in tree.find_clades()]
    })
    ## Polar plot setup
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10, 10))
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2)
    ax.axis("off")
    ## Draw branches
    for clade in tree.find_clades(order="preorder"):
        if clade.branch_length is None:
            continue
        ## Get coordinates
        theta, r = coord_dict[clade]
        path = tree.get_path(clade)
        parent = path[-2] if len(path) > 1 else tree.root
        ptheta, pr = coord_dict[parent]
        ## Draw line with width/color based on n_peaks_frac
        lw = max(0.1, clade.n_peaks_frac * width_scale)
        print(f"Clade: {clade.name}, n_peaks: {clade.n_peaks}, frac: {clade.n_peaks_frac:.3f}, lw: {lw:.2f}")
        if clade.n_peaks_frac <= 0.25:
            color = "lightgray"
            print(f"Note: {clade.name} has less than 0.25 aligned peaks.")
        else:
            if color_by_species:
                color = create_gradient_cmap(species_colors[species])(norm(clade.n_peaks_frac))
            else:
                color = cmap_obj(norm(clade.n_peaks_frac))
        ax.plot([ptheta, theta], [pr, r], color=color, lw=lw, solid_capstyle="round")
    ## Add colorbar
    if color_by_species:
        sm = mpl.cm.ScalarMappable(cmap=create_gradient_cmap(species_colors[species]), norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, fraction=0.03, pad=0.04)
        cbar.set_label("Fraction of aligned peaks (≥90%)", rotation=270, labelpad=15)
        cbar.ax.tick_params(labelsize=8)
        ## Save
        ax.set_title(f"{species.capitalize()} — Branch width ∝ # of peaks aligned ≥90%")
        plt.savefig(os.path.join(analysis_dir, "figures", f"phylo_tree_{species}_colors_branch_widths.pdf"))
    else:
        sm = mpl.cm.ScalarMappable(cmap=cmap_obj, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, fraction=0.03, pad=0.04)
        cbar.set_label("Normalized # of aligned peaks", rotation=270, labelpad=15)
        ## Save
        ax.set_title(f"{species.capitalize()} — Branch width ∝ # of peaks aligned ≥90%")
        plt.savefig(os.path.join(analysis_dir, "figures", f"phylo_tree_{species}_cmap_branch_widths.pdf"))
    ## Clean up
    plt.close(fig)

## Write out alignment summaries
for species, df in species_alignment_df.items():
    df.to_csv(os.path.join(analysis_dir, f"phylo_tree_{species}_alignment_summary.csv"), index=False)



# ## --- Attempt at matching itol style circular tree ---
# import toytree
# import matplotlib.pyplot as plt

# # Load your tree (same Newick file)
# ttree = toytree.tree(os.path.join(cactus_path, "zoonomia_447_no_length.nh"))

# # Ignore branch lengths (make all edges equal)
# # mod_tree = toytree.mod.edges_extend_tips_to_align(ttree).mod.edges_slider(prop=0.5)
# # Get the maximum root height (tree height)
# tree_data = ttree.get_node_data()
# tree_height = ttree.get_node_data().height.max()

# tree_data.height = tree_data.height / tree_height
# tree_data.dist = 1-tree_data.height

# # Scale all edge lengths so that the total height = 1
# unit_tree = ttree.set_node_data("height", dict(zip(tree_data.name, tree_data.height)))
# unit_tree = unit_tree.set_node_data("dist", dict(zip(tree_data.name, tree_data.dist)))

# # Draw iTOL-style circular tree
# canvas, axes, mark = unit_tree.draw(
#     width=900,
#     height=900,
#     layout="c",          # circular layout
#     node_sizes=0,        # hide node dots
#     edge_type="c",       # 'a' = arcs instead of straight lines
#     node_labels=False,
#     tip_labels=False,
#     #tip_labels_align=True,
#     #use_edge_lengths=True,
#     edge_widths=0.8,
# )
# toytree.save(canvas, os.path.join(analysis_dir, "figures", "phylo_tree_circular_itol_style.pdf"))
