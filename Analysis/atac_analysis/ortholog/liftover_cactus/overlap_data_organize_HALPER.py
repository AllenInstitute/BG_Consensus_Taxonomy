import pandas as pd
import numpy as np
import pyranges as pr
from tqdm import tqdm
import umap
import matplotlib.pyplot as plt
import gzip
import glob 
import os
import anndata as ad
from scipy.sparse import csr_matrix

## alias
chrom_alias = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/references/marmoset/ncbi/mcalja1.2.pat.x/genome/chromAlias_marmoset.tsv", sep="\t")
chrom_alias["name"] = "chr" + chrom_alias["name"].astype(str)

##
cactus_path = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/cactus/analysis"

## 
species_peak_files = {
    species: f"/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/data/{species}/peaks/merged_peaks.bed"
    for species in ["marmoset"]
}

path_to_species_overlaps = {
    "human": "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/cactus/hal/human",
    "macaque": "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/cactus/hal/macaque",
    "marmoset": "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/cactus/hal/marmoset"
}

species_conversion = {
    "human": "Homo_sapiens",
    "macaque": "Macaca_mulatta",
    "marmoset": "Callithrix_jacchus",
    "mouse": "Mus_musculus"
}

## Narrow peak file
names = ["Chromosome", "Start", "End", "summit_pos", "ID", "target_length", "query_length", "target_summit_to_target_ortholog_start_length", "target_summit_to_target_ortholog_end_length"]

##
for species, species_peak_file in species_peak_files.items():
    ## Create save file path
    save_file = os.path.join(cactus_path, f"{species}_zoonomia_overlap_HALPER.h5ad")
    # if os.path.exists(save_file):
    #     print(f"File {save_file} already exists. Skipping...")
    #     continue
    ## Read in species peaks
    species_bed = pd.read_csv(species_peak_file, sep="\t", header=None)
    species_bed.columns = ["chrom", "start", "end"]
    species_bed["peak_length"] = species_bed["end"] - species_bed["start"] - 1
    if species == "marmoset":
        species_bed["orig_id"] = species_bed["chrom"] + ":" + species_bed["start"].astype(str) + "-" + species_bed["end"].astype(str)
        species_bed["chrom"] = species_bed["chrom"].map(dict(zip(chrom_alias["name"], chrom_alias["refseq"])))
        species_bed["chrom"] = species_bed["chrom"].fillna(species_bed["orig_id"].str.split(":").str[0])
        species_bed["chrom"].replace("v1_random", ".1", regex=True, inplace=True)
    ## Set ID
    species_bed["id"] = species_bed["chrom"] + ":" + species_bed["start"].astype(str) + "-" + species_bed["end"].astype(str)
    species_bed.set_index("id", inplace=True)
    ## Get all overlap files for this species
    speciesTo_overlapFiles = [f for f in  glob.glob(f"{path_to_species_overlaps[species]}/*HALPER_.narrowPeak.gz") if "Anc" not in f]
    ## Ignore ancestor nodes with "Anc" followed by numbers in name
    print(f"Found {len(speciesTo_overlapFiles)} overlap files for {species}.")
    ## Read in hits
    overlapDf = pd.DataFrame(index=species_bed.index.tolist())
    for f in tqdm(speciesTo_overlapFiles):
        if species == "marmoset":
            species_cactus = os.path.basename(f).replace(f"merged_peaks_updated_chrom.", "")
        else:
            species_cactus = os.path.basename(f).replace(f"merged_peaks.", "")
        species_cactus = species_cactus.replace(".HALPER_.narrowPeak.gz","")
        species_cactus = species_cactus.replace(f"{species_conversion[species]}To","")
        df = pd.read_csv(f, sep="\t", names=names, header=None, compression="gzip")
        df["ID"] = df["ID"].str.replace(r':\d+$', '', regex=True) ## Correction for allowing HALPER to set peak ids
        df.set_index("ID", inplace=True)
        ## Ensure order matches 
        species_bed_df = species_bed
        species_bed_df["overlap"] = 0
        species_bed_df.loc[species_bed_df.index.intersection(df.index), "overlap"] = df["target_length"].values
        if np.all(species_bed.index == species_bed_df.index):
            overlapDf[species_cactus] = species_bed_df.loc[:,"overlap"].values ## df holds other metrics FYI
        else:
            print(f"Issue: {species_cactus}")
    ## Create AnnData object
    overlap_adata = ad.AnnData(X = csr_matrix(overlapDf.values),
                                obs = pd.DataFrame(index=overlapDf.index.tolist()),
                                var = pd.DataFrame(index=overlapDf.columns.tolist()))
    overlap_adata.obs["species"] = species_bed["species"].values if "species" in species_bed.columns else species
    ## Save some time and store the original peak coordinates as well
    if species == "marmoset":
        overlap_adata.obs["orig_id"] = species_bed["orig_id"].values
    ##
    overlap_adata.write_h5ad(save_file)



## --------------------------------------------------
## Plot phylogenetic tree of species in cactus HAL
## --------------------------------------------------
# from Bio import Phylo
# import matplotlib.pyplot as plt

# tree = Phylo.read(os.path.join(cactus_path, "447-mammalian-2022v1.nh"), "newick")  # or .newick
# tree.ladderize()  # optional tidy ordering

# ## Get all species (leaf names) from your adata
# species_set = set(overlap_adata.var.index)

# ## Dictionary to store leaves under each internal node
# node_dict = {}

# ## Helper to traverse the tree
# def get_leaves(node):
#     """Return all leaf names under this node."""
#     if node.is_terminal():
#         return [node.name]
#     leaves = []
#     for clade in node.clades:
#         leaves.extend(get_leaves(clade))
#     return leaves

# ## Iterate over all clades (nodes)
# for i, clade in enumerate(tree.find_clades(order="preorder")):
#     if not clade.is_terminal():
#         leaves = get_leaves(clade)
#         leaves_in_var = [s for s in leaves if s in species_set]
#         node_name = clade.name or f"Node_{i}"
#         node_dict[node_name] = leaves_in_var

# ## Convert to DataFrame: columns = internal nodes, rows = species, values = 0/1
# annotation_df = pd.DataFrame(False, index=overlap_adata.var.index, columns=node_dict.keys())
# for node, leaves in node_dict.items():
#     annotation_df.loc[leaves, node] = True

# ## Add to adata.var
# overlap_adata.var = pd.concat([overlap_adata.var, annotation_df], axis=1)

# ## Plot tree
# plt.figure(figsize=(16,10), dpi=300)
# ax = plt.gca()

# # Draw in circular layout
# Phylo.draw(tree, do_show=False, axes=ax, branch_labels=None, label_func=lambda x: x.name, circular=True)

# # Adjust label font size
# for text in ax.texts:
#     text.set_fontsize(2)

# # Adjust line thickness
# for line in ax.get_lines():
#     line.set_linewidth(0.25)

# ##
# plt.tight_layout()
# plt.savefig("/home/nelson.johansen/human_peaks_phylo_tree.png", dpi=300)

# overlap_adata.write_h5ad(os.path.join(path_to_species_overlaps, "human_cactus_alignment_basepair_overlap_capped.h5ad"))


