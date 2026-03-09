import anndata as ad
import pandas as pd
import numpy as np
import umap
import os
from tqdm import tqdm

from scipy.sparse import csr_matrix

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

##
data_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/cactus/analysis"

## ------------------------------------------------
## Handle ancestral alignments first
## ------------------------------------------------

## Load in ancestral CACTUS alignments
adata_files = {}
for species in ["human", "macaque", "marmoset"]:
    file_path = os.path.join(data_dir, f"{species}_zoonomia_overlap_HALPER_ancestral.h5ad")
    if os.path.exists(file_path):
        adata = ad.read_h5ad(file_path)
        adata.obs_names_make_unique()
        adata_files[species] = adata
    else:
        print(f"File for {species} not found at {file_path}")

## Merge all adatas
anc_adata = ad.concat(adata_files.values(), join='outer', label='source', keys=adata_files.keys(), fill_value=0)
anc_adata.obs_names_make_unique()

anc_adata.write_h5ad(os.path.join(data_dir, "all_species_zoonomia_overlap_HALPER_ancestral.h5ad"))

## ------------------------------------------------
## Now handle extant species alignments
## ------------------------------------------------

## Load in per-species CACTUS alignments to 447 species
adata_files = {}
for species in ["human", "macaque", "marmoset"]:
    file_path = os.path.join(data_dir, f"{species}_zoonomia_overlap_HALPER.h5ad")
    if os.path.exists(file_path):
        adata = ad.read_h5ad(file_path)
        adata.obs_names_make_unique()
        adata_files[species] = adata
        ## Ensure observations are unique
        ## adata_files[species] = adata_files[species][~adata_files[species].obs.index.duplicated(keep='first')]
    else:
        print(f"File for {species} not found at {file_path}")

## Merge all adatas
merged_adata = ad.concat(adata_files.values(), join='outer', label='source', keys=adata_files.keys(), fill_value=0)
merged_adata.obs_names_make_unique()

print(merged_adata.obs_names)

## Add basic stats
merged_adata.obs["n_species"] = (merged_adata.X > 0).sum(axis=1).A1
merged_adata.obs["sum_aligned_bp"] = merged_adata.X.sum(axis=1).A1
merged_adata.obs["mean_aligned_bp"] = merged_adata.obs["sum_aligned_bp"] / merged_adata.obs["n_species"]

## Create a UMAP
reducer = umap.UMAP(n_neighbors=15, min_dist=0.1)
embedding = reducer.fit_transform(merged_adata.X.toarray())
embedding_df = pd.DataFrame(embedding, columns=["UMAP1", "UMAP2"], index=merged_adata.obs_names)

## Add to anndata
merged_adata.obsm["X_umap"] = embedding_df.values

## --------------------------------------------------
## Determine levels of conservation
## --------------------------------------------------

## Indicator peaks with >= 90% alignment
N1_alignment_summary = csr_matrix((merged_adata.X.shape[0], merged_adata.X.shape[1]), dtype=bool)
N1_alignment_summary[((merged_adata.X/501) >= 0.9)] = True

## Count up for each peak the number of True values
merged_adata.obs["n_species_>=90%_aligned"] = N1_alignment_summary.sum(axis=1)
merged_adata.layers["N1_>=90%_aligned"] = N1_alignment_summary

## Indicator peaks with <= 10% alignment
N2_alignment_summary = csr_matrix((merged_adata.X.shape[0], merged_adata.X.shape[1]), dtype=bool)
N2_alignment_summary[((merged_adata.X/501) <= 0.1)] = True

## Count up for each peak the number of True values
merged_adata.obs["n_species_<=10%_aligned"] = N2_alignment_summary.sum(axis=1)
merged_adata.layers["N2_<=10%_aligned"] = N2_alignment_summary

## --------------------------------------------------
## Bring in annotations
## --------------------------------------------------
work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/"
anno_dir = os.path.join(work_dir, "analysis", "annotations")

## Ensure we are mapping the correct region per-species in the rare case of overlap between species
merged_adata.obs["region-id"] = merged_adata.obs_names + "-" + merged_adata.obs["species"].astype(str)

## Load annotations
for species in ["human", "macaque", "marmoset"]: 
    for anno_type in ["annotated_phyloP", "annotated_TE", "promoter"]:
        anno_path = os.path.join(anno_dir, f"{species}_peaks_{anno_type}.csv")
        print(f"Loading {anno_path}...")
        if os.path.exists(anno_path): #& f"{species}_{col}" not in merged_adata.obs.columns:
            anno_df = pd.read_csv(anno_path)
            anno_df["Name"] = anno_df["Chromosome"].astype(str) + ":" + anno_df["Start"].astype(str) + "-" + anno_df["End"].astype(str)
            if (species in ["marmoset", "macaque"]) & (anno_type == "annotated_phyloP"):
                ## Map chrom aliases back to original
                anno_df["Name"] = anno_df["ID"]
                anno_df["Name"] = anno_df["Name"].str.replace(r':\d+$', '', regex=True) ## Correction for allowing HALPER to set peak ids
            ##
            anno_df["Name"] = anno_df["Name"] + "-" + species
            anno_df.set_index("Name", inplace=True)
            ## Ensure index is unique
            anno_df = anno_df[~anno_df.index.duplicated(keep='first')]
            ## Subset to only those in merged_adata
            common_indices = merged_adata.obs.loc[
                merged_adata.obs["region-id"].isin(anno_df.index), "region-id"
            ].unique()
            anno_df = anno_df.loc[common_indices]
            print(anno_df.head())
            ## Subset to only those region IDs present in anno_df
            common_mask = merged_adata.obs["region-id"].isin(anno_df.index)
            common_ids = merged_adata.obs.loc[common_mask, "region-id"]
            ## Loop over annotation columns and map values by region-id
            for col in anno_df.columns:
                value_map = anno_df[col].to_dict()
                merged_adata.obs.loc[common_mask, f"{species}_{col}"] = common_ids.map(value_map)
        else:
            print(f"Annotation file for {species} not found at {anno_path}")

## Create a cross-species annotation for promoters
species = ["human", "macaque", "marmoset"]
fields = ["promoter", "phyloP_mean", "TE", "repName", "repClass", "repFamily"]

for field in fields:
    for sp in species:
        col_name = f"{sp}_{field}"
        if col_name in merged_adata.obs.columns:
            merged_adata.obs.loc[merged_adata.obs.species == sp, field] = merged_adata.obs.loc[merged_adata.obs.species == sp, col_name]
    merged_adata.obs[f"{field}"].value_counts(dropna=False)

## For repName, repClass and repFamily remove -1 from value with a ","
for field in ["repName", "repClass", "repFamily"]:
    merged_adata.obs[field] = merged_adata.obs[field].astype(str).str.replace("-1,", "")
    merged_adata.obs[field] = merged_adata.obs[field].astype(str).str.replace(",-1", "")

merged_adata.obs["promoter"].fillna(False, inplace=True)

## -----------------------------------------------------
## Compute species mixing (entropy) in UMAP space
## -----------------------------------------------------
# from sklearn.neighbors import NearestNeighbors
# from scipy.stats import entropy

# nbrs = NearestNeighbors(n_neighbors=30, algorithm='ball_tree').fit(merged_adata.obsm["X_umap"])
# distances, indices = nbrs.kneighbors(merged_adata.obsm["X_umap"])
# species_array = merged_adata.obs["species"].values
# entropy_values = []
# for idx in range(merged_adata.n_obs):
#     neighbor_species = species_array[indices[idx]]
#     species_counts = pd.Series(neighbor_species).value_counts(normalize=True)
#     ent = entropy(species_counts)
#     entropy_values.append(ent)
# merged_adata.obs["species_entropy"] = entropy_values

## -------------------------------------------------
## Add NCBI taxa information
## -------------------------------------------------
from ete3 import NCBITaxa

## Pull taxa from NCBI
ncbi = NCBITaxa()

## Add in taxanomic info for each species
species_list = list(adata.var.index)  # e.g., 'Homo_sapiens', 'Mus_musculus'
taxid_dict = {}

## Gather taxa info
for sp in tqdm(species_list):
    try:
        taxid = ncbi.get_name_translator([sp.replace("_", " ")])[sp.replace("_", " ")][0]
        lineage = ncbi.get_lineage(taxid)              # list of taxids up to root
        names = ncbi.get_taxid_translator(lineage)     # map taxid -> name
        ranks = ncbi.get_rank(lineage)                 # taxid -> rank
        taxid_dict[sp] = {ranks[t]: names[t] for t in lineage}
    except:
        taxid_dict[sp] = {}

## Convert to DataFrame
taxa_df = pd.DataFrame.from_dict(taxid_dict, orient='index')

## Remove no rank column
taxa_df = taxa_df.drop(columns=["no rank", "cellular root"])

## Add taxa info
merged_adata.var = pd.concat([merged_adata.var, taxa_df], axis=1)

## ----
## Add in Morgan's Primate calls
## ----
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/AnnoTables"
custom_anno = pd.read_csv(os.path.join(anno_dir, "primates_CACTUS_447_with_ncbi_taxa_custom_names.csv"), index_col=0)
custom_anno = custom_anno.loc[:,["primates_custom", "primates_shorthand"]]

merged_adata.var = merged_adata.var.merge(custom_anno, left_index=True, right_index=True, how="left")

## -------------------------------------------------
## Load in evo distance calculations
## -------------------------------------------------
evo_distance = pd.read_csv(os.path.join(data_dir, "zoonomia_overlap_HALPER_ancestral_species_evodist.csv"), index_col=0)

merged_adata.obs["evo_distance"] = evo_distance.loc[merged_adata.obs_names, "evo_distance"]

## --------------------------------------------------
## Include marmoset chrom alias
## --------------------------------------------------

## Marmoset chrom alias
chrom_alias = pd.read_csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/references/marmoset/ncbi/mcalja1.2.pat.x/genome/chromAlias_marmoset.tsv", sep="\t")
chrom_alias["name"] = "chr" + chrom_alias["name"].astype(str)

## --------------------------------------------------
## Load in gini scores
## --------------------------------------------------

gini_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/analysis/atac/specificity"

## Load in gini scores
gini_scores = pd.read_csv(os.path.join(gini_dir, "gini_scores_combined.csv"), index_col=0)

## Create region-id for mapping
gini_scores["region-id"] = (
    gini_scores.index.astype(str) + "-" + gini_scores["species"].astype(str)
)

## Update marmoset chrom names to NCBI names
marmoset_mask = gini_scores["species"] == "marmoset"
gini_scores_marmoset = gini_scores.loc[marmoset_mask, :].copy()
gini_scores_marmoset["Chromosome"] = gini_scores_marmoset.index.str.split(":").str[0]
gini_scores_marmoset["Start"] = gini_scores_marmoset.index.str.split(":").str[1].str.split("-").str[0].astype(int)
gini_scores_marmoset["End"] = gini_scores_marmoset.index.str.split(":").str[1].str.split("-").str[1].astype(int)

gini_scores_marmoset["Chromosome"] = gini_scores_marmoset["Chromosome"].map(
    dict(zip(chrom_alias["name"], chrom_alias["refseq"]))
)

## Update gini_scores with new chrom names
gini_scores_marmoset["new_index"] = (
    gini_scores_marmoset["Chromosome"].astype(str) + ":" +
    gini_scores_marmoset["Start"].astype(str) + "-" +
    gini_scores_marmoset["End"].astype(str) + "-marmoset"
)
gini_scores.loc[marmoset_mask, "region-id"] = gini_scores_marmoset["new_index"].values

## -----
## Now merge
## -----
merged_adata.obs["gini_scores"] = merged_adata.obs["region-id"].map(
    dict(zip(gini_scores["region-id"], gini_scores["gini_scores"]))
)
merged_adata.obs["gini_scores"].fillna(0, inplace=True)

merged_adata.obs.groupby("species")["gini_scores"].describe()

## -------------------------------------------------
## Load liftover annotations
## -------------------------------------------------
work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/BasalGanglia_paper_package/"
anno_dir = os.path.join(work_dir, "analysis", "annotations")

## Save for CACTUS analysis
liftover = pd.read_csv(os.path.join(anno_dir, "human_ref_liftover_HALPER_minMatch_0-5.tsv"), sep="\t", index_col=0)
liftover["region-id"] = liftover["human_ID"] + "-human"

## Map to merged_adata
merged_adata.obs = merged_adata.obs.merge(
    liftover,
    left_on = "region-id",
    right_on = "region-id",
    how="left"
)

## 
merged_adata.obs.index = merged_adata.obs['region-id'].str.replace("-human", "").str.replace("-macaque", "").str.replace("-marmoset", "").copy()

## --------------------------------------------------
## Save merged anndata
## --------------------------------------------------

## Remove region-id
merged_adata.obs = merged_adata.obs.drop(columns=["region-id"])

## Fix any obs column with mixed types by converting to string
for col in merged_adata.obs.columns:
    if merged_adata.obs[col].dtype == 'object':
        merged_adata.obs[col] = merged_adata.obs[col].astype(str)

##
merged_adata.write_h5ad(os.path.join(data_dir, "all_species_zoonomia_overlap_HALPER_all_anno.h5ad"))
