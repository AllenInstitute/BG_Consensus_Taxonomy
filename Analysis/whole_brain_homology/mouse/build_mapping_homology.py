import anndata as ad
import pandas as pd
import numpy as np
import sciduck
import h5py

##
def grep(l, s):
    return [i for i in l if s in i]

species = "xspecies"

## Helpful locations which are assumed to already exist
anno_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/AnnoTables"
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia"
cxgdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/xspecies"

## ----------------------------------------
with h5py.File(cxgdir + "/HMBA_Human_Macaque_Marmoset_Mouse_snRNA-seq_BG_alignment_v5.h5ad") as f:
    adata_obs = ad.experimental.read_elem(f['obs'])
    adata_obsm = ad.experimental.read_elem(f['obsm'])

## Add in mapping results from MapMyCells to primates
ait_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/"
label_df = pd.DataFrame()
for species in ["Human", "Macaque", "Marmoset"]:
    print(species)
    if species in ["Human", "Macaque"]:
        with h5py.File(f"{ait_dir}/{species}/BasalGanglia/{species}_basalganglia_AIBS_BICAN_MapMyCells.h5ad") as f:
            species_obs = ad.experimental.read_elem(f['obs'])
    elif species == "Marmoset":
        with h5py.File(f"{ait_dir}/{species}/BasalGanglia/{species}_basalganglia_anno_latest_MapMyCells.h5ad") as f:
            species_obs = ad.experimental.read_elem(f['obs'])
        species_obs = species_obs.loc[:, ~species_obs.columns.str.contains(r'Flat|Siletti_|Siletti-')]
    filtered_columns = [col for col in species_obs.columns if col.startswith(('Siletti', 'ABCmouse', 'HMBA_WB'))]
    label_df = pd.concat([label_df, species_obs.loc[species_obs.index.isin(adata_obs.index), filtered_columns]], axis=0)
        
## Now merge adata_obs and label_df based on index
merged_obs = adata_obs.merge(label_df, left_index=True, right_index=True, how='left')
primate_obs = merged_obs.loc[merged_obs.organism.isin(["Human", "Macaque", "Marmoset"]), :]

## Taxonomy
consensus_taxonomy = pd.read_csv(os.path.join(anno_dir, "consensus_anno_pre-print.csv"), index_col=0)
consensus_taxonomy = consensus_taxonomy.loc[:,["Group", "accession_group", "Subclass", "accession_subclass", "Class", "accession_class", "Neighborhood", "accession_neighborhood"]]

## Lets figure out the best mapped labels for each level of the consenus taxonomy
groupby = "organism"
min_proportion = 0.1
top_n = 3

##
mapped_labels = {
    'ABCmouse_MapMyCells': {
        "Class": "CLAS",
        "Subclass": "SUBC",
        "Group": "SUPT"
    },
    'HMBA_WB_Human_MapMyCells': {
        "Class": "Class",
        "Subclass": "Subclass",
    },
    'Silettisub_MapMyCells': {
        "Class": "Supercluster",
        "Subclass": "Cluster",
        "Group": "Subcluster"
    },
}

## Intersect non-empty sets across species
def safe_intersection(row):
    sets = [s for s in row if isinstance(s, set) and len(s) > 0]
    return set.intersection(*sets) if sets else set()

##
for mapping in mapped_labels:
    for anno in ["Class", "Subclass", "Group"]:
        if anno not in mapped_labels[mapping]:
            continue
        mapped_label = f"{mapping}_{mapped_labels[mapping][anno]}_label"
        ## Filter out low mapping cells
        correlation_low = primate_obs[f"{mapping}_{mapped_labels[mapping][anno]}_avg_correlation"].quantile(0.25)
        filtered_obs = primate_obs.loc[(primate_obs[f"{mapping}_{mapped_labels[mapping][anno]}_avg_correlation"] >= correlation_low) | 
                                        (primate_obs[f"{mapping}_{mapped_labels[mapping][anno]}_bootstrapping_probability"] >= 0.6),:].copy()
        ## Group by cell_type, species, and mapped_label to count occurrences
        grouped = filtered_obs.groupby([anno, groupby, mapped_label]).size().reset_index(name='count')
        grouped['total'] = grouped.groupby([anno, groupby])['count'].transform('sum')
        grouped['proportion'] = grouped['count'] / grouped['total']
        grouped = grouped[grouped['proportion'] >= min_proportion]
        ## Get top_n mapped_labels per (cell_type, species)
        top_hits = (
            grouped.sort_values([anno, "organism", 'proportion'], ascending=[True, True, False])
            .groupby([anno, "organism"])
            .head(top_n)
        )
        top_hits['label_info'] = top_hits[mapped_label]
        top_hits['label_score'] = top_hits['proportion'].apply(lambda x: f"{x:.2f}") 
        ## Aggregate into a single string per cell_type and species
        summary = (
            top_hits.groupby([anno, "organism"])['label_info']
            .apply(lambda x: " | ".join(x))
            .unstack(fill_value='')
        )
        ## Group to get sets of mapped_label per cell_type and species
        label_sets = (
            top_hits.groupby([anno, groupby])[mapped_label]
            .apply(set)
            .unstack(fill_value=set())
        )
        ## Replace NaNs in label_sets with empty sets
        label_sets = label_sets.applymap(lambda x: x if isinstance(x, set) else set())
        ## Get the intersection of sets across species for each cell_type
        common_hits = label_sets.apply(safe_intersection, axis=1)
        summary[f'reciprocal_hits'] = common_hits.apply(lambda x: " | ".join(sorted(x)) if x else '')
        summary.columns = [col + f'_{mapped_label}' for col in summary.columns]
        ##
        consensus_taxonomy = consensus_taxonomy.merge(summary, left_on=anno, right_index=True, how='left')

##
consensus_taxonomy.to_csv(os.path.join(anno_dir, f"consensus_taxonomy_mapping_homologies_Siletti.csv"))