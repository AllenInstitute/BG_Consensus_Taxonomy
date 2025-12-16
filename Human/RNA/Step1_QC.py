import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scvi
import os
import sys
import sciduck as sd
from datetime import datetime


##
species = "Human"
pipeline = "AIBS"

## Helpful locations which are assumed to already exist
data_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}"
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia"
figure_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia/figures"
model_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia/models"
cxg_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/HMBA/{species}"

os.chdir(work_dir)



## Load expression data
adata = ad.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_BICAN.h5ad"))

adata.X = adata.X.astype(np.int32) ## Data is float64 from NeMO which is excessive for counts.

for col in adata.obs.columns:
    if adata.obs[col].dtype == 'int64':
        adata.obs[col] = adata.obs[col].astype('int32')

for col in adata.var.columns:
    if adata.var[col].dtype == 'int64':
        adata.var[col] = adata.var[col].astype('int32')
        

adata.obs['ROI'] = adata.obs['roi'].copy()
adata.obs['ROI'] = adata.obs['ROI'].apply(lambda x: f'Human {x}' if not x.startswith('Human') else x)


adata.obs['merged_ROIs'] = adata.obs['ROI'].copy()
adata.obs['merged_ROIs'].replace(['Human CaB', 'Human CaH', 'Human CaT'],'Human Ca', inplace=True)
adata.obs['merged_ROIs'].replace(['Human GPeR', 'Human GPeC'],'Human GPe', inplace=True)
adata.obs['merged_ROIs'].replace(['Human NACc', 'Human NACs'],'Human NAC', inplace=True)
adata.obs['merged_ROIs'].replace(['Human PuC','Human PuPV', 'Human PuR', ],'Human Pu', inplace=True)
adata.obs['merged_ROIs'].replace(['Human SN'],'Human SN-VTA', inplace=True)


## Normalize the count matrix after storing the raw counts in the raw slot
adata.raw = adata.copy()

## Calculate some basic QC metrics around UMI and gene counts
sc.pp.calculate_qc_metrics(adata, inplace=True)
adata.var['mt'] = adata.var.index.str.startswith('MT-') 
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

##
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)


## Add manual QC metrics
adata.obs["total_genes"] = (adata.X > 0).sum(axis=1)

## Plot QC metrics
sc.pl.violin(
    adata,
    ["total_genes", "total_counts", "doublet_score","pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    save=f"_{species}_basic_qc.png",
)


## Perform cell level basic QC
sd.qc.add_range_constraint(adata, column="total_counts", gt=2e3, lt=1e5)
sd.qc.add_range_constraint(adata, column="total_genes", gt=1000, lt=13000)
sd.qc.add_range_constraint(adata, column="doublet_score", lt=0.35)
sd.qc.add_range_constraint(adata, column="pct_counts_mt", lt=3)


## ---- FIRST SHAPE CHANGE ----
## Remove flagged cells
adata = sd.qc.apply_constraints(adata, inplace=True)
adata.obs.keeper_cells.value_counts()

#adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_basic_qc.h5ad"))

## -----------------------------
## Run scVI to remove donor effects
#adata = sc.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_basic_qc.h5ad"))

##
adata.layers["UMIs"] = adata.raw.X ## This is a shallow copy, no extra space is used!
sc.pp.highly_variable_genes(adata, n_top_genes=4000, layer="UMIs", subset=False, flavor="seurat_v3", batch_key="donor_id")

##
adata_hvg = adata[:, adata.var.highly_variable].copy()

#import torch
#torch.set_float32_matmul_precision("medium")

#from torch.utils.data import DataLoader
#train_dataloader = DataLoader(adata_hvg, batch_size=32, shuffle=True, num_workers=7)


## Run scVI with known confounders
scvi.model.SCVI.setup_anndata(
    adata_hvg,
    layer="UMIs",
    batch_key="donor_id", ## This holds species for Macaque (rhesus, nemestrina). Human use donor_id here?
    categorical_covariate_keys=["donor_id"],
)



model = scvi.model.SCVI(adata_hvg, 
                        dispersion="gene-batch", 
                        n_hidden=256, 
                        n_latent=64, 
                        n_layers=3)
model.train(max_epochs=200)
 
 
## Save the latent space from scVI for downstream analysis
adata.obsm["X_scVI"] = model.get_latent_representation()

del adata_hvg
del adata.layers["UMIs"]


## UMAP from scVI latent space
sc.pp.neighbors(adata, use_rep = "X_scVI")
adata = sc.tl.umap(adata, min_dist=0.3, copy=True)


## -----------------------------
## Run leiden for basic cluster-wise qc
#sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.leiden(adata, resolution=20, key_added="leiden_scVI") #10

## Save data to explore in cellxgene 
#adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_scVI_cluster.h5ad")) ## Overwrites previous save.
#adata.write_h5ad(os.path.join(cxg_dir, f"{species}_basalganglia_{pipeline}_HMBA_scVI_cluster.h5ad"))

# transfer AIT19.3 lables
AIT193 = sc.read_h5ad("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/Human/allHuman_preHMBA_HMBA1stOne/0.RawData/1.RNA/merge/analysisv0.5/Merge3ClassafterGeneQCperClass/AIT19.3/AIT19.3_4khvg_updated.h5ad")
AIT193.obs['cluster_AIT19.3'] = AIT193.obs['cluster']
AIT193.obs['Neighborhood_AIT19.3'] = AIT193.obs['Neighborhood']
AIT193.obs['Class_AIT19.3'] = AIT193.obs['Class']
AIT193.obs['Subclass_AIT19.3'] = AIT193.obs['Subclass']
AIT193.obs['Group_AIT19.3'] = AIT193.obs['Group']

AIT193.obs["new_cellindex"] = AIT193.obs["barcodes"].astype(str) + "-" + AIT193.obs["load_name"].astype(str)
adata.obs["new_cellindex"] = adata.obs_names


columns_to_transfer =['cluster_AIT19.3', 'Neighborhood_AIT19.3', 'Class_AIT19.3','Subclass_AIT19.3','Group_AIT19.3', 'new_cellindex']
df_transfer = AIT193.obs[columns_to_transfer].copy()
del AIT193

adata.obs = adata.obs.merge(df_transfer, on="new_cellindex", how="left")
adata.obs[columns_to_transfer] = adata.obs[columns_to_transfer].astype("category")
#adata.obs['Neighborhood_AIT19.3'].value_counts()
#adata
adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_scVI_LeidenCluster_tags.h5ad")) 

cols_to_keep_in_h5ad = [ 'umi.counts', 'library_prep',  'gene_counts_0', 'doublet_score', 'exclude', 'exclude2', 'genome', 'pipeline_version', 
                         'gex_fraction_of_transcriptomic_reads_in_cells', 'gex_mean_raw_reads_per_cell', 'gex_median_umi_counts_per_cell', 'gex_median_genes_per_cell', 'gex_percent_duplicates', 
                        'gex_q30_bases_in_umi', 'gex_q30_bases_in_barcode', 'gex_q30_bases_in_read_2', 'gex_q30_bases_in_sample_index_i1', 'gex_q30_bases_in_sample_index_i2', 
                        'gex_reads_mapped_antisense_to_gene', 'gex_reads_mapped_confidently_to_exonic_regions', 'gex_reads_mapped_confidently_to_genome', 'gex_reads_mapped_confidently_to_intergenic_regions',
                        'gex_reads_mapped_confidently_to_intronic_regions', 'gex_reads_mapped_confidently_to_transcriptome', 'gex_reads_mapped_to_genome', 'gex_reads_with_tso', 'gex_sequenced_read_pairs', 
                        'gex_total_genes_detected', 'keeper_median_genes', 'keeper_cells', 'percent_keeper', 'percent_doublet', 'percent_usable', 'tso_percent', 'pipeline', 'pass_fail',  'batch', 'batch_vendor_name',
                        'tube', 'tube_internal_name', 'tube_contents_nm', 'tube_contents_nm_from_vendor', 'tube_avg_size_bp', 'tube_input_fmol', 'facs_container', 'sample_name', 
                        'patched_cell_container', 'studies', 'hemisphere_name', 'sample_quantity_count', 'sample_quantity_pg', 'donor_name', 'external_donor_name', 'age', 'species', 'sex',
                        'control', 'full_genotype', 'facs_population_plan', 'cre_line', 'reporter', 'injection_roi', 'injection_method', 'injection_materials', 'roi', 'patchseq_roi', 
                        'medical_conditions', 'slice_min_pos', 'slice_max_pos', 'rna_amplification_set', 'rna_amplification', 'method', 'amp_date', 'pcr_cycles', 'percent_cdna_longer_than_400bp',
                        'rna_amplification_pass_fail', 'amplified_quantity_ng', 'port_well', 'lib_method', 'lib_date', 'avg_size_bp', 'quantification2_ng', 'quantification_fmol', 'quantification2_nm',
                        'exp_cluster_density_thousands_per_mm2', 'lane_read_count', 'vendor_read_count', 'experiment_component_failed', 'library_avg_size_bp', 'library_concentration_nm', 'library_input_ng',
                        'library_prep_pass_fail', 'library_prep_set', 'library_quantification_fmol', 'library_quantification_ng', 'barcoded_cell_sample_local_name', 'barcoded_cell_sample_tag_local_name', 
                        'barcoded_cell_sample_number_of_expected_cells', 'assay', 'barcoded_cell_input_quantity_count', 'tissue_sample_local_name', 'brain_region_ontology_term_id', 'organism_ontology_term_id', 
                        'project_identifier', 'donor_id', 'donor_age', 'age_at_death_unit', 'self_reported_sex', 'organism', 'anatomical_region', 'nemo_bucket', 'ROI', 'merged_ROIs', 'n_genes_by_counts', 
                        'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes',
                        'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_genes', 'leiden_scVI', 'cluster_AIT19.3', 'Neighborhood_AIT19.3', 'Class_AIT19.3', 'Subclass_AIT19.3', 'Group_AIT19.3']

adata.obs = adata.obs[cols_to_keep_in_h5ad]
adata.write_h5ad(os.path.join(cxg_dir, f"{species}_basalganglia_{pipeline}_HMBA_scVI_LeidenCluster_tags.h5ad"))


## ----------------------------- 
## cluster-wise QC
print(adata.uns['qc_constraints'])
del adata.uns['qc_constraints']

sd.qc.add_group_level_constraint(adata, column="doublet_score", groupby = "leiden_scVI", lt=0.23, agg_func = "median")
sd.qc.add_group_level_constraint(adata, column="pct_counts_mt", groupby = "leiden_scVI", lt=2, agg_func = "median")

adata = sd.qc.apply_constraints(adata, inplace=True)
print(adata.obs.keeper_cells.value_counts())

## Remove clusters with size <100 && dominant_lib >0.9
adata = adata[~adata.obs['leiden_scVI'].isin(['370', '369', '368', '366', '365', '364', '362', '360', '359', '358', '353', '351', '350'])].copy()
print(adata.obs.keeper_cells.value_counts())

adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_scVI_LeidenCluster_filtered.h5ad")) 

#### Re-run scVI for Leiden_cluster-filtered h5ad
## -----------------------------
## Re-run scVI to remove donor effects

##
adata.layers["UMIs"] = adata.raw.X ## This is a shallow copy, no extra space is used!
sc.pp.highly_variable_genes(adata, n_top_genes=4000, layer="UMIs", subset=False, flavor="seurat_v3", batch_key="donor_id")

##
adata_hvg = adata[:, adata.var.highly_variable].copy()

#import torch
#torch.set_float32_matmul_precision("medium")

#from torch.utils.data import DataLoader
#train_dataloader = DataLoader(adata_hvg, batch_size=32, shuffle=True, num_workers=7)


## Run scVI with known confounders
scvi.model.SCVI.setup_anndata(
    adata_hvg,
    layer="UMIs",
    batch_key="donor_id", ## This holds species for Macaque (rhesus, nemestrina). Human use donor_id here?
    categorical_covariate_keys=["donor_id"],
)



model = scvi.model.SCVI(adata_hvg, 
                        dispersion="gene-batch", 
                        n_hidden=256, 
                        n_latent=64, 
                        n_layers=3)
model.train(max_epochs=200)
 
 
## Save the latent space from scVI for downstream analysis
adata.obsm["X_scVI"] = model.get_latent_representation()

del adata_hvg
del adata.layers["UMIs"]


## UMAP from scVI latent space
sc.pp.neighbors(adata, use_rep = "X_scVI")
adata = sc.tl.umap(adata, min_dist=0.3, copy=True)


adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_LeidenCluster_filtered_scVI.h5ad")) 

adata.obs = adata.obs[cols_to_keep_in_h5ad]
adata.write_h5ad(os.path.join(cxg_dir, f"{species}_basalganglia_{pipeline}_HMBA_LeidenCluster_filtered_scVI.h5ad"))
