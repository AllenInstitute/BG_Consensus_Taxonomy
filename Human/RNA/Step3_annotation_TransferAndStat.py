import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os
import sys
import sciduck as sd

##
species = "Human"
pipeline = "AIBS"

## Helpful locations which are assumed to already exist
work_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/{species}/BasalGanglia"
cxg_dir = f"/allen/programs/celltypes/workgroups/rnaseqanalysis/hct_ux3_cellxgene/anndata_080/BasalGanglia/HMBA/{species}"
os.chdir(work_dir)

adata = sc.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_iterCluster.h5ad"))
adata.obs.index = adata.obs['new_cellindex'].copy() #cellID-libID

meta = sc.read_h5ad(work_dir + "/Human_basalganglia_AIBS_BICAN_MapMyCells.h5ad") #HANN mapping results (MapMyCells)
#meta
mapmycells_cols = [col for col in meta.obs.columns if "MapMyCells" in col]
meta_filtered = meta.obs[mapmycells_cols]
meta_filtered_sorted = meta_filtered.loc[adata.obs['new_cellindex']]
adata.obs = adata.obs.join(meta_filtered_sorted, how="left")

adata.obs.drop(columns=['new_cellindex'], inplace=True)

adata.obs["cluster_doublet_score"] = adata.obs.groupby("cluster")["doublet_score"].transform("median")
adata.obs["cluster_umi_counts"] = adata.obs.groupby("cluster")["umi.counts"].transform("median")
adata.obs["cluster_number_genes"] = adata.obs.groupby("cluster")["n_genes_by_counts"].transform("median")

adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_AIT19-4_anno.h5ad"))

cols_to_keep_in_h5ad = ['umi.counts', 'library_prep',  'gene_counts_0', 'doublet_score', 'exclude', 'exclude2', 'genome', 'pipeline_version', 'gex_fraction_of_transcriptomic_reads_in_cells',
                        'gex_mean_raw_reads_per_cell', 'gex_median_umi_counts_per_cell', 'gex_median_genes_per_cell', 'gex_percent_duplicates', 'gex_q30_bases_in_umi', 'gex_q30_bases_in_barcode',
                        'gex_q30_bases_in_read_2', 'gex_q30_bases_in_sample_index_i1', 'gex_q30_bases_in_sample_index_i2', 'gex_reads_mapped_antisense_to_gene', 'gex_reads_mapped_confidently_to_exonic_regions', 
                        'gex_reads_mapped_confidently_to_genome', 'gex_reads_mapped_confidently_to_intergenic_regions', 'gex_reads_mapped_confidently_to_intronic_regions',
                        'gex_reads_mapped_confidently_to_transcriptome', 'gex_reads_mapped_to_genome', 'gex_reads_with_tso', 'gex_sequenced_read_pairs', 'gex_total_genes_detected', 
                        'keeper_median_genes', 'keeper_cells', 'percent_keeper', 'percent_doublet', 'percent_usable', 'tso_percent', 'pipeline', 'pass_fail',  'batch', 'batch_vendor_name', 
                        'tube', 'tube_internal_name', 'tube_contents_nm', 'tube_contents_nm_from_vendor', 'tube_avg_size_bp', 'tube_input_fmol', 'facs_container', 'sample_name', 'patched_cell_container', 
                        'studies', 'hemisphere_name', 'sample_quantity_count', 'sample_quantity_pg', 'donor_name', 'external_donor_name', 'age', 'species', 'sex', 'control', 'full_genotype', 'facs_population_plan',
                        'cre_line', 'reporter', 'injection_roi', 'injection_method', 'injection_materials', 'roi', 'patchseq_roi', 'medical_conditions', 'slice_min_pos', 'slice_max_pos', 'rna_amplification_set', 
                        'rna_amplification', 'method', 'amp_date', 'pcr_cycles', 'percent_cdna_longer_than_400bp','rna_amplification_pass_fail', 'amplified_quantity_ng', 'port_well', 'lib_method', 'lib_date', 
                        'avg_size_bp', 'quantification2_ng', 'quantification_fmol', 'quantification2_nm', 'exp_cluster_density_thousands_per_mm2', 'lane_read_count', 'vendor_read_count', 'experiment_component_failed', 
                        'library_avg_size_bp', 'library_concentration_nm', 'library_input_ng', 'library_prep_pass_fail', 'library_prep_set', 'library_quantification_fmol', 'library_quantification_ng', 
                        'barcoded_cell_sample_local_name', 'barcoded_cell_sample_tag_local_name', 'barcoded_cell_sample_number_of_expected_cells', 'assay', 'barcoded_cell_input_quantity_count', 'tissue_sample_local_name', 
                        'brain_region_ontology_term_id', 'organism_ontology_term_id', 'project_identifier', 'donor_id', 'donor_age', 'age_at_death_unit', 'self_reported_sex', 'organism', 'anatomical_region', 'nemo_bucket', 
                        'ROI', 'merged_ROIs', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 
                        'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_genes', 'leiden_scVI', 'cluster_AIT19.3', 'Neighborhood_AIT19.3', 'Class_AIT19.3', 'Subclass_AIT19.3', 
                        'Group_AIT19.3', 'cluster', 'cluster_doublet_score', 'cluster_umi_counts', 'cluster_number_genes']
cols_to_keep_in_h5ad = cols_to_keep_in_h5ad + list(meta_filtered_sorted.columns)
adata.obs = adata.obs[cols_to_keep_in_h5ad]
adata.write_h5ad(os.path.join(cxg_dir, f"{species}_basalganglia_{pipeline}_HMBA_AIT19-4_anno.h5ad"))
del adata

########------------------- annotation stats  -------------------########
adata = ad.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_AIT19-4_anno.h5ad"), backed="r")

## Gather columns to put into annotation table
anno = adata.obs.columns[adata.obs.columns.str.contains("|".join(["AIT117_.*_label", "AIT193_.*_label", "Silettisub_.*_label", "ABCmouse_.*_label"]), regex=True)].tolist()
numeric_anno = adata.obs.columns[adata.obs.columns.str.contains("|".join([".*_avg_correlation", ".*_bootstrapping_probability"]), regex=True)].tolist()

## Previous work
anno_prev = ['Neighborhood_AIT19.3', 'Class_AIT19.3', 'Subclass_AIT19.3', 'Group_AIT19.3', 'cluster_AIT19.3']

## Additional metadata
meta = ["donor_name", "barcoded_cell_sample_label", "sex", "age", "anatomical_region", "ROI", "merged_ROIs"]

##
anno_table = sd.anno.build_annotation_table(adata, 
                                    group_by="cluster", 
                                    categorical_annotations=anno_prev + anno + meta,
                                    numeric_annotations=["doublet_score", "n_genes_by_counts", "umi.counts", "percent_keeper", "tso_percent"] + numeric_anno, 
                                    min_percent=0.05, 
                                    annotation_alerts={
                                        "donor_name": 0.90,
                                        "barcoded_cell_sample_label": 0.90,
                                        "sex": 0.85,
                                        "merged_ROIs": 0.9,
                                        "anatomical_region": 0.8    
                                    })
anno_table.to_csv(os.path.join(work_dir, f"{species}_{pipeline}_consensus_anno_table.csv"))
