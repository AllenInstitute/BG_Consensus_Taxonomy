import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scvi
import os
import sys
from datetime import datetime
import importlib
import math
import json
import transcriptomic_clustering as tc
from transcriptomic_clustering.iterative_clustering import (
    build_cluster_dict, iter_clust, OnestepKwargs, summarize_final_clusters
)
import pickle

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

adata = sc.read_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_LeidenCluster_filtered_scVI.h5ad"))


########---------------- iterative clustering ----------------########
##transcriptomic clustering following scrattch.hicat

print("iter_Clustering Start Time : " + datetime.now().strftime("%H:%M:%S"))

## Assign keyword arguments (kwargs). Any unassigned args will be set to their respective function defaults.
def setup_transcriptomic_clustering(): 
    means_vars_kwargs = {
        'low_thresh': 0.6931472,
        'min_cells': 4
    }
    highly_variable_kwargs = {
        'max_genes': 4000
    }
    pca_kwargs = {
        'cell_select': 30000,
        'n_comps': 50,
        'svd_solver': 'randomized'
    }
    filter_pcs_kwargs = {
        'known_components': None,
        'similarity_threshold': 0.7,
        'method': 'zscore', # or elbow
        'zth': 2,
        'max_pcs': None,
    }
    ## Leave empty if you don't want to use known_modes
    filter_known_modes_kwargs = {
        # 'known_modes': known_modes_df,
        'similarity_threshold': 0.7
    }
    ## !!NEW!! Original method: "PCA", allows the user to select any obsm latent space such as "X_scVI" for leiden clustering.
    latent_kwargs = {
        'latent_component': "X_scVI"
    }
    cluster_louvain_kwargs = {
        'k': 15, # number of nn
        'nn_measure': 'euclidean',
        'knn_method': 'annoy',
        'louvain_method': 'taynaud', #'vtraag',
        'weighting_method': 'jaccard',
        'n_jobs': 36, # cpus
        'resolution': 1.0, # resolution of louvain for taynaud method
    }
    merge_clusters_kwargs = {
        'thresholds': {
            'q1_thresh': 0.5,
            'q2_thresh': None,
            'cluster_size_thresh': 10,
            'qdiff_thresh': 0.7,
            'padj_thresh': 0.05,
            'lfc_thresh': 0.6931472,
            'score_thresh': 100,
            'low_thresh': 0.6931472,
            'min_genes': 5
        },
        'k': 4, # number of nn for de merge 
        'de_method': 'ebayes'
    }
    onestep_kwargs = OnestepKwargs(
        means_vars_kwargs = means_vars_kwargs,
        highly_variable_kwargs = highly_variable_kwargs,
        pca_kwargs = pca_kwargs,
        filter_pcs_kwargs = filter_pcs_kwargs,
        filter_known_modes_kwargs = filter_known_modes_kwargs,
        latent_kwargs = latent_kwargs,
        cluster_louvain_kwargs = cluster_louvain_kwargs,
        merge_clusters_kwargs = merge_clusters_kwargs
    )
    return onestep_kwargs

##
onestep_kwargs = setup_transcriptomic_clustering()


###  Run clustering
clusters, markers = iter_clust(
    adata,
    onestep_kwargs=onestep_kwargs,
    min_samples=4,
    random_seed=123,
    tmp_dir= work_dir + "/tmp2"
) 


# save the results in .pkl format, which will be used for final merging
with open(work_dir + 'hmba_human_bg_clustering_results.pkl', 'wb') as f:
    pickle.dump(clusters, f)

with open(work_dir + 'hmba_human_bg_markers.pkl', 'wb') as f:
    pickle.dump(markers, f)


cluster_dict = build_cluster_dict(clusters)
with open(work_dir + '/hmba_human_bg_iter_clusters_beforeMerge.json', 'w') as f:
    json.dump(cluster_dict, f)

adata.obs["clusterBeforeMerge"] = ""
for cluster in cluster_dict.keys():
    adata.obs.clusterBeforeMerge[cluster_dict[cluster]] = cluster

adata.obs.clusterBeforeMerge = adata.obs.clusterBeforeMerge.astype('category')


## write the h5ad with adding a new column, called cluster. that is the ground truth for downstream analysis. For example, mapping annotation transfer.
adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_iterCluster_beforeMerge.h5ad")) 

print("iter_Clustering End Time : " + datetime.now().strftime("%H:%M:%S"))



#########----------------------------- Final Merge of iterative clustering -----------------------------########

from transcriptomic_clustering.final_merging import final_merge, FinalMergeKwargs

def setup_merging(): 
    merge_clusters_kwargs = {
        'thresholds': {
            'q1_thresh': 0.5,
            'q2_thresh': None,
            'cluster_size_thresh': 10, 
            'qdiff_thresh': 0.7, 
            'padj_thresh': 0.05, 
            'lfc_thresh': 0.6931472, 
            'score_thresh': 100, 
            'low_thresh': 0.6931472, 
            'min_genes': 5
        },
        'k': 4,
        'de_method': 'ebayes',
        'n_markers': None,
    }
    latent_kwargs = {  
        "latent_component": "X_scVI"
    }

    merge_kwargs = FinalMergeKwargs(
        merge_clusters_kwargs=merge_clusters_kwargs,
        latent_kwargs=latent_kwargs  
    )
    return merge_kwargs

merge_kwargs = setup_merging()


clusters_as_lists = [array.tolist() for array in clusters]

# Run the final merging
clusters_after_merging, markers = final_merge(
    adata, 
    clusters_as_lists, 
    n_samples_per_clust=20, 
    random_seed=2025, 
    n_jobs = 36, # modify this to the number of cores you want to use
    return_markers_df = False, # return the pair-wise DE results for each cluster pair. If False (default), only return a set of markers (top 20 of up and down regulated genes in each pair comparison)
    final_merge_kwargs=merge_kwargs
)


with open(os.path.join(work_dir, "hmba_human_bg_clustering_results_after_merging.pkl"), 'wb') as f:
        pickle.dump(clusters_after_merging, f)

cluster_after_merging_dict = build_cluster_dict(clusters_after_merging)
with open(work_dir + '/hmba_human_bg_iter_clusters_afterMerge.json', 'w') as f:
    json.dump(cluster_after_merging_dict, f)

adata.obs["cluster"] = ""
for cluster in cluster_after_merging_dict.keys():
    adata.obs.cluster[cluster_after_merging_dict[cluster]] = cluster

adata.obs.cluster = adata.obs.cluster.astype('category')


## write the h5ad with adding a new column, called cluster. that is the ground truth for downstream analysis. For example, mapping annotation transfer.
adata.write_h5ad(os.path.join(work_dir, f"{species}_basalganglia_{pipeline}_HMBA_iterCluster.h5ad")) 


cols_to_keep_in_h5ad = ['umi.counts', 'library_prep',  'gene_counts_0', 'doublet_score', 'exclude', 'exclude2', 'genome', 'pipeline_version', 'gex_fraction_of_transcriptomic_reads_in_cells',
                        'gex_mean_raw_reads_per_cell', 'gex_median_umi_counts_per_cell', 'gex_median_genes_per_cell', 'gex_percent_duplicates', 'gex_q30_bases_in_umi', 'gex_q30_bases_in_barcode',
                        'gex_q30_bases_in_read_2', 'gex_q30_bases_in_sample_index_i1', 'gex_q30_bases_in_sample_index_i2', 'gex_reads_mapped_antisense_to_gene', 'gex_reads_mapped_confidently_to_exonic_regions',
                        'gex_reads_mapped_confidently_to_genome', 'gex_reads_mapped_confidently_to_intergenic_regions', 'gex_reads_mapped_confidently_to_intronic_regions', 'gex_reads_mapped_confidently_to_transcriptome',
                        'gex_reads_mapped_to_genome', 'gex_reads_with_tso', 'gex_sequenced_read_pairs', 'gex_total_genes_detected', 'keeper_median_genes', 'keeper_cells', 'percent_keeper', 'percent_doublet',
                        'percent_usable', 'tso_percent', 'pipeline', 'pass_fail',  'batch', 'batch_vendor_name', 'tube', 'tube_internal_name', 'tube_contents_nm', 'tube_contents_nm_from_vendor', 'tube_avg_size_bp', 
                        'tube_input_fmol', 'facs_container', 'sample_name', 'patched_cell_container', 'studies', 'hemisphere_name', 'sample_quantity_count', 'sample_quantity_pg', 'donor_name', 'external_donor_name',
                        'age', 'species', 'sex', 'control', 'full_genotype', 'facs_population_plan', 'cre_line', 'reporter', 'injection_roi', 'injection_method', 'injection_materials', 'roi', 'patchseq_roi', 
                        'medical_conditions', 'slice_min_pos', 'slice_max_pos', 'rna_amplification_set', 'rna_amplification', 'method', 'amp_date', 'pcr_cycles', 'percent_cdna_longer_than_400bp','rna_amplification_pass_fail',
                        'amplified_quantity_ng', 'port_well', 'lib_method', 'lib_date', 'avg_size_bp', 'quantification2_ng', 'quantification_fmol', 'quantification2_nm', 'exp_cluster_density_thousands_per_mm2', 
                        'lane_read_count', 'vendor_read_count', 'experiment_component_failed', 'library_avg_size_bp', 'library_concentration_nm', 'library_input_ng', 'library_prep_pass_fail', 'library_prep_set',
                        'library_quantification_fmol', 'library_quantification_ng', 'barcoded_cell_sample_local_name', 'barcoded_cell_sample_tag_local_name', 'barcoded_cell_sample_number_of_expected_cells', 'assay',
                        'barcoded_cell_input_quantity_count', 'tissue_sample_local_name', 'brain_region_ontology_term_id', 'organism_ontology_term_id', 'project_identifier', 'donor_id', 'donor_age', 'age_at_death_unit',
                        'self_reported_sex', 'organism', 'anatomical_region', 'nemo_bucket', 'ROI', 'merged_ROIs', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts',
                        'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 
                        'total_genes', 'leiden_scVI', 'cluster_AIT19.3', 'Neighborhood_AIT19.3', 'Class_AIT19.3', 'Subclass_AIT19.3', 'Group_AIT19.3', 'clusterBeforeMerge', 'cluster']
adata.obs = adata.obs[cols_to_keep_in_h5ad]
adata.write_h5ad(os.path.join(cxg_dir, f"{species}_basalganglia_{pipeline}_HMBA_iterCluster.h5ad"))

print("iter_Clustering End Time : " + datetime.now().strftime("%H:%M:%S"))


