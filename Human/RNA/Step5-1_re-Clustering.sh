#!/bin/bash
#SBATCH --job-name=hicat_mpi
#SBATCH --partition=celltypes
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00
#SBATCH --mem=500G
#SBATCH --mail-type=END,FAIL

#  please refer to https://github.com/dyuan1111/hicatMPI/tree/main for accelerated verion

adata_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalGanglia/Human_basalganglia_AIBS_HMBA_AIT19-4_anno_filtered_scVI.h5ad" # norm or counts in X is fine.
out_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/HMBA/Aim1_Atlases/Human/BasalGanglia/" # output directory


clust_kwargs="{
    'means_vars_kwargs': {
        'low_thresh': 0.6931472, # lowest value required for a gene to pass filtering. set to 1 originally, change to 0.6931472 to match to the bigcat default
        'min_cells': 4 # minimum number of cells expressed required for a gene to pass filtering
    },
    'highly_variable_kwargs': {
        'max_genes': 4000 # originally 3000, change to 4000 to match to the bigcat default
    },
    'pca_kwargs': {
        'cell_select': 30000, # originally 500000 cells
        'n_comps': 50,
        'svd_solver': 'randomized'
    },
    'filter_pcs_kwargs': {
        'known_components': None,
        'similarity_threshold': 0.7,
        'method': 'zscore', # or elbow
        'zth': 2,
        'max_pcs': None,
    },
    # if not using known_modes, set filter_known_modes_kwargs to None or an empty dict. Only applies for PCA
    'filter_known_modes_kwargs': {
        'known_modes': 'log2ngene', 
        'similarity_threshold': 0.7
    },
    ## !!NEW!! Original method: "PCA", allows the user to select any obsm latent space such as "X_scVI" for leiden clustering.
    'latent_kwargs': {
        'latent_component': 'X_scVI' # None (default to run PCA) or obsm key such as "X_scVI"
    },
    'cluster_louvain_kwargs': {
        'k': 15, # number of nn, originally 150, change to 15 to match to the bigcat default
        'nn_measure': 'euclidean',
        'knn_method': 'annoy',
        'louvain_method': 'taynaud', #'vtraag',
        'weighting_method': 'jaccard',
        'n_jobs': 30, # cpus, originally 8
        'resolution': 1.0 # resolution of louvain for taynaud method
    },
    'merge_clusters_kwargs': {
        'thresholds': {
            'q1_thresh': 0.5,
            'q2_thresh': None,
            'cluster_size_thresh': 10, ## originally uses 50, change to 10 to match to the bigcat default
            'qdiff_thresh': 0.7, 
            'padj_thresh': 0.05, 
            'lfc_thresh': 0.6931472, # log2 fold change threshold for DE genes
            'score_thresh': 100, # originally uses 200, change to 100 to match to the bigcat default
            'low_thresh': 0.6931472, # originally uses 1 # applied to log2(cpm+1) to determine if a gene is expressed or not, change to 0.6931472 to match to the bigcat default
            'min_genes': 5
        },
        'k': 4, # number of nn for de merge, originaly 2, change to 4 to match to the bigcat default
        'de_method': 'ebayes'
    }
}"



manager_script="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/custom_packages/mpi_tc/hicatMPI/iterative_clustering_mpi_manager.py"
worker_script="/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/custom_packages/mpi_tc/hicatMPI/iterative_clustering_mpi_worker.py"

source /allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/miniconda3/etc/profile.d/conda.sh
conda activate /allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/miniconda3/envs/rapids_singlecell
export PYTHONPATH=/allen/programs/celltypes/workgroups/rnaseqanalysis/dyuan/tool/transcriptomic_clustering:$PYTHONPATH

module load mpi/mpich-3.2-x86_64 #only in old Slurm
#module load mpi/mpich-x86_64

# Navigate to the output directory, where the log and out files will be saved
# Check if the out directory exists. if not, create it
if [ ! -d "$out_dir" ]; then
    # If the directory doesn't exist, create it
    mkdir -p "$out_dir"
    echo "Directory '$out_dir' created."
fi
cd "$out_dir"

time mpiexec -n 1 sh -c "python \"$manager_script\" \"$adata_dir\" \"$latent_dir\" \"$out_dir\" \"$clust_kwargs\"> manager_output.log 2> manager_error.log" : -n 5 sh -c "python \"$worker_script\" > worker_output.log 2> worker_error.log" # somehow worker_output.log is trucated, but its ok.
