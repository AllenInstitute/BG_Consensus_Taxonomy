#!/bin/bash
#SBATCH --job-name=Full_MAGMA_pipeline
#SBATCH --partition celltypes 
#SBATCH --time=1:00:00
#SBATCH --mem=20G

source "$1"

Rscript create_top_results_matrix_step3.R "$filename_marg_result_R" "$filename_output_R" "$genome_built" "$base_dir_BG_dataset"
