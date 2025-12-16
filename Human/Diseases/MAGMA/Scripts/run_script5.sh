#!/bin/bash
#SBATCH --job-name=Full_MAGMA_pipeline
#SBATCH --partition celltypes 
#SBATCH --time=1:00:00
#SBATCH --mem=20G

source "$1"

Rscript forward_selection_condition_results_step5.R "$filename_cond" "$filename_marg" "$pheno" "$ancestry"
