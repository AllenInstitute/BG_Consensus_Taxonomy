#!/usr/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --job-name=MAGMA_step2
#SBATCH --partition celltypes 

source "$1"

ml gcc
outfile=newest_BG_l2_conti-spe_${pheno}

covar_file_conti="${base_dir_BG_dataset}${genome_built}/conti_specificity_matrix.txt"
MAGMA_v1.10/magma --gene-results ${folder}/${file}.step2.genes.raw --model direction=greater --gene-covar ${covar_file_conti} --out ${folder}/${outfile}
