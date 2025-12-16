#!/usr/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --job-name=MAGMA_step4
#SBATCH --partition celltypes 

source "$1"

ml gcc
gwas_name=${pheno}

covar_file_spe="newest_BG_l2_specificity_matrix_sig-only" 
out_file_spe="newest_BG_l2_joint_spe_sig-only"


MAGMA_v1.10/magma --gene-results ${folder}/${file}.step2.genes.raw --gene-covar ${folder}/${covar_file_spe}_${gwas_name}.txt --model joint-pairs direction=greater --out ${folder}/${out_file_spe}_${gwas_name}
