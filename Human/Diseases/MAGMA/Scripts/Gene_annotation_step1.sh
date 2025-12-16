#!/usr/bin/bash
#SBATCH --job-name=MAGMA_step1
#SBATCH --partition celltypes 
#SBATCH --time=2:00:00
#SBATCH --mem=16G

source "$1"

echo "GWAS file : ${folder}/${file}"

ml gcc

snp_col="ID"
p_col="PVAL"
ncol="NEFF"

MAGMA_v1.10/magma --annotate window=35,10 \
--snp-loc ${folder}/snp_loc_${file} \
--gene-loc External_data/${genome_built}.3.gene.loc \
--out ${folder}/${file}.step1

MAGMA_v1.10/magma --bfile External_data/${ancestry}/g1000_${ancestry_MAGMA} \ 
--pval ${folder}/${file} use=${snp_col},${p_col} ncol=${ncol} \
--gene-annot ${folder}/${file}.step1.genes.annot \
--out ${folder}/${file}.step2
