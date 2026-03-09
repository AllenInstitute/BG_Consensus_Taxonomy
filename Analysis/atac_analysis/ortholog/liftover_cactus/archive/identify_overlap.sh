#!/bin/bash

#SBATCH --job-name=halliftover_overlap
#SBATCH --array=0-446
#SBATCH --time 24:00:00
#SBATCH --ntasks=4
#SBATCH --mem=64GB
#SBATCH --partition=celltypes
#SBATCH --error=logs/halliftover_%A_%a.err.txt
#SBATCH --output=logs/halliftover_%A_%a.out.txt

source ~/.bash_profile
conda activate py311

python /home/nelson.johansen/Analysis/HMBA_Genomics/SpinalCord/Analysis/atac_analysis/ortholog/liftover_cactus/identify_overlap.py --species_id $SLURM_ARRAY_TASK_ID