#!/bin/bash
#SBATCH --job-name hu_bg_atac
#SBATCH --output ./log/hu_bg_atac_sa_%A_%a.out # %A is job ID, %a is array index
#SBATCH --error ./log/hu_bg_atac_sa_%A_%a.err
#SBATCH --time 6:00:00
#SBATCH --array=1-204 # 204 libs in human_BG; in .py (-1) 
#SBATCH --partition celltypes
#SBATCH --mem 48gb
#SBATCH --ntasks 1

source ~/.bash_profile
conda activate snapatac2_hpc

python /home/nelson.johansen/Analysis/HMBA_Genomics/BasalGanglia/Human/ATAC/0_snapatac_preprocess.py $SLURM_ARRAY_TASK_ID
