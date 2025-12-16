#!/bin/bash
#SBATCH --job-name m1_atac_preprocess
#SBATCH --output /home/nelson.johansen/log/v1_sa_%A_%a.out # %A is replaced by job ID, %a by array index
#SBATCH --error /home/nelson.johansen/log/v1_sa_%A_%a.err
#SBATCH --time 3:00:00
#SBATCH --array=0-51 # Create an array job for each line in notebooks.txt
#SBATCH --partition celltypes
#SBATCH --mem 48gb
#SBATCH --ntasks 2

source ~/.bashrc

conda activate snapatac2_hpc

python /home/nelson.johansen/Analysis/BG/NHP/AIT117/ATAC/snapatac2_preprocess.py $SLURM_ARRAY_TASK_ID