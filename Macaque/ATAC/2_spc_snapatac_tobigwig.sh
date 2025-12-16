#!/bin/bash
#SBATCH --job-name spc_sa
#SBATCH --output /home/nelson.johansen/log/rheMac_bg_sa_%A_%a.out # %A is job ID, %a is array index
#SBATCH --error /home/nelson.johansen/log/rheMac_bg_sa_%A_%a.err
#SBATCH --time 24:00:00
#SBATCH --partition celltypes
#SBATCH --mem 1000gb
#SBATCH --ntasks 1

source ~/.bashrc

conda activate atac

python ~/Matthew/code/HMBA_Genomics/SpinalCord/preprocessing/atac/2_spc_snapatac_tobigwig.py $SLURM_ARRAY_TASK_ID
