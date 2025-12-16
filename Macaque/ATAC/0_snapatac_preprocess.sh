#!/bin/bash
#SBATCH --job-name rheMac_bg
#SBATCH --output /home/nelson.johansen/log/rheMac_bg_sa_%A_%a.out # %A is job ID, %a is array index
#SBATCH --error /home/nelson.johansen/log/rheMac_bg_sa_%A_%a.err
#SBATCH --time 10:00:00
#SBATCH --array=1-90 
#SBATCH --partition celltypes
#SBATCH --mem 64gb
#SBATCH --ntasks 4

source ~/.bash_profile

conda activate snapatac2_hpc

python /home/nelson.johansen/Analysis/HMBA_Genomics/BasalGanglia/Macaque/ATAC/0_spc_snapatac_preprocess.py $SLURM_ARRAY_TASK_ID
