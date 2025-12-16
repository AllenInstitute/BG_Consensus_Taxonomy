#!/bin/bash
#SBATCH --job-name hu_bg_atac
#SBATCH --output ./log/hu_bg_atac_sa_%A_%a.out # %A is job ID, %a is array index
#SBATCH --error ./log/hu_bg_atac_sa_%A_%a.err
#SBATCH --time 96:00:00 
#SBATCH --partition celltypes
#SBATCH --mem 1000gb
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=24

source ~/.bashrc
conda activate human_bg_atac

python ./2_snapatac_tobigwig.py
