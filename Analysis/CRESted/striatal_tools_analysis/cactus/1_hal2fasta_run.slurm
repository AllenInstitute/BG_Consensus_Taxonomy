#!/bin/bash
#SBATCH --job-name=hal2maf_array
#SBATCH --output=logs/hal2maf_%A_%a.out
#SBATCH --error=logs/hal2maf_%A_%a.err
#SBATCH --array=1-34        # <-- set this to the number of lines in your BED
#SBATCH --partition=celltypes
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

source ~/.bash_profile
conda activate hal

BED="striatal_tools_enhancers.bed"
HAL="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/genome/hal/zoonomia/447-mammalian-2022v1_replaceAIBSexceptHuman.hal"

# Get the line corresponding to this array index
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BED")
read -r chr start end length species enhancer_id <<< "$LINE"

# Map species to proper names for hal
if [[ "$species" == "human" ]]; then
    refGenome="Homo_sapiens"
elif [[ "$species" == "mouse" ]]; then
    refGenome="Mus_musculus"
else
    refGenome="$species"
fi

# Output file
OUTFILE="${enhancer_id}.maf"

# Run hal2maf
hal2maf \
  --refGenome "$refGenome" \
  --refSequence "$chr" \
  --start "$start" \
  --length "$length" \
  --noAncestors \
  "$HAL" \
  "$OUTFILE"
