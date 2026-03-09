#!/bin/bash

#SBATCH --job-name=halliftover
#SBATCH --time 24:00:00
#SBATCH --ntasks=4
#SBATCH --mem=64GB
#SBATCH --partition=celltypes
#SBATCH --error=logs/halliftover_%A_%a.err.txt
#SBATCH --output=logs/halliftover_%A_%a.out.txt
#SBATCH --array=1-892

MAX_FRAC=2
MIN_LEN=50
PROTECT_DIST=10
SPECIES="macaque"
SPECIES_SCI="Macaca_mulatta"
HAL_DIR="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/${SPECIES}"
UNIQUEBED="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/${SPECIES}/Consensus_peaks/merged_peaks_with_names.bed"
TARGETS="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/liftOver/species.txt"

# Get array of targets, either from comma-delimited string, or by reading from a file.
if [ ! -f "$TARGETS" ]; then
    # interpret TARGETS as a comma-delimited string and read elements into array
    TARGETS=( $(echo $TARGETS | tr "," "\n") )
else
    # TARGETS is a file, read lines from it and store as array
    echo "Reading in targets line-by-line from $TARGETS"
    unset -v lines
    while IFS= read -r; do
        lines+=("$REPLY")
    done <$TARGETS
    [[ $REPLY ]] && lines+=("$REPLY")
    TARGETS=("${lines[@]}")
fi

## array job
TARGET=${TARGETS[$(($SLURM_ARRAY_TASK_ID - 1))]}
echo "Processing target: $TARGET"

HALLIFTEDTFILE="merged_peaks_with_names.${SPECIES_SCI}To${TARGET}.halLiftover.tFile.bed"
HALLIFTEDSFILE="merged_peaks_with_names.${SPECIES_SCI}To${TARGET}.halLiftover.sFile.bed"
zcat "${HAL_DIR}/$HALLIFTEDTFILE.gz" > "${HAL_DIR}/$HALLIFTEDTFILE"
zcat "${HAL_DIR}/$HALLIFTEDSFILE.gz" > "${HAL_DIR}/$HALLIFTEDSFILE"

OUTFILE="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/SpinalCord/manuscript/ATAC/hal/${SPECIES}/merged_peaks_with_names.${SPECIES_SCI}To${TARGET}.HALPER_${MAX_FRAC}.narrowPeak"

##
args="-max_frac $MAX_FRAC -min_len $MIN_LEN -protect_dist $PROTECT_DIST \
    -tFile ${HAL_DIR}/$HALLIFTEDTFILE -sFile ${HAL_DIR}/$HALLIFTEDSFILE \
    -qFile $UNIQUEBED -oFile $OUTFILE"
##
python -m orthologFind $args

## Compress
gzip $OUTFILE

## Remove non gz files
rm -f "${HAL_DIR}/$HALLIFTEDTFILE" "${HAL_DIR}/$HALLIFTEDSFILE"
