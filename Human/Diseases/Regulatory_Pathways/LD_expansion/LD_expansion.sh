#!/bin/bash
#SBATCH --job-name=LD_expand
#SBATCH --array=0-99
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --mem=60G
#SBATCH --partition celltypes

PLINK=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Alma/tools/plink/plink

#BFILE=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Alma/MAGMA_analysis/data/MAGMA_external/european/g1000_eur   # Hg37 Panel de référence (PLINK .bed/.bim/.fam)
BFILE=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Alma/CRESTED_analysis/external/hg38_g1000/g1000_hg38_biallelic #Hg38
CHUNKDIR=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/external/Disease_related/Alzheimer/cleaned_files/LD_chunks
OUTDIR=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/external/Disease_related/Alzheimer/cleaned_files/LD_temp
TEMP_DIR=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/BasalGanglia/external/Disease_related/Alzheimer/cleaned_files/LD_temp_files
R2=0.8
KB=1000

grep -Ff $SNPS_FILE ${BFILE}.bim

mkdir -p $TEMP_DIR
mkdir -p $OUTDIR

# Selection of the file for this job
CHUNK=$(printf "%02d" ${SLURM_ARRAY_TASK_ID})
SNPS_FILE=$CHUNKDIR/snps_chunk_$CHUNK
OUTFILE=$OUTDIR/LD_chunk_$CHUNK.ld


touch $OUTFILE

while read SNP; do
    echo "[$SLURM_ARRAY_TASK_ID] Expanding $SNP"
    "$PLINK" \
      --bfile $BFILE \
      --ld-snp $SNP \
      --r2 \
      --ld-window-kb $KB \
      --ld-window 99999 \
      --ld-window-r2 $R2 \
      --allow-extra-chr \
      --out $TEMP_DIR/temp_ld_${SNP} && \
      tail -n +2 temp_ld_${SNP}.ld >> $OUTFILE

# Check if the LD file was generated and concatenate it to the output file.
    if [ -s $TEMP_DIR/temp_ld_${SNP}.ld ]; then
        echo "Appending LD data for $SNP"
        tail -n +2 $TEMP_DIR/temp_ld_${SNP}.ld >> $OUTFILE
    else
        echo "No LD data found for $SNP"
    fi

    # Remove temporary files
    rm -f $TEMP_DIR/temp_ld_${SNP}.*
done < $SNPS_FILE
