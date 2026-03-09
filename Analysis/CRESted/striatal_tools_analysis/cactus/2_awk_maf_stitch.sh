#!/bin/bash

## For all .maf files apply .awk script to stitch together blocks that are close together
# Usage:
#   awk -f stitch_from_maf.awk input.maf > stitched.fasta

for maf in *.maf; do
    base=$(basename "$maf" .maf)
    awk -f stitch_from_maf.awk "$maf" > "${base}_stitched.fasta"
    mafft --globalpair --maxiterate 1000 --thread -1 "${base}_stitched.fasta" > "${base}_stitched_and_aligned.fasta"
done