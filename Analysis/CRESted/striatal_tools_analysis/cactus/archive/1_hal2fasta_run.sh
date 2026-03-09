#!/bin/bash

BED="striatal_tools_enhancers.bed"
HAL="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/genome/hal/zoonomia/447-mammalian-2022v1_replaceAIBSexceptHuman.hal"

while IFS=$'\t' read -r chr start end length species enhancer_id; do
    # Map species names
    if [[ "$species" == "human" ]]; then
        refGenome="Homo_sapiens"
    elif [[ "$species" == "mouse" ]]; then
        refGenome="Mus_musculus"
    else
        refGenome="$species"
    fi
    ##
    OUTFILE="${enhancer_id}.maf"
    ##
    hal2maf \
      --refGenome "$refGenome" \
      --refSequence "$chr" \
      --start "$start" \
      --length "$length" \
      --noAncestors \
      "$HAL" \
      "$OUTFILE"
done < "$BED"

